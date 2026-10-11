/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"

#include <fstream>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::constitutive::multifluid;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

static void writeFile( string const & fileName, string const & content )
{
  std::ofstream os( fileName );
  ASSERT_TRUE( os.is_open() );
  os << content;
}

static void writeKValueTables()
{
  writeFile( "kv_pressure.txt", "1.0e5\n5.0e7\n" );
  writeFile( "kv_temperature.txt", "290.15\n400.15\n" );
  auto const writeK = []( string const & name, real64 const k0 )
  {
    writeFile( name,
               GEOS_FMT( "{0}\n{1}\n{2}\n{3}\n", k0, k0 * 1.02, k0 * 0.98, k0 * 1.01 ) );
  };
  writeK( "kv_CO2_vapour.txt", 3.0 );
  writeK( "kv_C1C2_vapour.txt", 8.0 );
  writeK( "kv_H2O_vapour.txt", 0.2 );
  writeK( "kv_CO2_aqueous.txt", 0.8 );
  writeK( "kv_C1C2_aqueous.txt", 0.01 );
  writeK( "kv_H2O_aqueous.txt", 60.0 );
}

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers>
      <CompositionalMultiphaseFVM name="compflow"
                                  logLevel="1"
                                  discretization="fluidTPFA"
                                  temperature="368.15"
                                  useMass="1"
                                  isThermal="1"
                                  initialDt="1000"
                                  maxCompFractionChange="0.5"
                                  targetRegions="{ region }">
        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="100"
                                   lineSearchAction="None"
                                   maxTimeStepCuts="5" />
        <LinearSolverParameters directParallel="0" />
      </CompositionalMultiphaseFVM>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{ C3D8 }"
                    xCoords="{ 0, 20 }"
                    yCoords="{ 0, 1 }"
                    zCoords="{ 0, 1 }"
                    nx="{ 5 }"
                    ny="{ 1 }"
                    nz="{ 1 }"
                    cellBlockNames="{ cb }" />
    </Mesh>
    <Geometry>
      <Box name="sink"
           xMin="{ -0.01, -0.01, -0.01 }"
           xMax="{ 4.01, 1.01, 1.01 }" />
      <Box name="source"
           xMin="{ -0.01, -0.01, -0.01 }"
           xMax="{ 4.01, 1.01, 1.01 }" />
    </Geometry>
    <Events maxTime="1000">
      <PeriodicEvent name="solverApplications"
                     maxEventDt="1000"
                     target="/Solvers/compflow" />
    </Events>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA" />
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region"
                         cellBlocks="{ * }"
                         materialList="{ fluid, rock, relperm, thermalCond }" />
    </ElementRegions>
    <Constitutive>
      <CompressibleSolidConstantPermeability name="rock"
                                             solidModelName="nullSolid"
                                             porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"
                                             solidInternalEnergyModelName="rockInternalEnergy"/>
      <NullModel name="nullSolid" />
      <PressurePorosity name="rockPorosity"
                        defaultReferencePorosity="0.2"
                        referencePressure="0.0"
                        compressibility="1.0e-9" />
      <SolidInternalEnergy name="rockInternalEnergy"
                           referenceVolumetricHeatCapacity="1.95e6"
                           referenceTemperature="368.15"
                           referenceInternalEnergy="0" />
      <ConstantPermeability name="rockPerm"
                            permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }" />
      <CompositionalThreePhaseKValueThermalFluidLohrenzBrayClark
        name="fluid"
        phaseNames="{ oil, gas, water }"
        equationsOfState="{ PengRobinson, PengRobinson, PengRobinson }"
        componentNames="{ CO2, C1C2, H2O }"
        componentMolarWeight="{ 44.01e-3, 22.0e-3, 18.015e-3 }"
        componentCriticalPressure="{ 73.77e5, 45.99e5, 220.64e5 }"
        componentCriticalTemperature="{ 304.13, 190.56, 647.10 }"
        componentAcentricFactor="{ 0.2239, 0.0114, 0.3443 }"
        kValueTables="{ kv_CO2_vapour, kv_C1C2_vapour, kv_H2O_vapour, kv_CO2_aqueous, kv_C1C2_aqueous, kv_H2O_aqueous }"
        componentHeatCapacityCoefficients="{ { 35.0, 0.0, 0.0, 0.0, 0.0 },
                                             { 35.0, 0.0, 0.0, 0.0, 0.0 },
                                             { 35.0, 0.0, 0.0, 0.0, 0.0 } }" />
      <BrooksCoreyRelativePermeability name="relperm"
                                       phaseNames="{ oil, gas, water }"
                                       phaseMinVolumeFraction="{ 0.0, 0.0, 0.0 }"
                                       phaseRelPermExponent="{ 1.5, 1.5, 1.5 }"
                                       phaseRelPermMaxValue="{ 0.9, 0.9, 0.9 }" />
      <MultiPhaseConstantThermalConductivity name="thermalCond"
                                             phaseNames="{ oil, gas, water }"
                                             thermalConductivityComponents="{ 0.6, 0.6, 0.6 }" />
    </Constitutive>
    <Functions>
      <TableFunction name="kv_CO2_vapour"
                     inputVarNames="{ pressure, temperature }"
                     interpolation="linear"
                     coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }"
                     voxelFile="kv_CO2_vapour.txt"/>
      <TableFunction name="kv_C1C2_vapour"
                     inputVarNames="{ pressure, temperature }"
                     interpolation="linear"
                     coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }"
                     voxelFile="kv_C1C2_vapour.txt"/>
      <TableFunction name="kv_H2O_vapour"
                     inputVarNames="{ pressure, temperature }"
                     interpolation="linear"
                     coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }"
                     voxelFile="kv_H2O_vapour.txt"/>
      <TableFunction name="kv_CO2_aqueous"
                     inputVarNames="{ pressure, temperature }"
                     interpolation="linear"
                     coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }"
                     voxelFile="kv_CO2_aqueous.txt"/>
      <TableFunction name="kv_C1C2_aqueous"
                     inputVarNames="{ pressure, temperature }"
                     interpolation="linear"
                     coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }"
                     voxelFile="kv_C1C2_aqueous.txt"/>
      <TableFunction name="kv_H2O_aqueous"
                     inputVarNames="{ pressure, temperature }"
                     interpolation="linear"
                     coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }"
                     voxelFile="kv_H2O_aqueous.txt"/>
    </Functions>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="9e6" />
      <FieldSpecification name="initialTemperature"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="temperature"
                          scale="368.15" />
      <FieldSpecification name="initialComposition_co2"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="0"
                          scale="0.1" />
      <FieldSpecification name="initialComposition_c1c2"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="1"
                          scale="0.2" />
      <FieldSpecification name="initialComposition_water"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="2"
                          scale="0.7" />
      <FieldSpecification name="sinkPressure"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="7e6" />
      <FieldSpecification name="sinkTemperature"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="temperature"
                          scale="368.15" />
      <FieldSpecification name="sinkTermComposition_co2"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="0"
                          scale="0.1" />
      <FieldSpecification name="sinkTermComposition_c1c2"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="1"
                          scale="0.2" />
      <FieldSpecification name="sinkTermComposition_water"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="2"
                          scale="0.7" />
      <FieldSpecification name="sourcePressure"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1.45e7" />
      <FieldSpecification name="sourceTemperature"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="temperature"
                          scale="300.15" />
      <FieldSpecification name="sourceTermComposition_co2"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="0"
                          scale="0.8" />
      <FieldSpecification name="sourceTermComposition_c1c2"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="1"
                          scale="0.1" />
      <FieldSpecification name="sourceTermComposition_water"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="2"
                          scale="0.1" />
    </FieldSpecifications>
  </Problem>
  )xml";

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseFVM & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );

  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( hostMemorySpace, false );

  array1d< real64 > residualOrig( residual );

  jacobian.move( hostMemorySpace );
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  fillCellCenteredNumericalJacobian( solver,
                                     domain,
                                     true,
                                     perturbParameter,
                                     residual.toView(),
                                     residualOrig.toView(),
                                     jacobian.toView(),
                                     jacobianFD.toView(),
                                     assembleFunction );

  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol, 1e-6 );
}

class ThermalKValueCompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  ThermalKValueCompositionalMultiphaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    writeKValueTables();
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "compflow" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

    solver->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  CompositionalMultiphaseFVM * solver;
};

real64 constexpr ThermalKValueCompositionalMultiphaseFlowTest::time;
real64 constexpr ThermalKValueCompositionalMultiphaseFlowTest::dt;
real64 constexpr ThermalKValueCompositionalMultiphaseFlowTest::eps;

TEST_F( ThermalKValueCompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-4;

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testCompositionNumericalDerivatives( *solver, domain, perturb, tol );
}

TEST_F( ThermalKValueCompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1;

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
