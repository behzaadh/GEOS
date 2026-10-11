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

#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"

#include <fstream>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
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
<?xml version="1.0" ?>
<Problem>
  <Solvers>
    <CompositionalMultiphaseReservoir
      name="reservoirSystem"
      flowSolverName="compflow"
      wellSolverName="compositionalMultiphaseWell"
      logLevel="1"
      initialDt="1e4"
      targetRegions="{ region, injwell }">
      <NonlinearSolverParameters
        newtonTol="1.0e-5"
        lineSearchAction="None"
        newtonMaxIter="40"/>
      <LinearSolverParameters
        directParallel="0"/>
    </CompositionalMultiphaseReservoir>
    <CompositionalMultiphaseFVM
      name="compflow"
      logLevel="1"
      discretization="fluidTPFA"
      temperature="368.15"
      useMass="0"
      isThermal="1"
      initialDt="1e4"
      maxCompFractionChange="0.5"
      targetRegions="{ region }">
    </CompositionalMultiphaseFVM>
    <WellManager
      name="compositionalMultiphaseWell"
      targetRegions="{ injwell }"
      isThermal="1"
      logLevel="1"
      initialDt="1e4"
      useMass="0">
      <CompositionalMultiphaseWell
        name="WC_CO2_INJ"
        type="injector"
        useSurfaceConditions="1"
        surfacePressure="101325"
        surfaceTemperature="288.706"
        enableCrossflow="0"
        control="totalVolRate">
        <MaximumBHPConstraint
          name="maxbhp"
          targetBHP="1e8"
          referenceElevation="-0.01"/>
        <InjectionVolumeRateConstraint
          name="maxvolrateinj"
          volumeRate="0.001"
          injectionStream="{ 0.9, 0.0, 0.1 }"
          injectionTemperature="306"/>
      </CompositionalMultiphaseWell>
    </WellManager>
  </Solvers>
  <Mesh>
    <InternalMesh
      name="mesh1"
      elementTypes="{ C3D8 }"
      xCoords="{ 0, 20 }"
      yCoords="{ 0, 20 }"
      zCoords="{ 0, 1 }"
      nx="{ 2 }"
      ny="{ 2 }"
      nz="{ 1 }"
      cellBlockNames="{ cb }">
      <InternalWell
        name="inj11"
        wellRegionName="injwell"
        wellControlsName="WC_CO2_INJ"
        logLevel="1"
        polylineNodeCoords="{{5.0,5.0,1.01},{5.0,5.0,0.01}}"
        polylineSegmentConn="{{0,1}}"
        radius="0.1"
        numElementsPerSegment="1">
        <Perforation
          name="inj1_perf1"
          distanceFromHead="0.5"/>
      </InternalWell>
    </InternalMesh>
  </Mesh>
  <Events maxTime="1.0e4">
    <PeriodicEvent
      name="solverApplications"
      maxEventDt="1.0e4"
      target="/Solvers/reservoirSystem"/>
  </Events>
  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="fluidTPFA"/>
    </FiniteVolume>
  </NumericalMethods>
  <ElementRegions>
    <CellElementRegion
      name="region"
      cellBlocks="{ cb }"
      materialList="{ fluid, rock, relperm, thermalCond }"/>
    <WellElementRegion
      name="injwell"
      materialList="{ fluid, relperm }"/>
  </ElementRegions>
  <Constitutive>
    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"
      solidInternalEnergyModelName="rockInternalEnergy"/>
    <NullModel name="nullSolid"/>
    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="0.2"
      referencePressure="0.0"
      compressibility="1.0e-9"/>
    <SolidInternalEnergy
      name="rockInternalEnergy"
      referenceVolumetricHeatCapacity="1.95e6"
      referenceTemperature="368.15"
      referenceInternalEnergy="0"/>
    <ConstantPermeability
      name="rockPerm"
      permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>
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
                                           { 35.0, 0.0, 0.0, 0.0, 0.0 } }"/>
    <BrooksCoreyRelativePermeability
      name="relperm"
      phaseNames="{ oil, gas, water }"
      phaseMinVolumeFraction="{ 0.0, 0.0, 0.0 }"
      phaseRelPermExponent="{ 1.5, 1.5, 1.5 }"
      phaseRelPermMaxValue="{ 0.9, 0.9, 0.9 }"/>
    <MultiPhaseConstantThermalConductivity
      name="thermalCond"
      phaseNames="{ oil, gas, water }"
      thermalConductivityComponents="{ 0.6, 0.6, 0.6 }"/>
  </Constitutive>
  <Functions>
    <TableFunction name="kv_CO2_vapour" inputVarNames="{ pressure, temperature }" interpolation="linear"
                   coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }" voxelFile="kv_CO2_vapour.txt"/>
    <TableFunction name="kv_C1C2_vapour" inputVarNames="{ pressure, temperature }" interpolation="linear"
                   coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }" voxelFile="kv_C1C2_vapour.txt"/>
    <TableFunction name="kv_H2O_vapour" inputVarNames="{ pressure, temperature }" interpolation="linear"
                   coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }" voxelFile="kv_H2O_vapour.txt"/>
    <TableFunction name="kv_CO2_aqueous" inputVarNames="{ pressure, temperature }" interpolation="linear"
                   coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }" voxelFile="kv_CO2_aqueous.txt"/>
    <TableFunction name="kv_C1C2_aqueous" inputVarNames="{ pressure, temperature }" interpolation="linear"
                   coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }" voxelFile="kv_C1C2_aqueous.txt"/>
    <TableFunction name="kv_H2O_aqueous" inputVarNames="{ pressure, temperature }" interpolation="linear"
                   coordinateFiles="{ kv_pressure.txt, kv_temperature.txt }" voxelFile="kv_H2O_aqueous.txt"/>
  </Functions>
  <FieldSpecifications>
    <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }"
                        objectPath="ElementRegions/region/cb" fieldName="pressure" scale="9e6"/>
    <FieldSpecification name="initialTemperature" initialCondition="1" setNames="{ all }"
                        objectPath="ElementRegions/region/cb" fieldName="temperature" scale="368.15"/>
    <FieldSpecification name="initialComposition_co2" initialCondition="1" setNames="{ all }"
                        objectPath="ElementRegions/region/cb" fieldName="globalCompFraction" component="0" scale="0.1"/>
    <FieldSpecification name="initialComposition_c1c2" initialCondition="1" setNames="{ all }"
                        objectPath="ElementRegions/region/cb" fieldName="globalCompFraction" component="1" scale="0.2"/>
    <FieldSpecification name="initialComposition_water" initialCondition="1" setNames="{ all }"
                        objectPath="ElementRegions/region/cb" fieldName="globalCompFraction" component="2" scale="0.7"/>
  </FieldSpecifications>
</Problem>
)xml";

class ThermalKValueReservoirWellTest : public ::testing::Test
{
public:
  ThermalKValueReservoirWellTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:
  void SetUp() override
  {
    writeKValueTables();
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > >( "reservoirSystem" );

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

  GeosxState state;
  CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > * solver;
};

real64 constexpr ThermalKValueReservoirWellTest::time;
real64 constexpr ThermalKValueReservoirWellTest::dt;

TEST_F( ThermalKValueReservoirWellTest, assembleCoupledSystem )
{
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  CRSMatrix< real64, globalIndex > const & jacobian = solver->getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );
  residual.zero();
  jacobian.zero();
  solver->assembleSystem( time, dt, domain, solver->getDofManager(), jacobian.toViewConstSizes(), residual.toView() );
  EXPECT_GT( jacobian.numRows(), 0 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
