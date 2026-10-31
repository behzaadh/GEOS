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

/**
 * @file testMultiFluidKValueThermal.cpp
 */

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/HeatCapacityCoefficients.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/KValueFlashParameters.hpp"
#include "common/DataLayouts.hpp"
#include "common/initializeEnvironment.hpp"
#include "common/Units.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NUM_COMP >
static void populateArray( arraySlice1d< real64 > array, std::array< real64 const, NUM_COMP > const data )
{
  for( integer i = 0; i < NUM_COMP; ++i )
  {
    array[i] = data[i];
  }
}

static void fillHeatCapacity( dataRepository::Group & fluid, integer const numComps )
{
  using Keys = HeatCapacityCoefficients::viewKeyStruct;
  array2d< real64 > & coefficients = fluid.getReference< array2d< real64 > >( Keys::componentHeatCapacityCoefficientsString() );
  coefficients.resize( numComps, 5 );
  coefficients.zero();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    coefficients( ic, 0 ) = 35.0;
  }
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
class KValueThermalFluidTest : public FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >
{
public:
  using Base = FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >;
  using FluidModel = FLUID_TYPE;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-3;

  void createKValueTable( string const & tableName, real64 const k0 )
  {
    array1d< real64_array > coordinates( 2 );
    coordinates[0].emplace_back( 1.0e5 );
    coordinates[0].emplace_back( 6.0e7 );
    coordinates[1].emplace_back( 280.0 );
    coordinates[1].emplace_back( 420.0 );

    real64_array values( 4 );
    values[0] = k0;
    values[1] = k0 * 1.05;
    values[2] = k0 * 0.95;
    values[3] = k0 * 1.02;
    Base::createTable( tableName, coordinates, values );
  }

  void fillComponentProperties( FluidModel & fluid, std::unique_ptr< TestFluid< NUM_COMP > > const & testFluid )
  {
    using Keys = typename FluidModel::viewKeyStruct;
    string_array & componentNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    componentNames = testFluid->componentNames;
    array1d< real64 > & molarWeight = fluid.template getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    TestFluid< NUM_COMP >::createArray( molarWeight, testFluid->molecularWeight );
    array1d< real64 > & criticalPressure = fluid.template getReference< array1d< real64 > >( Keys::componentCriticalPressureString() );
    TestFluid< NUM_COMP >::createArray( criticalPressure, testFluid->criticalPressure );
    array1d< real64 > & criticalTemperature = fluid.template getReference< array1d< real64 > >( Keys::componentCriticalTemperatureString() );
    TestFluid< NUM_COMP >::createArray( criticalTemperature, testFluid->criticalTemperature );
    array1d< real64 > & acentricFactor = fluid.template getReference< array1d< real64 > >( Keys::componentAcentricFactorString() );
    TestFluid< NUM_COMP >::createArray( acentricFactor, testFluid->acentricFactor );
  }
};

class TwoPhaseKValueThermalTest : public KValueThermalFluidTest< CompositionalTwoPhaseKValueThermalFluidLohrenzBrayClark, 3, 2 >
{
public:
  TwoPhaseKValueThermalTest()
  {
    createKValueTable( "KV_CO2", 2.5 );
    createKValueTable( "KV_CH4", 8.0 );
    createKValueTable( "KV_C10", 0.05 );

    createFluid( "twoPhaseKValueThermal", [&]( FluidModel & fluid )
    {
      string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
      phaseNames = {"oil", "gas"};

      string_array & equationOfState = fluid.getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
      equationOfState = {"PengRobinson", "PengRobinson"};

      fillComponentProperties( fluid, TestFluid< 3 >::create( {Fluid::CO2, Fluid::CH4, Fluid::C10H22} ) );
      fillHeatCapacity( fluid, 3 );

      string_array & kValueTables = fluid.getReference< string_array >( KValueFlashParameters< 2 >::viewKeyStruct::kValueTablesString() );
      kValueTables = {"KV_CO2", "KV_CH4", "KV_C10"};
    } );
  }
};

class ThreePhaseKValueThermalTest : public KValueThermalFluidTest< CompositionalThreePhaseKValueThermalFluidLohrenzBrayClark, 3, 3 >
{
public:
  ThreePhaseKValueThermalTest()
  {
    createKValueTable( "KV_H2O", 0.2 );
    createKValueTable( "KV_CO2_V", 3.0 );
    createKValueTable( "KV_C10_V", 0.04 );
    createKValueTable( "KW_H2O", 60.0 );
    createKValueTable( "KW_CO2", 0.8 );
    createKValueTable( "KW_C10", 1.0e-3 );

    createFluid( "threePhaseKValueThermal", [&]( FluidModel & fluid )
    {
      string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
      phaseNames = {"oil", "gas", "water"};

      string_array & equationOfState = fluid.getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
      equationOfState = {"PengRobinson", "PengRobinson", "PengRobinson"};

      fillComponentProperties( fluid, TestFluid< 3 >::create( {Fluid::H2O, Fluid::CO2, Fluid::C10H22} ) );
      fillHeatCapacity( fluid, 3 );

      string_array & kValueTables = fluid.getReference< string_array >( KValueFlashParameters< 3 >::viewKeyStruct::kValueTablesString() );
      kValueTables = {"KV_H2O", "KV_CO2_V", "KV_C10_V", "KW_H2O", "KW_CO2", "KW_C10"};
    } );
  }
};

TEST_F( TwoPhaseKValueThermalTest, isThermal )
{
  EXPECT_TRUE( getFluid( "twoPhaseKValueThermal" )->isThermal() );
}

TEST_F( TwoPhaseKValueThermalTest, numericalDerivatives )
{
  FluidModel * fluid = getFluid( "twoPhaseKValueThermal" );
  fluid->setMassFlag( false );
  TestPoint const data( 5.0e6, units::convertCToK( 40.0 ), Feed< 3 >{0.3, 0.3, 0.4} );
  testNumericalDerivatives( fluid, data, 1.0e-7, relTol, absTol );
}

TEST_F( ThreePhaseKValueThermalTest, isThermal )
{
  EXPECT_TRUE( getFluid( "threePhaseKValueThermal" )->isThermal() );
}

TEST_F( ThreePhaseKValueThermalTest, numericalDerivatives )
{
  FluidModel * fluid = getFluid( "threePhaseKValueThermal" );
  fluid->setMassFlag( false );
  TestPoint const data( 5.0e6, units::convertCToK( 50.0 ), Feed< 3 >{0.4, 0.2, 0.4} );
  testNumericalDerivatives( fluid, data, 1.0e-7, relTol, absTol );
}

TEST_F( ThreePhaseKValueThermalTest, co2DissolvesInWater )
{
  FluidModel * fluid = getFluid( "threePhaseKValueThermal" );
  integer constexpr size = 1;
  m_parent.resize( size );
  fluid->allocateConstitutiveData( m_parent, 1 );

  array2d< real64, compflow::LAYOUT_PHASE > compositionArray( size, 3 );
  populateArray< 3 >( compositionArray[0], {0.4, 0.2, 0.4} );

  FluidWrapper fluidWrapper = fluid->createKernelWrapper();
  fluidWrapper.update( 0, 0, 5.0e6, units::convertCToK( 50.0 ), compositionArray[0] );

  arrayView3d< real64 const, multifluid::USD_PHASE > const phaseFraction = fluid->phaseFraction();
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseCompFraction = fluid->phaseCompFraction();
  integer constexpr aqueousIndex = 2;
  integer constexpr co2Index = 1;
  EXPECT_GT( phaseFraction[0][0][aqueousIndex], 1.0e-4 );
  EXPECT_GT( phaseCompFraction[0][0][aqueousIndex][co2Index], 1.0e-6 );
}

} // namespace testing
} // namespace geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::setupEnvironment( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::cleanupEnvironment();
  return result;
}
