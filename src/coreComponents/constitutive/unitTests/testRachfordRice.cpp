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


// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

static constexpr real64 relTol = 1e-4;

TEST( RachfordRiceTest, testRachfordRiceTwoComponents )
{
  constexpr integer numComps = 2;

  array1d< real64 > kValues( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  ////////////////////////////////////////

  kValues[0] = 1.12223;
  kValues[1] = 1.12223;

  feed[0] = 0.5;
  feed[1] = 0.5;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction1 = 1;
  real64 const vaporFraction1 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction1, expectedVaporFraction1, relTol );

  ////////////////////////////////////////

  kValues[0] = 56.1091;
  kValues[1] = 56.1091;

  feed[0] = 0.5;
  feed[1] = 0.5;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction2 = 1;
  real64 const vaporFraction2 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction2, expectedVaporFraction2, relTol );

  ////////////////////////////////////////

  kValues[0] = 54.866;
  kValues[1] = 54.866;

  feed[0] = 0.1;
  feed[1] = 0.9;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction3 = 1;
  real64 const vaporFraction3 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction3, expectedVaporFraction3, relTol );

  ////////////////////////////////////////

  kValues[0] = 1.09733;
  kValues[1] = 1.09733;

  feed[0] = 0.1;
  feed[1] = 0.9;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction4 = 1;
  real64 const vaporFraction4 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction4, expectedVaporFraction4, relTol );

  ////////////////////////////////////////

  kValues[0] = 0.9;
  kValues[1] = 1.09733;

  feed[0] = 1.0e-10;
  feed[1] = 1.0 - feed[0];

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction5 = 1;
  real64 const vaporFraction5 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction5, expectedVaporFraction5, relTol );

  ////////////////////////////////////////

  kValues[0] = 1.09733;
  kValues[1] = 0.9;

  feed[0] = 1.0e-10;
  feed[1] = 1.0 - feed[0];

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction6 = 0.0;
  real64 const vaporFraction6 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction6, expectedVaporFraction6, relTol );
}

TEST( RachfordRiceTest, testRachfordRiceFourComponents )
{
  constexpr integer numComps = 4;

  array1d< real64 > kValues( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  ////////////////////////////////////////

  kValues[0] = 537.526;
  kValues[1] = 0.00297518;
  kValues[2] = 6.19123e-08;
  kValues[3] = 1.44696;

  feed[0] = 0.099;
  feed[1] = 0.3;
  feed[2] = 0.6;
  feed[3] = 0.001;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction1 = 0.0975568;
  real64 const vaporFraction1 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction1, expectedVaporFraction1, relTol );

  ////////////////////////////////////////

  kValues[0] = 17.4329;
  kValues[1] = 0.000770753;
  kValues[2] = 1.18694e-06;
  kValues[3] = 0.00968376;

  feed[0] = 0.1;
  feed[1] = 0.1;
  feed[2] = 0.1;
  feed[3] = 0.7;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction2 = 0.045999;
  real64 const vaporFraction2 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction2, expectedVaporFraction2, relTol );

  ////////////////////////////////////////

  kValues[0] = 11.5506;
  kValues[1] = 0.000210682;
  kValues[2] = 1.57686e-08;
  kValues[3] = 0.0442361;

  feed[0] = 0.0984186;
  feed[1] = 0.297297;
  feed[2] = 0.593142;
  feed[3] = 0.0111427;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction3 = 0.0130259;
  real64 const vaporFraction3 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction3, expectedVaporFraction3, relTol );

}

TEST( RachfordRiceTest, testRachfordRiceThreePhase )
{
  constexpr integer numComps = 3;
  array1d< real64 > kVapour( numComps );
  array1d< real64 > kAqueous( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  // Constructed three-phase point: L=0.3, V=0.3, W=0.4
  kVapour[0] = 2.5;
  kVapour[1] = 0.4 / 0.6;
  kVapour[2] = 0.5;
  kAqueous[0] = 0.25;
  kAqueous[1] = 0.05 / 0.6;
  kAqueous[2] = 4.5;
  feed[0] = 0.23;
  feed[1] = 0.32;
  feed[2] = 0.45;
  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;

  real64 vapourFraction = 0.0;
  real64 aqueousFraction = 0.0;
  RachfordRice::solveThreePhase( kVapour.toSliceConst(),
                                 kAqueous.toSliceConst(),
                                 feed.toSliceConst(),
                                 presentComponentIds.toSliceConst(),
                                 vapourFraction,
                                 aqueousFraction );

  checkRelativeError( vapourFraction, 0.3, relTol );
  checkRelativeError( aqueousFraction, 0.4, relTol );

  real64 fV = 0.0;
  real64 fW = 0.0;
  real64 j11 = 0.0;
  real64 j12 = 0.0;
  real64 j21 = 0.0;
  real64 j22 = 0.0;
  RachfordRice::evaluateThreePhaseResidualAndJacobian( kVapour.toSliceConst(),
                                                       kAqueous.toSliceConst(),
                                                       feed.toSliceConst(),
                                                       presentComponentIds.toSliceConst(),
                                                       vapourFraction,
                                                       aqueousFraction,
                                                       fV, fW, j11, j12, j21, j22 );
  EXPECT_NEAR( fV, 0.0, 1.0e-8 );
  EXPECT_NEAR( fW, 0.0, 1.0e-8 );

  real64 const liquidFraction = 1.0 - vapourFraction - aqueousFraction;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const D = 1.0 + vapourFraction * ( kVapour[ic] - 1.0 ) + aqueousFraction * ( kAqueous[ic] - 1.0 );
    real64 const xi = feed[ic] / D;
    real64 const yi = kVapour[ic] * xi;
    real64 const wi = kAqueous[ic] * xi;
    checkRelativeError( liquidFraction * xi + vapourFraction * yi + aqueousFraction * wi, feed[ic], relTol );
  }
}

TEST( RachfordRiceTest, testRachfordRiceThreePhaseOilGasFallback )
{
  constexpr integer numComps = 2;
  array1d< real64 > kVapour( numComps );
  array1d< real64 > kAqueous( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  kVapour[0] = 5.0;
  kVapour[1] = 0.2;
  kAqueous[0] = 1.0e-3;
  kAqueous[1] = 1.0e-3;
  feed[0] = 0.4;
  feed[1] = 0.6;
  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 vapourFraction = 0.0;
  real64 aqueousFraction = 0.0;
  RachfordRice::solveThreePhase( kVapour.toSliceConst(),
                                 kAqueous.toSliceConst(),
                                 feed.toSliceConst(),
                                 presentComponentIds.toSliceConst(),
                                 vapourFraction,
                                 aqueousFraction );

  EXPECT_NEAR( aqueousFraction, 0.0, 1.0e-8 );
  real64 const twoPhaseVapour =
    RachfordRice::solve( kVapour.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );
  checkRelativeError( vapourFraction, twoPhaseVapour, relTol );
}

// Water-rich CO2/H2O mixture: the liquid (oil) phase must vanish continuously as the feed moves
// across the gas-water edge of the phase simplex. The returned phase fractions must match the
// two-phase gas-water Rachford-Rice solution once the oil phase disappears.
TEST( RachfordRiceTest, testRachfordRiceThreePhaseGasWaterFallback )
{
  constexpr integer numComps = 2;
  array1d< real64 > kVapour( numComps );
  array1d< real64 > kAqueous( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  // CO2: volatile, slightly soluble in water. H2O: low vapour K, strongly aqueous.
  kVapour[0] = 2.5;
  kVapour[1] = 8.0e-3;
  kAqueous[0] = 0.03;
  kAqueous[1] = 80.0;
  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  array1d< real64 > kGasWater( numComps );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    kGasWater[ic] = kVapour[ic] / kAqueous[ic];
  }

  real64 previousVapour = -1.0;
  real64 previousAqueous = -1.0;
  for( integer step = 1; step < 50; ++step )
  {
    feed[0] = 0.02 * step;
    feed[1] = 1.0 - feed[0];

    real64 vapourFraction = 0.0;
    real64 aqueousFraction = 0.0;
    RachfordRice::solveThreePhase( kVapour.toSliceConst(),
                                   kAqueous.toSliceConst(),
                                   feed.toSliceConst(),
                                   presentComponentIds.toSliceConst(),
                                   vapourFraction,
                                   aqueousFraction );

    // Fractions must stay inside the simplex
    EXPECT_GE( vapourFraction, 0.0 );
    EXPECT_GE( aqueousFraction, 0.0 );
    EXPECT_LE( vapourFraction + aqueousFraction, 1.0 + 1.0e-12 );

    // With no oil phase present the split must coincide with the gas-water Rachford-Rice solution
    real64 const liquidFraction = 1.0 - vapourFraction - aqueousFraction;
    if( liquidFraction < 1.0e-10 )
    {
      real64 const twoPhaseVapour =
        RachfordRice::solve( kGasWater.toSliceConst(),
                             feed.toSliceConst(),
                             presentComponentIds.toSliceConst() );
      real64 const expectedVapour = LvArray::math::min( LvArray::math::max( twoPhaseVapour, 0.0 ), 1.0 );
      EXPECT_NEAR( vapourFraction, expectedVapour, 1.0e-8 );
    }

    // The solution must be continuous in the feed (no jumps larger than the feed increment scale)
    if( previousVapour >= 0.0 )
    {
      EXPECT_LT( LvArray::math::abs( vapourFraction - previousVapour ), 0.1 );
      EXPECT_LT( LvArray::math::abs( aqueousFraction - previousAqueous ), 0.1 );
    }
    previousVapour = vapourFraction;
    previousAqueous = aqueousFraction;
  }
}
