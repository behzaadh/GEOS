/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "common/DataTypes.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/WaterSaturationPressure.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

// Reference values from IAPWS-IF97 (Table 35 of the release): T = 300, 500, 600 K
TEST( WaterSaturationPressureTest, referenceValues )
{
  real64 dPsat = 0.0;
  EXPECT_NEAR( WaterSaturationPressure::compute( 300.0, dPsat ), 0.353658941e4, 1.0e-6 * 0.353658941e4 );
  EXPECT_NEAR( WaterSaturationPressure::compute( 500.0, dPsat ), 0.263889776e7, 1.0e-6 * 0.263889776e7 );
  EXPECT_NEAR( WaterSaturationPressure::compute( 600.0, dPsat ), 0.123443146e8, 1.0e-6 * 0.123443146e8 );
}

TEST( WaterSaturationPressureTest, derivative )
{
  for( real64 const T : {280.0, 300.0, 373.15, 450.0, 550.0, 620.0} )
  {
    real64 dPsat = 0.0;
    real64 dummy = 0.0;
    WaterSaturationPressure::compute( T, dPsat );
    real64 const h = 1.0e-3;
    real64 const fd = (WaterSaturationPressure::compute( T + h, dummy ) - WaterSaturationPressure::compute( T - h, dummy )) / (2.0 * h);
    EXPECT_NEAR( dPsat, fd, 1.0e-5 * std::abs( fd ) );
  }
}

TEST( WaterSaturationPressureTest, outOfRangeIsClamped )
{
  real64 dPsat = 1.0;
  real64 dummy = 0.0;
  EXPECT_DOUBLE_EQ( WaterSaturationPressure::compute( 250.0, dPsat ), WaterSaturationPressure::compute( 273.15, dummy ) );
  EXPECT_DOUBLE_EQ( dPsat, 0.0 );
}

} // testing

} // geos
