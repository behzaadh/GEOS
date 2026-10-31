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

/**
 * @file WaterSaturationPressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_WATERSATURATIONPRESSURE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_WATERSATURATIONPRESSURE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief Water saturation pressure from the IAPWS-IF97 Region 4 equation (273.15 K to critical point).
 * Outside that range the temperature is clamped and the derivative is zero.
 */
struct WaterSaturationPressure
{
  /**
   * @brief Calculate the saturation pressure of pure water
   * @param[in] temperature temperature [K]
   * @param[out] dPsat_dT derivative of the saturation pressure with respect to temperature [Pa/K]
   * @return saturation pressure [Pa]
   */
  GEOS_HOST_DEVICE
  static real64 compute( real64 const temperature, real64 & dPsat_dT )
  {
    constexpr real64 n1 = 0.11670521452767e4;
    constexpr real64 n2 = -0.72421316703206e6;
    constexpr real64 n3 = -0.17073846940092e2;
    constexpr real64 n4 = 0.12020824702470e5;
    constexpr real64 n5 = -0.32325550322333e7;
    constexpr real64 n6 = 0.14915108613530e2;
    constexpr real64 n7 = -0.48232657361591e4;
    constexpr real64 n8 = 0.40511340542057e6;
    constexpr real64 n9 = -0.23855557567849;
    constexpr real64 n10 = 0.65017534844798e3;
    constexpr real64 minTemperature = 273.15;
    constexpr real64 maxTemperature = 647.096;
    constexpr real64 MPaToPa = 1.0e6;

    dPsat_dT = 0.0;
    bool const inRange = (minTemperature < temperature) && (temperature < maxTemperature);
    real64 const T = inRange ? temperature : (temperature <= minTemperature ? minTemperature : maxTemperature);

    real64 const dT = T - n10;
    real64 const theta = T + n9 / dT;
    real64 const dTheta = 1.0 - n9 / (dT * dT);

    real64 const A = theta * theta + n1 * theta + n2;
    real64 const B = n3 * theta * theta + n4 * theta + n5;
    real64 const C = n6 * theta * theta + n7 * theta + n8;
    real64 const dA = 2.0 * theta + n1;
    real64 const dB = 2.0 * n3 * theta + n4;
    real64 const dC = 2.0 * n6 * theta + n7;

    real64 const D = sqrt( B * B - 4.0 * A * C );
    real64 const dD = (B * dB - 2.0 * (dA * C + A * dC)) / D;

    real64 const denominator = D - B;
    real64 const f = 2.0 * C / denominator;
    real64 const f2 = f * f;

    if( inRange )
    {
      real64 const df = (2.0 * dC * denominator - 2.0 * C * (dD - dB)) / (denominator * denominator);
      dPsat_dT = MPaToPa * 4.0 * f2 * f * df * dTheta;
    }
    return MPaToPa * f2 * f2;
  }
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_WATERSATURATIONPRESSURE_HPP_
