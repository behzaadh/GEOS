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
 * @file RachfordRice.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

struct RachfordRice
{
public:
  /// Tolerance of the SSI loop
  static constexpr real64 SSITolerance = MultiFluidConstants::SSITolerance;
  /// Tolerance of the Newton loop
  static constexpr real64 newtonTolerance = MultiFluidConstants::newtonTolerance;
  /// Max number of SSI iterations
  static constexpr integer maxSSIIterations = MultiFluidConstants::maxSSIIterations;
  /// Max number of Newton iterations
  static constexpr integer maxNewtonIterations = MultiFluidConstants::maxNewtonIterations;
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
  /// Max number of step halvings in the three-phase backtracking
  static constexpr integer maxLineSearchCuts = 60;
  /// Bounding box for the three-phase negative-flash iterate
  static constexpr real64 maxNegativeFlashFraction = 2.0;
  /// Smallest phase fraction for a phase to be considered present in the three-phase solve
  static constexpr real64 minPhaseFraction = 1.0e-12;

  /**
   * @brief Function solving the Rachford-Rice equation
   * @input[in] kValues the array of K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @return the gas mole fraction
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  solve( arraySlice1d< real64 const, USD2 > const & kValues,
         arraySlice1d< real64 const, USD1 > const & feed,
         arraySlice1d< integer const > const & presentComponentIds )
  {
    real64 gasPhaseMoleFraction = 0;

    // min and max Kvalues for non-zero composition
    real64 minK, maxK;
    findKValueRange( kValues, presentComponentIds, minK, maxK );

    // check for trivial solutions.
    // this corresponds to bad Kvalues
    if( maxK < 1.0 )
    {
      return 0.0;
    }
    if( minK > 1.0 )
    {
      return 1.0;
    }

    // start the solver loop

    // step 1: find solution window
    real64 xMin = 1.0 / ( 1 - maxK );
    real64 xMax = 1.0 / ( 1 - minK );
    real64 const sqrtEpsilon = sqrt( epsilon );
    xMin += sqrtEpsilon * ( LvArray::math::abs( xMin ) + sqrtEpsilon );
    xMax -= sqrtEpsilon * ( LvArray::math::abs( xMax ) + sqrtEpsilon );

    real64 currentError = 1 / epsilon;

    // step 2: start the SSI loop
    // Evaluate at the bounds
    real64 funcXMin = evaluate( kValues, feed, presentComponentIds, xMin );
    real64 funcXMax = evaluate( kValues, feed, presentComponentIds, xMax );
    real64 funcXMid = 0.0;

    // If the bound values are the same sign then we have a trivial solution
    if( 0.0 < funcXMin * funcXMax )
    {
      gasPhaseMoleFraction = (0.0 < funcXMin) ? 1.0 : 0.0;
      return gasPhaseMoleFraction;
    }

    integer SSIIteration = 0;

    while( ( currentError > SSITolerance ) && ( SSIIteration < maxSSIIterations ) )
    {
      real64 const xMid = 0.5 * ( xMin + xMax );
      funcXMid = evaluate( kValues, feed, presentComponentIds, xMid );

      if( 0.0 < funcXMax * funcXMid )
      {
        xMax = xMid;
        funcXMax = funcXMid;
      }
      else if( 0.0 < funcXMin * funcXMid )
      {
        xMin = xMid;
        funcXMin = funcXMid;
      }

      currentError = LvArray::math::min( LvArray::math::abs( funcXMid ),
                                         LvArray::math::abs( xMax - xMin ) );
      SSIIteration++;

      // TODO: add warning if max number of SSI iterations is reached
    }

    gasPhaseMoleFraction = 0.5 * ( xMax + xMin );

    // step 3: start the Newton loop
    integer newtonIteration = 0;
    real64 newtonValue = gasPhaseMoleFraction;
    real64 funcNewton = evaluate( kValues, feed, presentComponentIds, newtonValue );

    while( ( currentError > newtonTolerance ) && ( newtonIteration < maxNewtonIterations ) )
    {
      real64 deltaNewton = -funcNewton / evaluateDerivative( kValues, feed, presentComponentIds, newtonValue );

      // test if we are stepping out of the [xMin;xMax] interval
      if( newtonValue + deltaNewton < xMin )
      {
        deltaNewton = 0.5 * ( xMin - newtonValue );
      }
      else if( newtonValue + deltaNewton > xMax )
      {
        deltaNewton = 0.5 * ( xMax - newtonValue );
      }

      newtonValue = newtonValue + deltaNewton;

      funcNewton = evaluate( kValues, feed, presentComponentIds, newtonValue );

      currentError = LvArray::math::min( LvArray::math::abs( funcNewton ),
                                         LvArray::math::abs( deltaNewton ) );
      newtonIteration++;

      // TODO: add warning if max number of Newton iterations is reached
    }
    return gasPhaseMoleFraction = newtonValue;
  }

  /**
   * @brief Solve the three-phase Rachford-Rice equations for vapour and aqueous mole fractions.
   * @details K-values are relative to the liquid (oil) phase: K^V_i = y_i / x_i, K^W_i = w_i / x_i.
   *          Falls back to a two-phase Rachford-Rice solve when a phase is absent.
   * @param[in] kVapour K-values of the vapour phase relative to liquid
   * @param[in] kAqueous K-values of the aqueous phase relative to liquid
   * @param[in] feed overall (feed) component mole fractions
   * @param[in] presentComponentIds indices of components with a non-zero feed
   * @param[out] vapourFraction vapour mole fraction V
   * @param[out] aqueousFraction aqueous mole fraction W
   */
  template< integer USD1, integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static void
  solveThreePhase( arraySlice1d< real64 const, USD2 > const & kVapour,
                   arraySlice1d< real64 const, USD3 > const & kAqueous,
                   arraySlice1d< real64 const, USD1 > const & feed,
                   arraySlice1d< integer const > const & presentComponentIds,
                   real64 & vapourFraction,
                   real64 & aqueousFraction )
  {
    // The three-phase Rachford-Rice equations f_V = f_W = 0 are the stationarity conditions of
    //   F(V,W) = -sum_i z_i ln( D_i ),   D_i = 1 + V (K_i^V - 1) + W (K_i^W - 1),
    // which is convex on the polytope { D_i > 0 } (negative-flash domain). The phase split is
    // therefore the minimiser of F over the closed simplex { V >= 0, W >= 0, V + W <= 1 }:
    //  1. Newton on F from the simplex centroid (always feasible) with feasibility backtracking.
    //     If the unconstrained minimiser is strictly inside the simplex, three phases coexist.
    //  2. Otherwise the constrained minimiser lies on the simplex boundary. Each edge is a
    //     two-phase Rachford-Rice problem, solved with the robust 1D solver and clamped to the
    //     edge; the edge with the lowest F is selected.
    // Selecting the constrained minimiser of a convex function makes the phase fractions
    // continuous across phase appearance/disappearance, which is required for the Newton
    // solver of the flow problem to behave.

    real64 V = 1.0 / 3.0;
    real64 W = 1.0 / 3.0;
    bool converged = false;

    real64 fV = 0.0;
    real64 fW = 0.0;
    real64 j11 = 0.0;
    real64 j12 = 0.0;
    real64 j21 = 0.0;
    real64 j22 = 0.0;

    real64 objective = evaluateThreePhaseObjective( kVapour, kAqueous, feed, presentComponentIds, V, W );

    for( integer newtonIteration = 0; newtonIteration < maxNewtonIterations; ++newtonIteration )
    {
      evaluateThreePhaseResidualAndJacobian( kVapour, kAqueous, feed, presentComponentIds,
                                             V, W, fV, fW, j11, j12, j21, j22 );
      real64 const currentError = LvArray::math::max( LvArray::math::abs( fV ), LvArray::math::abs( fW ) );
      if( currentError <= newtonTolerance )
      {
        converged = true;
        break;
      }

      real64 const det = j11 * j22 - j12 * j21;
      if( LvArray::math::abs( det ) < epsilon )
      {
        break;
      }
      real64 const invDet = 1.0 / det;
      real64 const dV = ( -j22 * fV + j12 * fW ) * invDet;
      real64 const dW = ( j21 * fV - j11 * fW ) * invDet;

      // Backtrack to stay strictly inside the negative-flash domain, inside a bounding box,
      // and to guarantee a decrease of the convex objective.
      real64 alpha = 1.0;
      bool accepted = false;
      real64 trialV = V;
      real64 trialW = W;
      real64 trialObjective = objective;
      for( integer cut = 0; cut < maxLineSearchCuts; ++cut )
      {
        trialV = V + alpha * dV;
        trialW = W + alpha * dW;
        if( isThreePhaseFeasible( kVapour, kAqueous, presentComponentIds, trialV, trialW ) )
        {
          trialObjective = evaluateThreePhaseObjective( kVapour, kAqueous, feed, presentComponentIds, trialV, trialW );
          if( trialObjective <= objective + newtonTolerance )
          {
            accepted = true;
            break;
          }
        }
        alpha *= 0.5;
      }
      if( !accepted )
      {
        break;
      }

      real64 const stepSize = LvArray::math::max( LvArray::math::abs( trialV - V ), LvArray::math::abs( trialW - W ) );
      V = trialV;
      W = trialW;
      objective = trialObjective;
      if( stepSize <= newtonTolerance )
      {
        // Require a true RR root before accepting an interior three-phase state.
        // SSITolerance (1e-3) is far too loose for the subsequent IFT derivatives in
        // KValueFlashModel, which assume f_V = f_W = 0 and will otherwise poison the
        // flow Newton Jacobian (period-2 stalls in Rvol/Renergy).
        evaluateThreePhaseResidualAndJacobian( kVapour, kAqueous, feed, presentComponentIds,
                                               V, W, fV, fW, j11, j12, j21, j22 );
        converged = LvArray::math::max( LvArray::math::abs( fV ), LvArray::math::abs( fW ) ) <= newtonTolerance;
        break;
      }
    }

    real64 const L = 1.0 - V - W;
    if( converged && V > minPhaseFraction && W > minPhaseFraction && L > minPhaseFraction )
    {
      vapourFraction = V;
      aqueousFraction = W;
      return;
    }

    // The constrained minimiser is on the boundary of the simplex: compare the three edges.

    // Edge W = 0: liquid-vapour split
    real64 const vOilGas = clampUnit( solve( kVapour, feed, presentComponentIds ) );
    real64 const fOilGas = evaluateThreePhaseObjective( kVapour, kAqueous, feed, presentComponentIds, vOilGas, 0.0 );

    // Edge V = 0: liquid-aqueous split
    real64 const wOilWater = clampUnit( solve( kAqueous, feed, presentComponentIds ) );
    real64 const fOilWater = evaluateThreePhaseObjective( kVapour, kAqueous, feed, presentComponentIds, 0.0, wOilWater );

    // Edge L = 0: vapour-aqueous split with K^{GW} = K^V / K^W
    integer const numComps = kVapour.size();
    stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > kGasWater( numComps );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      kGasWater[ic] = kVapour[ic] / LvArray::math::max( kAqueous[ic], epsilon );
    }
    real64 const vGasWater = clampUnit( solve( kGasWater.toSliceConst(), feed, presentComponentIds ) );
    real64 const fGasWater = evaluateThreePhaseObjective( kVapour, kAqueous, feed, presentComponentIds, vGasWater, 1.0 - vGasWater );

    if( fOilGas <= fOilWater && fOilGas <= fGasWater )
    {
      vapourFraction = vOilGas;
      aqueousFraction = 0.0;
    }
    else if( fOilWater <= fGasWater )
    {
      vapourFraction = 0.0;
      aqueousFraction = wOilWater;
    }
    else
    {
      vapourFraction = vGasWater;
      aqueousFraction = 1.0 - vGasWater;
    }
  }

  /**
   * @brief Evaluate the convex objective F(V,W) = -sum_i z_i ln(D_i) whose gradient is the
   *        (negated) three-phase Rachford-Rice residual.
   */
  template< integer USD1, integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static real64
  evaluateThreePhaseObjective( arraySlice1d< real64 const, USD2 > const & kVapour,
                               arraySlice1d< real64 const, USD3 > const & kAqueous,
                               arraySlice1d< real64 const, USD1 > const & feed,
                               arraySlice1d< integer const > const & presentComponentIds,
                               real64 const vapourFraction,
                               real64 const aqueousFraction )
  {
    real64 value = 0.0;
    integer const numPresent = presentComponentIds.size();
    for( integer i = 0; i < numPresent; ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const D = 1.0 + vapourFraction * ( kVapour[ic] - 1.0 ) + aqueousFraction * ( kAqueous[ic] - 1.0 );
      value -= feed[ic] * LvArray::math::log( LvArray::math::max( D, epsilon ) );
    }
    return value;
  }

  /**
   * @brief Evaluate the three-phase Rachford-Rice residuals and 2x2 Jacobian at (V, W).
   */
  template< integer USD1, integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static void
  evaluateThreePhaseResidualAndJacobian( arraySlice1d< real64 const, USD2 > const & kVapour,
                                         arraySlice1d< real64 const, USD3 > const & kAqueous,
                                         arraySlice1d< real64 const, USD1 > const & feed,
                                         arraySlice1d< integer const > const & presentComponentIds,
                                         real64 const vapourFraction,
                                         real64 const aqueousFraction,
                                         real64 & fV,
                                         real64 & fW,
                                         real64 & j11,
                                         real64 & j12,
                                         real64 & j21,
                                         real64 & j22 )
  {
    fV = 0.0;
    fW = 0.0;
    j11 = 0.0;
    j12 = 0.0;
    j21 = 0.0;
    j22 = 0.0;
    integer const numPresent = presentComponentIds.size();
    for( integer i = 0; i < numPresent; ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const dV = kVapour[ic] - 1.0;
      real64 const dW = kAqueous[ic] - 1.0;
      real64 const D = LvArray::math::max( 1.0 + vapourFraction * dV + aqueousFraction * dW, epsilon );
      real64 const invD = 1.0 / D;
      real64 const xTerm = feed[ic] * invD;
      fV += dV * xTerm;
      fW += dW * xTerm;
      real64 const invD2 = invD * invD;
      real64 const zInvD2 = feed[ic] * invD2;
      j11 -= zInvD2 * dV * dV;
      j12 -= zInvD2 * dV * dW;
      j21 -= zInvD2 * dW * dV;
      j22 -= zInvD2 * dW * dW;
    }
  }

private:
  /**
   * @brief Clamp a phase fraction to [0, 1].
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 clampUnit( real64 const value )
  {
    return LvArray::math::min( LvArray::math::max( value, 0.0 ), 1.0 );
  }

  /**
   * @brief Check that (V, W) is strictly inside the negative-flash domain { D_i > 0 } and
   *        inside the bounding box used to keep the Newton iterate finite.
   */
  template< integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static bool isThreePhaseFeasible( arraySlice1d< real64 const, USD2 > const & kVapour,
                                    arraySlice1d< real64 const, USD3 > const & kAqueous,
                                    arraySlice1d< integer const > const & presentComponentIds,
                                    real64 const vapourFraction,
                                    real64 const aqueousFraction )
  {
    if( LvArray::math::abs( vapourFraction ) > maxNegativeFlashFraction ||
        LvArray::math::abs( aqueousFraction ) > maxNegativeFlashFraction )
    {
      return false;
    }
    integer const numPresent = presentComponentIds.size();
    for( integer i = 0; i < numPresent; ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const D = 1.0 + vapourFraction * ( kVapour[ic] - 1.0 ) + aqueousFraction * ( kAqueous[ic] - 1.0 );
      if( D <= epsilon )
      {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Function evaluating the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the Rachford-Rice function is evaluated
   * @return the value of the Rachford-Rice function at x
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  evaluate( arraySlice1d< real64 const, USD2 > const & kValues,
            arraySlice1d< real64 const, USD1 > const & feed,
            arraySlice1d< integer const > const & presentComponentIds,
            real64 const & x )
  {
    real64 value = 0.0;
    for( integer i = 0; i < presentComponentIds.size(); ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const k = ( kValues[ic] - 1.0 );
      value += feed[ic] * k / ( 1.0 + x * k );
    }
    return value;
  }

  /**
   * @brief Function evaluating the derivative of the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the derivative of the Rachford-Rice function is evaluated
   * @return the value of the derivative of the Rachford-Rice function at x
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  evaluateDerivative( arraySlice1d< real64 const, USD2 > const & kValues,
                      arraySlice1d< real64 const, USD1 > const & feed,
                      arraySlice1d< integer const > const & presentComponentIds,
                      real64 const & x )
  {
    real64 value = 0.0;
    for( integer i = 0; i < presentComponentIds.size(); ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const k = ( kValues[ic] - 1.0 );
      real64 const r = k / ( 1.0 + x * k );
      value -= feed[ic] * r * r;
    }
    return value;
  }

  /**
   * @brief Calculate the minimum and maximum k-value
   * @input[in] kValues the array fo K-values
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[out] minK the minimum k-value for non-zero components
   * @input[out] maxK the maximum k-value for non-zero components
   **/
  template< integer USD >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void
  static
  findKValueRange( arraySlice1d< real64 const, USD > const & kValues,
                   arraySlice1d< integer const > const & presentComponentIds,
                   real64 & minK,
                   real64 & maxK )
  {
    minK = 1.0 / epsilon;
    maxK = 0.0;
    for( integer const ic : presentComponentIds )
    {
      if( kValues[ic] > maxK )
      {
        maxK = kValues[ic];
      }
      if( kValues[ic] < minK )
      {
        minK = kValues[ic];
      }
    }
  }

};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
