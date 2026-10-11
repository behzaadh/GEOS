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
 * @file KValueFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/Utilities.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
class KValueFlashParameters;

template< integer NUM_PHASE >
class KValueFlashModelUpdate final : public FunctionBaseUpdate
{
  static constexpr integer numPhases = NUM_PHASE;
  // Must be 2-phase or 3-phase
  static_assert( NUM_PHASE == 2 || NUM_PHASE == 3, "KValue flash must be 2-phase or 3-phase" );

public:
  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using Deriv = constitutive::multifluid::DerivativeOffset;

  KValueFlashModelUpdate( integer const numComponents,
                          integer const liquidIndex,
                          integer const vapourIndex,
                          integer const aqueousIndex,
                          TableFunction const & pressureTable,
                          TableFunction const & temperatureTable,
                          arrayView1d< integer const > const & presentComponents,
                          arrayView4d< real64 const > const & kValues );

  // Mark as a 2-phase or 3-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return numPhases; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

private:
  /**
   * @brief Fill phase fractions and compositions for a two-phase K-value split.
   * @param[in] primaryIndex index of the liquid-like (K-value denominator) phase
   * @param[in] secondaryIndex index of the vapour-like (K-value numerator) phase
   * @param[in] absentIndex index of an absent third phase, or -1
   * @param[in] kValues K-values of the secondary phase relative to the primary
   * @param[in] dKdP pressure derivatives of the K-values
   * @param[in] dKdT temperature derivatives of the K-values
   * @param[in] compFraction overall composition
   * @param[in] secondaryFraction mole fraction of the secondary phase
   * @param[in,out] phaseFraction phase mole fractions and derivatives
   * @param[in,out] phaseCompFraction phase compositions and derivatives
   */
  template< int USD1, int USD2, int USD3, int USD4 >
  GEOS_HOST_DEVICE
  void computeTwoPhaseFromKValues( integer const primaryIndex,
                                   integer const secondaryIndex,
                                   integer const absentIndex,
                                   arraySlice1d< real64 const, USD1 > const & kValues,
                                   arraySlice1d< real64 const, USD2 > const & dKdP,
                                   arraySlice1d< real64 const, USD3 > const & dKdT,
                                   arraySlice1d< real64 const, USD4 > const & compFraction,
                                   real64 secondaryFraction,
                                   PhaseProp::SliceType const phaseFraction,
                                   PhaseComp::SliceType const phaseCompFraction ) const;

  /**
   * @brief Assign a single-phase state (all phase compositions equal the feed).
   */
  template< int USD1 >
  GEOS_HOST_DEVICE
  void assignSinglePhase( integer const presentPhaseIndex,
                          arraySlice1d< real64 const, USD1 > const & compFraction,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseComp::SliceType const phaseCompFraction ) const;

  integer const m_liquidIndex;
  integer const m_vapourIndex;
  integer const m_aquoesIndex;

  integer const m_numComponents{0};
  integer const m_numPressurePoints{0};
  integer const m_numTemperaturePoints{0};

  /// Table with pressure indices
  TableFunction::KernelWrapper m_pressureTable;

  /// Table with temperature indices
  TableFunction::KernelWrapper m_temperatureTable;

  /// Index of available components
  arrayView1d< integer const > m_presentComponents;

  /// Hypercube of k-values
  arrayView4d< real64 const > m_kValues;
};

template< integer NUM_PHASE >
class KValueFlashModel : public FunctionBase
{
public:
  KValueFlashModel( string const & name,
                    ComponentProperties const & componentProperties,
                    ModelParameters const & modelParameters,
                    arrayView1d< integer const > const phaseTypes );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = KValueFlashModelUpdate< NUM_PHASE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

  // Determine phase ordering
  static void calculatePhaseOrdering( arrayView1d< integer const > const & phaseTypes,
                                      arrayView1d< integer > const & phaseOrder );
private:
  /// Index of present component
  array1d< integer > m_presentComponents;

  KValueFlashParameters< NUM_PHASE > const * m_parameters{};
  arrayView1d< integer const > const m_phaseTypes;
};

template< integer NUM_PHASE >
template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void
KValueFlashModelUpdate< NUM_PHASE >::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                              real64 const & pressure,
                                              real64 const & temperature,
                                              arraySlice1d< real64 const, USD1 > const & compFraction,
                                              arraySlice2d< real64, USD2 > const & kValues,
                                              PhaseProp::SliceType const phaseFraction,
                                              PhaseComp::SliceType const phaseCompFraction ) const
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numDofs = 2 + m_numComponents;

  // This code block is responsible for retrieving the K-equilibrium coefficients
  // ($K_i = y_i / x_i$) for all components in all secondary phases, and calculating
  // their partial derivatives with respect to pressure (P) and temperature (T).
  //
  // Since the K-values are stored in a 4D hypercube (`m_kValues`) defined on a
  // non-uniform grid of pressure and temperature points, bilinear interpolation is used.
  //
  // The overall process is:
  // 1. Use pre-built index tables (`m_pressureTable`, `m_temperatureTable`) to map the physical
  // `pressure` and `temperature` values to normalized, fractional coordinates ($p_{\alpha}$, $t_{\alpha}$). These
  // tables return:
  // - The integer part of the coordinate, which is the required lower **grid index** ($p_n, t_n$).
  // - The fractional part of the coordinate, which is the **interpolation weight** ($p_a, t_a$).
  // - The derivative of the normalized coordinate with respect to the physical
  // coordinate, required for the **chain rule** in derivative calculation ($dp_{\alpha}/dP, dt_{\alpha}/dT$).
  // 2. Loop through all secondary phases and components, read the four bounding K-values from the
  // hypercube, and apply the bilinear formula to calculate the interpolated $K_i$ and its
  // derivatives $dK_i/dP$ and $dK_i/dT$.

  real64 pa = 0.0;
  real64 pa_dp = 0.0;
  pa = m_pressureTable.compute( &pressure, &pa_dp );
  // The integer part of the normalized coordinate is the lower pressure index (pn) for interpolation.
  integer const pn = LvArray::math::min( static_cast< integer >(pa), m_numPressurePoints - 2 );
  // The fractional part is the interpolation weight (pa).
  pa -= pn;

  real64 ta = 0.0;
  real64 ta_dt = 0.0;
  ta = m_temperatureTable.compute( &temperature, &ta_dt );
  // The integer part of the normalized coordinate is the lower temperature index (tn) for interpolation.
  integer const tn = LvArray::math::min( static_cast< integer >(ta), m_numTemperaturePoints - 2 );
  // The fractional part is the interpolation weight (ta).
  ta -= tn;

  // Allocate temporary space for derivatives: 2 (for P and T) x (NUM_PHASE-1) x (m_numComponents)
  stackArray3d< real64, 2*(NUM_PHASE-1)*MultiFluidConstants::MAX_NUM_COMPONENTS > kValueDerivatives( 2, NUM_PHASE-1, m_numComponents );
  arraySlice3d< real64 > dK = kValueDerivatives.toSlice();
  for( integer ip = 0; ip < NUM_PHASE-1; ++ip )
  {
    for( integer ic = 0; ic < m_numComponents; ic++ )
    {
      // Retrieve the four bounding K-values for the 2x2 stencil at the grid indices (pn, tn)
      real64 const k00 = m_kValues( ip, ic, pn, tn );     // K at (P_n, T_n)
      real64 const k01 = m_kValues( ip, ic, pn, tn+1 );   // K at (P_n, T_{n+1})
      real64 const k10 = m_kValues( ip, ic, pn+1, tn );   // K at (P_{n+1}, T_n)
      real64 const k11 = m_kValues( ip, ic, pn+1, tn+1 ); // K at (P_{n+1}, T_{n+1})

      kValues( ip, ic ) = (1.0-pa)*(1.0-ta)*k00 + (1.0-pa)*ta*k01 + pa*(1.0-ta)*k10 + pa*ta*k11;

      // Derivative with respect to pressure (P): dK/dP = dK/d(pa) * d(pa)/dP
      dK( Deriv::dP, ip, ic ) = pa_dp*( -(1.0-ta)*k00 - ta*k01 + (1.0-ta)*k10 + ta*k11 );

      // Derivative with respect to temperature (T): dK/dT = dK/d(ta) * d(ta)/dT
      dK( Deriv::dT, ip, ic ) = ta_dt*( -(1.0-pa)*k00 + (1.0-pa)*k01 - pa*k10 + pa*k11 );
    }
  }

  // Find out which components are active
  stackArray1d< integer, MultiFluidConstants::MAX_NUM_COMPONENTS > componentIndices( m_numComponents );
  calculatePresentComponents( m_numComponents, compFraction, componentIndices );
  auto const presentComponents = componentIndices.toSliceConst();

  if constexpr (NUM_PHASE == 2)
  {
    // Check if the Rachford-Rice will actually solve
    // If all the k-values are on one side of unity then the Rachford-Rice will return a trivial solution
    integer constexpr LIQUID_PRESENT = 1;
    integer constexpr VAPOUR_PRESENT = 2;
    integer phasesPresent = 0;
    for( integer const ic : presentComponents )
    {
      real64 const k = kValues( 0, ic ) - 1.0;
      phasesPresent |= ((k < 0.0) ? LIQUID_PRESENT : VAPOUR_PRESENT);
    }

    real64 vapourFraction = 0.0;
    if( phasesPresent == (VAPOUR_PRESENT|LIQUID_PRESENT) )
    {
      vapourFraction = RachfordRice::solve( kValues[0].toSliceConst(), compFraction, presentComponents );

      // Calculate derivatives implicitly from the Rachford-Rice equation
      real64 denominator = 0.0;
      real64 pressureNumerator = 0.0;
      real64 temperatureNumerator = 0.0;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const k = kValues( 0, ic ) - 1.0;
        real64 const r = 1.0 / ( 1.0 + vapourFraction*k );
        pressureNumerator += compFraction[ic] * dK( Deriv::dP, 0, ic ) * r * r;
        temperatureNumerator += compFraction[ic] * dK( Deriv::dT, 0, ic ) * r * r;
        denominator += compFraction[ic] * k * k * r * r;
        phaseFraction.derivs( m_vapourIndex, Deriv::dC + ic ) = k * r;
      }

      GEOS_ERROR_IF( denominator < MultiFluidConstants::epsilon,
                     "Failed to calculate derivatives for the Rachford-Rice equation." );
      real64 const invDenominator = 1.0 / denominator;
      phaseFraction.derivs( m_vapourIndex, Deriv::dP ) = pressureNumerator * invDenominator;
      phaseFraction.derivs( m_vapourIndex, Deriv::dT ) = temperatureNumerator * invDenominator;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseFraction.derivs( m_vapourIndex, Deriv::dC + ic ) *= invDenominator;
      }

      // Calculate phase compositions
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const k = kValues( 0, ic ) - 1.0;
        real64 const r = 1.0 / ( 1.0 + vapourFraction*k );
        real64 const xi = compFraction[ic] * r;
        real64 const yi = xi * kValues( 0, ic );
        phaseCompFraction.value( m_liquidIndex, ic ) = xi;
        phaseCompFraction.value( m_vapourIndex, ic ) = yi;

        for( integer const idof : {Deriv::dP, Deriv::dT} )
        {
          real64 const dV = phaseFraction.derivs( m_vapourIndex, idof );
          real64 const dKi = dK( idof, 0, ic );
          real64 const dxi = r*(-vapourFraction*xi*dKi + (xi-yi)*dV);
          phaseCompFraction.derivs( m_liquidIndex, ic, idof ) = dxi;
          phaseCompFraction.derivs( m_vapourIndex, ic, idof ) = xi*dKi + dxi*kValues( 0, ic );
        }
        for( integer jc = 0; jc < m_numComponents; ++jc )
        {
          integer const idof = Deriv::dC + jc;

          real64 const dV = phaseFraction.derivs( m_vapourIndex, idof );
          real64 const dz = (jc == ic) ? 1.0 : 0.0;
          real64 const dxi = r*(dz + (xi-yi)*dV);
          phaseCompFraction.derivs( m_liquidIndex, ic, idof ) = dxi;
          phaseCompFraction.derivs( m_vapourIndex, ic, idof ) = dxi*kValues( 0, ic );
        }
      }

      if( vapourFraction < MultiFluidConstants::epsilon || 1.0-vapourFraction < MultiFluidConstants::epsilon )
      {
        vapourFraction = LvArray::math::min( LvArray::math::max( vapourFraction, 0.0 ), 1.0 );
        phaseFraction.value[m_vapourIndex] = vapourFraction;
        LvArray::forValuesInSlice( phaseFraction.derivs[m_vapourIndex], setZero );

        integer const singlePhaseIndex = vapourFraction < MultiFluidConstants::epsilon ? m_liquidIndex : m_vapourIndex;
        LvArray::forValuesInSlice( phaseCompFraction.derivs[singlePhaseIndex], setZero );
        for( integer ic = 0; ic < m_numComponents; ++ic )
        {
          phaseCompFraction.value( singlePhaseIndex, ic ) = compFraction[ic];
          phaseCompFraction.derivs( singlePhaseIndex, ic, Deriv::dC + ic ) = 1.0;
        }
      }
    }
    else
    {
      LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
      LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );
      vapourFraction = (phasesPresent == LIQUID_PRESENT) ? 0.0 : 1.0;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseCompFraction.value( m_liquidIndex, ic ) = compFraction[ic];
        phaseCompFraction.value( m_vapourIndex, ic ) = compFraction[ic];
        phaseCompFraction.derivs( m_liquidIndex, ic, Deriv::dC + ic ) = 1.0;
        phaseCompFraction.derivs( m_vapourIndex, ic, Deriv::dC + ic ) = 1.0;
      }
    }

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_vapourIndex] = vapourFraction;
    phaseFraction.value[m_liquidIndex] = 1.0 - vapourFraction;
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }
  else
  {
    integer constexpr vapourKIndex = 0;
    integer constexpr aqueousKIndex = 1;
    integer const liquidIndex = m_liquidIndex;
    integer const vapourIndex = m_vapourIndex;
    integer const aqueousIndex = m_aquoesIndex;

    LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
    LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );

    real64 vapourFraction = 0.0;
    real64 aqueousFraction = 0.0;
    RachfordRice::solveThreePhase( kValues[vapourKIndex].toSliceConst(),
                                   kValues[aqueousKIndex].toSliceConst(),
                                   compFraction,
                                   presentComponents,
                                   vapourFraction,
                                   aqueousFraction );

    // Align presence with RachfordRice::minPhaseFraction so the derivative branch
    // matches the flash classification (machine-epsilon chatter caused 2P/3P flips).
    real64 const phasePresenceTol = RachfordRice::minPhaseFraction;

    vapourFraction = LvArray::math::min( LvArray::math::max( vapourFraction, 0.0 ), 1.0 );
    aqueousFraction = LvArray::math::min( LvArray::math::max( aqueousFraction, 0.0 ), 1.0 );
    if( vapourFraction + aqueousFraction > 1.0 )
    {
      real64 const scale = 1.0 / ( vapourFraction + aqueousFraction );
      vapourFraction *= scale;
      aqueousFraction *= scale;
    }
    // Snap numerically absent phases to exact zero so edge derivatives stay consistent.
    if( vapourFraction < phasePresenceTol )
    {
      vapourFraction = 0.0;
    }
    if( aqueousFraction < phasePresenceTol )
    {
      aqueousFraction = 0.0;
    }
    real64 liquidFraction = 1.0 - vapourFraction - aqueousFraction;
    if( liquidFraction < phasePresenceTol )
    {
      liquidFraction = 0.0;
      if( vapourFraction + aqueousFraction > MultiFluidConstants::epsilon )
      {
        real64 const scale = 1.0 / ( vapourFraction + aqueousFraction );
        vapourFraction *= scale;
        aqueousFraction *= scale;
      }
    }

    bool const vapourPresent = vapourFraction > phasePresenceTol;
    bool const aqueousPresent = aqueousFraction > phasePresenceTol;
    bool const liquidPresent = liquidFraction > phasePresenceTol;
    integer const numPresentPhases = ( vapourPresent ? 1 : 0 ) + ( aqueousPresent ? 1 : 0 ) + ( liquidPresent ? 1 : 0 );

    if( numPresentPhases <= 1 )
    {
      integer presentPhaseIndex = liquidIndex;
      if( vapourPresent )
      {
        presentPhaseIndex = vapourIndex;
      }
      else if( aqueousPresent )
      {
        presentPhaseIndex = aqueousIndex;
      }
      assignSinglePhase( presentPhaseIndex, compFraction, phaseFraction, phaseCompFraction );
    }
    else if( liquidPresent && vapourPresent && !aqueousPresent )
    {
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdP( m_numComponents );
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdT( m_numComponents );
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        dKdP[ic] = dK( Deriv::dP, vapourKIndex, ic );
        dKdT[ic] = dK( Deriv::dT, vapourKIndex, ic );
      }
      computeTwoPhaseFromKValues( liquidIndex, vapourIndex, aqueousIndex,
                                  kValues[vapourKIndex].toSliceConst(),
                                  dKdP.toSliceConst(),
                                  dKdT.toSliceConst(),
                                  compFraction,
                                  vapourFraction,
                                  phaseFraction,
                                  phaseCompFraction );
    }
    else if( liquidPresent && aqueousPresent && !vapourPresent )
    {
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdP( m_numComponents );
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdT( m_numComponents );
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        dKdP[ic] = dK( Deriv::dP, aqueousKIndex, ic );
        dKdT[ic] = dK( Deriv::dT, aqueousKIndex, ic );
      }
      computeTwoPhaseFromKValues( liquidIndex, aqueousIndex, vapourIndex,
                                  kValues[aqueousKIndex].toSliceConst(),
                                  dKdP.toSliceConst(),
                                  dKdT.toSliceConst(),
                                  compFraction,
                                  aqueousFraction,
                                  phaseFraction,
                                  phaseCompFraction );
    }
    else if( vapourPresent && aqueousPresent && !liquidPresent )
    {
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > kGasWater( m_numComponents );
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdP( m_numComponents );
      stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdT( m_numComponents );
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const kW = LvArray::math::max( kValues( aqueousKIndex, ic ), MultiFluidConstants::epsilon );
        kGasWater[ic] = kValues( vapourKIndex, ic ) / kW;
        dKdP[ic] = ( dK( Deriv::dP, vapourKIndex, ic ) * kW - kValues( vapourKIndex, ic ) * dK( Deriv::dP, aqueousKIndex, ic ) ) / ( kW * kW );
        dKdT[ic] = ( dK( Deriv::dT, vapourKIndex, ic ) * kW - kValues( vapourKIndex, ic ) * dK( Deriv::dT, aqueousKIndex, ic ) ) / ( kW * kW );
      }
      computeTwoPhaseFromKValues( aqueousIndex, vapourIndex, liquidIndex,
                                  kGasWater.toSliceConst(),
                                  dKdP.toSliceConst(),
                                  dKdT.toSliceConst(),
                                  compFraction,
                                  vapourFraction,
                                  phaseFraction,
                                  phaseCompFraction );
    }
    else
    {
      // Three-phase region: implicit function theorem on (f_V, f_W)
      real64 fV = 0.0;
      real64 fW = 0.0;
      real64 j11 = 0.0;
      real64 j12 = 0.0;
      real64 j21 = 0.0;
      real64 j22 = 0.0;
      RachfordRice::evaluateThreePhaseResidualAndJacobian( kValues[vapourKIndex].toSliceConst(),
                                                           kValues[aqueousKIndex].toSliceConst(),
                                                           compFraction,
                                                           presentComponents,
                                                           vapourFraction,
                                                           aqueousFraction,
                                                           fV, fW, j11, j12, j21, j22 );

      // IFT requires a true root. If the flash left a non-zero residual (e.g. early exit),
      // fall back to the best two-phase edge so derivatives stay consistent with values.
      real64 const flashResidual = LvArray::math::max( LvArray::math::abs( fV ), LvArray::math::abs( fW ) );
      if( flashResidual > RachfordRice::newtonTolerance )
      {
        real64 const vOilGas = RachfordRice::solve( kValues[vapourKIndex].toSliceConst(),
                                                    compFraction, presentComponents );
        real64 const wOilWater = RachfordRice::solve( kValues[aqueousKIndex].toSliceConst(),
                                                      compFraction, presentComponents );
        stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > kGasWater( m_numComponents );
        for( integer ic = 0; ic < m_numComponents; ++ic )
        {
          kGasWater[ic] = kValues( vapourKIndex, ic ) /
                          LvArray::math::max( kValues( aqueousKIndex, ic ), MultiFluidConstants::epsilon );
        }
        real64 const vGasWater = RachfordRice::solve( kGasWater.toSliceConst(),
                                                      compFraction, presentComponents );

        real64 const fOilGas = RachfordRice::evaluateThreePhaseObjective(
          kValues[vapourKIndex].toSliceConst(), kValues[aqueousKIndex].toSliceConst(),
          compFraction, presentComponents, vOilGas, 0.0 );
        real64 const fOilWater = RachfordRice::evaluateThreePhaseObjective(
          kValues[vapourKIndex].toSliceConst(), kValues[aqueousKIndex].toSliceConst(),
          compFraction, presentComponents, 0.0, wOilWater );
        real64 const fGasWater = RachfordRice::evaluateThreePhaseObjective(
          kValues[vapourKIndex].toSliceConst(), kValues[aqueousKIndex].toSliceConst(),
          compFraction, presentComponents, vGasWater, 1.0 - vGasWater );

        stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdP( m_numComponents );
        stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > dKdT( m_numComponents );
        if( fOilGas <= fOilWater && fOilGas <= fGasWater )
        {
          for( integer ic = 0; ic < m_numComponents; ++ic )
          {
            dKdP[ic] = dK( Deriv::dP, vapourKIndex, ic );
            dKdT[ic] = dK( Deriv::dT, vapourKIndex, ic );
          }
          computeTwoPhaseFromKValues( liquidIndex, vapourIndex, aqueousIndex,
                                      kValues[vapourKIndex].toSliceConst(),
                                      dKdP.toSliceConst(), dKdT.toSliceConst(),
                                      compFraction, vOilGas, phaseFraction, phaseCompFraction );
        }
        else if( fOilWater <= fGasWater )
        {
          for( integer ic = 0; ic < m_numComponents; ++ic )
          {
            dKdP[ic] = dK( Deriv::dP, aqueousKIndex, ic );
            dKdT[ic] = dK( Deriv::dT, aqueousKIndex, ic );
          }
          computeTwoPhaseFromKValues( liquidIndex, aqueousIndex, vapourIndex,
                                      kValues[aqueousKIndex].toSliceConst(),
                                      dKdP.toSliceConst(), dKdT.toSliceConst(),
                                      compFraction, wOilWater, phaseFraction, phaseCompFraction );
        }
        else
        {
          for( integer ic = 0; ic < m_numComponents; ++ic )
          {
            real64 const kW = LvArray::math::max( kValues( aqueousKIndex, ic ), MultiFluidConstants::epsilon );
            dKdP[ic] = ( dK( Deriv::dP, vapourKIndex, ic ) * kW
                         - kValues( vapourKIndex, ic ) * dK( Deriv::dP, aqueousKIndex, ic ) ) / ( kW * kW );
            dKdT[ic] = ( dK( Deriv::dT, vapourKIndex, ic ) * kW
                         - kValues( vapourKIndex, ic ) * dK( Deriv::dT, aqueousKIndex, ic ) ) / ( kW * kW );
          }
          computeTwoPhaseFromKValues( aqueousIndex, vapourIndex, liquidIndex,
                                      kGasWater.toSliceConst(),
                                      dKdP.toSliceConst(), dKdT.toSliceConst(),
                                      compFraction, vGasWater, phaseFraction, phaseCompFraction );
        }
      }
      else
      {
        real64 const det = j11 * j22 - j12 * j21;
        GEOS_ERROR_IF( LvArray::math::abs( det ) < MultiFluidConstants::epsilon,
                       "Failed to calculate derivatives for the three-phase Rachford-Rice equation." );
        real64 const invDet = 1.0 / det;

        phaseFraction.value[vapourIndex] = vapourFraction;
        phaseFraction.value[aqueousIndex] = aqueousFraction;
        phaseFraction.value[liquidIndex] = liquidFraction;

        for( integer idof = 0; idof < numDofs; ++idof )
        {
          real64 rhsV = 0.0;
          real64 rhsW = 0.0;
          if( idof == Deriv::dP || idof == Deriv::dT )
          {
            for( integer ic = 0; ic < m_numComponents; ++ic )
            {
              real64 const dV = kValues( vapourKIndex, ic ) - 1.0;
              real64 const dW = kValues( aqueousKIndex, ic ) - 1.0;
              real64 const D = LvArray::math::max( 1.0 + vapourFraction * dV + aqueousFraction * dW, MultiFluidConstants::epsilon );
              real64 const invD = 1.0 / D;
              real64 const dKiV = dK( idof, vapourKIndex, ic );
              real64 const dKiW = dK( idof, aqueousKIndex, ic );
              real64 const dD = vapourFraction * dKiV + aqueousFraction * dKiW;
              rhsV += compFraction[ic] * ( dKiV * invD - dV * dD * invD * invD );
              rhsW += compFraction[ic] * ( dKiW * invD - dW * dD * invD * invD );
            }
          }
          else
          {
            integer const jc = idof - Deriv::dC;
            real64 const dV = kValues( vapourKIndex, jc ) - 1.0;
            real64 const dW = kValues( aqueousKIndex, jc ) - 1.0;
            real64 const D = LvArray::math::max( 1.0 + vapourFraction * dV + aqueousFraction * dW, MultiFluidConstants::epsilon );
            rhsV = dV / D;
            rhsW = dW / D;
          }

          real64 const dVap = ( -j22 * rhsV + j12 * rhsW ) * invDet;
          real64 const dAq = ( j21 * rhsV - j11 * rhsW ) * invDet;
          phaseFraction.derivs( vapourIndex, idof ) = dVap;
          phaseFraction.derivs( aqueousIndex, idof ) = dAq;
          phaseFraction.derivs( liquidIndex, idof ) = -dVap - dAq;
        }

        for( integer ic = 0; ic < m_numComponents; ++ic )
        {
          real64 const dV = kValues( vapourKIndex, ic ) - 1.0;
          real64 const dW = kValues( aqueousKIndex, ic ) - 1.0;
          real64 const D = LvArray::math::max( 1.0 + vapourFraction * dV + aqueousFraction * dW, MultiFluidConstants::epsilon );
          real64 const invD = 1.0 / D;
          real64 const xi = compFraction[ic] * invD;
          real64 const yi = xi * kValues( vapourKIndex, ic );
          real64 const wi = xi * kValues( aqueousKIndex, ic );
          phaseCompFraction.value( liquidIndex, ic ) = xi;
          phaseCompFraction.value( vapourIndex, ic ) = yi;
          phaseCompFraction.value( aqueousIndex, ic ) = wi;

          for( integer idof = 0; idof < numDofs; ++idof )
          {
            real64 const dVap = phaseFraction.derivs( vapourIndex, idof );
            real64 const dAq = phaseFraction.derivs( aqueousIndex, idof );
            real64 dKiV = 0.0;
            real64 dKiW = 0.0;
            real64 dz = 0.0;
            if( idof == Deriv::dP || idof == Deriv::dT )
            {
              dKiV = dK( idof, vapourKIndex, ic );
              dKiW = dK( idof, aqueousKIndex, ic );
            }
            else if( idof == Deriv::dC + ic )
            {
              dz = 1.0;
            }
            real64 const dD = dVap * dV + vapourFraction * dKiV + dAq * dW + aqueousFraction * dKiW;
            real64 const dxi = ( dz - xi * dD ) * invD;
            phaseCompFraction.derivs( liquidIndex, ic, idof ) = dxi;
            phaseCompFraction.derivs( vapourIndex, ic, idof ) = dKiV * xi + kValues( vapourKIndex, ic ) * dxi;
            phaseCompFraction.derivs( aqueousIndex, ic, idof ) = dKiW * xi + kValues( aqueousKIndex, ic ) * dxi;
          }
        }
      }
    }
  }
}

template< integer NUM_PHASE >
template< int USD1, int USD2, int USD3, int USD4 >
GEOS_HOST_DEVICE
void
KValueFlashModelUpdate< NUM_PHASE >::computeTwoPhaseFromKValues( integer const primaryIndex,
                                                                 integer const secondaryIndex,
                                                                 integer const absentIndex,
                                                                 arraySlice1d< real64 const, USD1 > const & kValues,
                                                                 arraySlice1d< real64 const, USD2 > const & dKdP,
                                                                 arraySlice1d< real64 const, USD3 > const & dKdT,
                                                                 arraySlice1d< real64 const, USD4 > const & compFraction,
                                                                 real64 secondaryFraction,
                                                                 PhaseProp::SliceType const phaseFraction,
                                                                 PhaseComp::SliceType const phaseCompFraction ) const
{
  integer const numDofs = 2 + m_numComponents;

  real64 denominator = 0.0;
  real64 pressureNumerator = 0.0;
  real64 temperatureNumerator = 0.0;
  for( integer ic = 0; ic < m_numComponents; ++ic )
  {
    real64 const k = kValues[ic] - 1.0;
    real64 const r = 1.0 / ( 1.0 + secondaryFraction * k );
    pressureNumerator += compFraction[ic] * dKdP[ic] * r * r;
    temperatureNumerator += compFraction[ic] * dKdT[ic] * r * r;
    denominator += compFraction[ic] * k * k * r * r;
    phaseFraction.derivs( secondaryIndex, Deriv::dC + ic ) = k * r;
  }

  if( denominator < MultiFluidConstants::epsilon )
  {
    assignSinglePhase( secondaryFraction > 0.5 ? secondaryIndex : primaryIndex,
                       compFraction, phaseFraction, phaseCompFraction );
    return;
  }

  real64 const invDenominator = 1.0 / denominator;
  phaseFraction.derivs( secondaryIndex, Deriv::dP ) = pressureNumerator * invDenominator;
  phaseFraction.derivs( secondaryIndex, Deriv::dT ) = temperatureNumerator * invDenominator;
  for( integer ic = 0; ic < m_numComponents; ++ic )
  {
    phaseFraction.derivs( secondaryIndex, Deriv::dC + ic ) *= invDenominator;
  }

  for( integer ic = 0; ic < m_numComponents; ++ic )
  {
    real64 const k = kValues[ic] - 1.0;
    real64 const r = 1.0 / ( 1.0 + secondaryFraction * k );
    real64 const xi = compFraction[ic] * r;
    real64 const yi = xi * kValues[ic];
    phaseCompFraction.value( primaryIndex, ic ) = xi;
    phaseCompFraction.value( secondaryIndex, ic ) = yi;

    for( integer const idof : {Deriv::dP, Deriv::dT} )
    {
      real64 const dV = phaseFraction.derivs( secondaryIndex, idof );
      real64 const dKi = ( idof == Deriv::dP ) ? dKdP[ic] : dKdT[ic];
      real64 const dxi = r * ( -secondaryFraction * xi * dKi + ( xi - yi ) * dV );
      phaseCompFraction.derivs( primaryIndex, ic, idof ) = dxi;
      phaseCompFraction.derivs( secondaryIndex, ic, idof ) = xi * dKi + dxi * kValues[ic];
    }
    for( integer jc = 0; jc < m_numComponents; ++jc )
    {
      integer const idof = Deriv::dC + jc;
      real64 const dV = phaseFraction.derivs( secondaryIndex, idof );
      real64 const dz = ( jc == ic ) ? 1.0 : 0.0;
      real64 const dxi = r * ( dz + ( xi - yi ) * dV );
      phaseCompFraction.derivs( primaryIndex, ic, idof ) = dxi;
      phaseCompFraction.derivs( secondaryIndex, ic, idof ) = dxi * kValues[ic];
    }
  }

  if( secondaryFraction < MultiFluidConstants::epsilon || 1.0 - secondaryFraction < MultiFluidConstants::epsilon )
  {
    secondaryFraction = LvArray::math::min( LvArray::math::max( secondaryFraction, 0.0 ), 1.0 );
    integer const singlePhaseIndex = secondaryFraction < MultiFluidConstants::epsilon ? primaryIndex : secondaryIndex;
    assignSinglePhase( singlePhaseIndex, compFraction, phaseFraction, phaseCompFraction );
    return;
  }

  phaseFraction.value[secondaryIndex] = secondaryFraction;
  phaseFraction.value[primaryIndex] = 1.0 - secondaryFraction;
  for( integer ic = 0; ic < numDofs; ic++ )
  {
    phaseFraction.derivs[primaryIndex][ic] = -phaseFraction.derivs[secondaryIndex][ic];
  }

  if( absentIndex >= 0 )
  {
    phaseFraction.value[absentIndex] = 0.0;
    for( integer ic = 0; ic < numDofs; ++ic )
    {
      phaseFraction.derivs( absentIndex, ic ) = 0.0;
    }
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      phaseCompFraction.value( absentIndex, ic ) = compFraction[ic];
      for( integer idof = 0; idof < numDofs; ++idof )
      {
        phaseCompFraction.derivs( absentIndex, ic, idof ) = 0.0;
      }
      phaseCompFraction.derivs( absentIndex, ic, Deriv::dC + ic ) = 1.0;
    }
  }
}

template< integer NUM_PHASE >
template< int USD1 >
GEOS_HOST_DEVICE
void
KValueFlashModelUpdate< NUM_PHASE >::assignSinglePhase( integer const presentPhaseIndex,
                                                        arraySlice1d< real64 const, USD1 > const & compFraction,
                                                        PhaseProp::SliceType const phaseFraction,
                                                        PhaseComp::SliceType const phaseCompFraction ) const
{
  LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );
  for( integer ip = 0; ip < NUM_PHASE; ++ip )
  {
    phaseFraction.value[ip] = ( ip == presentPhaseIndex ) ? 1.0 : 0.0;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      phaseCompFraction.value( ip, ic ) = compFraction[ic];
      phaseCompFraction.derivs( ip, ic, Deriv::dC + ic ) = 1.0;
    }
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
