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
 * @file ThermalCompositionalMultiphaseReservoirFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief ThermalCompositionalMultiphaseReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *               0 = reservoir pressure
 *               1 = reservoir density (component 1 )
 *             ... = reservoir densities (# components , NC)
 *          NC + 1 = reservoir temperature
 *
 * numResLabels - 1 = reservoir temperature
 *     numResLabels = well pressure
 * numResLabels + 1 = well density
 *              ... = ... (well densities)
 *   isothermal well (WellManager isThermal="0"):
 *     numResLabels + numWellLabels - 2 = last well density
 *     numResLabels + numWellLabels - 1 = well rate
 *   thermal well (WellManager isThermal="1"):
 *     numResLabels + numWellLabels - 3 = last well density
 *     numResLabels + numWellLabels - 2 = well rate
 *     numResLabels + numWellLabels - 1 = well temperature
 *
 * Level 0 eliminates the entire well block, so isothermal wells (no well
 * temperature DOF) and thermal wells are both supported. Reservoir temperature
 * remains a coarse unknown (CPTR).
 *
 * 3-level MGR reduction strategy
 *   - 1st level: eliminate the well block
 *   - 2nd level: eliminate the reservoir density associated with the volume constraint
 *   - 3rd level: eliminate the remaining the reservoir densities
 *   - The coarse grid (pressure and temperature system) is solved with BoomerAMG with numFunctions==2.
 *
 */

class ThermalCompositionalMultiphaseReservoirFVM : public MGRStrategyBase< 3 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit ThermalCompositionalMultiphaseReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );

    // Level 0: eliminate the well block (isothermal or thermal well DOF layout)
    m_labels[0].resize( numResLabels );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );

    // Level 1: eliminate last density which corresponds to the volume constraint equation
    m_labels[1].resize( numResLabels - 2 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 0 );
    m_labels[1].push_back( numResLabels-1 ); // keep temperature
    // Level 2: eliminate the other densities
    m_labels[2].push_back( 0 ); // keep pressure
    m_labels[2].push_back( numResLabels-1 ); // keep temperature

    setupLabels();

    // Level 0: eliminate well block (same recipe as CompositionalMultiphaseReservoirFVM).
    // gsElimWInverse is required: Jacobi F-relax leaves zero-diagonal well rows and
    // hypre then fails ("Zero diagonal found"). areWellsShut switches to Jacobi below.
    m_levelFRelaxType[0]          = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    m_levelFRelaxType[1]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::jacobi; // Diagonal scaling (Jacobi)
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin; // Standard Galerkin
    // blockGaussSeidel + col-lumped restriction both segfault in hypre 3.2.0 on
    // this thermal-reservoir + wells system (isolated by A/B runs on HE_3Comp).
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;

    m_levelFRelaxType[2]          = MGRFRelaxationType::none;
    m_levelInterpType[2]          = MGRInterpolationType::injection;
    // Prefer injection over partialColLumped/blockColLumped: both col-lumped
    // options segfault in hypre_MGRSetup after the well-elimination level.
    m_levelRestrictType[2]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[2] = 1;
  }

  /**
   * @brief Apply MGR parameters that affect the reduction hierarchy.
   * @param mgrParams MGR configuration parameters
   */
  void configure( LinearSolverParameters::MGR const & mgrParams )
  {
    // if the wells are shut, using Gaussian elimination as F-relaxation for the well block is an overkill
    // (and can be singular); in that case, we just use Jacobi
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[0] = MGRFRelaxationType::jacobi;
    }
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams parameters for the configuration of the MGR recipe
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    configure( mgrParams );
    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure/temperature reduced system
    setPressureTemperatureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_*/
