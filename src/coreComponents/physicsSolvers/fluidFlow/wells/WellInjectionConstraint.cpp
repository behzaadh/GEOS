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

/*
 * @file WellInjectionConstraint.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellInjectionConstraint.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

#include "WellMassRateConstraint.hpp"
#include "WellPhaseVolumeRateConstraint.hpp"
#include "WellVolumeRateConstraint.hpp"

namespace geos
{
using namespace dataRepository;

template< typename ConstraintRateType >
InjectionConstraint< ConstraintRateType >::InjectionConstraint( string const & name, Group * const parent )
  : ConstraintRateType( name, parent )
{
  // set rate sign for injectors (base class member)
  this->m_rateSign = 1.0;
  classtype::registerWrapper( injectionStreamKey::injectionStreamString(), &m_injectionStream ).
    setDefaultValue( -1 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Global component densities of the injection stream [moles/m^3 or kg/m^3]" );

  InjectionConstraint< ConstraintRateType >::registerWrapper( injectionStreamKey::injectionTemperatureString(), &m_injectionTemperature ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature of the injection stream [K]" );

  InjectionConstraint< ConstraintRateType >::registerWrapper( injectionStreamKey::injectionTemperatureTableNameString(), &m_injectionTemperatureTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the table of injection stream temperature [K] versus time. Replaces "
                    + string( injectionStreamKey::injectionTemperatureString() ) );
}
template< typename ConstraintRateType >
InjectionConstraint< ConstraintRateType >::~InjectionConstraint()
{}
template< typename ConstraintRateType >
void InjectionConstraint< ConstraintRateType >::postInputInitialization()
{
  // Validate value and table options
  ConstraintRateType::postInputInitialization();

  if( !m_injectionTemperatureTableName.empty() )
  {
    GEOS_THROW_IF( m_injectionTemperature >= 0,
                   this->getDataContext() << ": " << injectionStreamKey::injectionTemperatureString() << " and "
                                          << injectionStreamKey::injectionTemperatureTableNameString() << " cannot be specified together",
                   InputError );
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_injectionTemperatureTable = &(functionManager.getGroup< TableFunction const >( m_injectionTemperatureTableName ));
  }

// Validate the injection stream and temperature
  validateInjectionStream(  );

}

template< typename ConstraintRateType >
void InjectionConstraint< ConstraintRateType >::setNextDtFromTables( real64 const currentTime, real64 & nextDt ) const
{
  ConstraintRateType::setNextDtFromTables( currentTime, nextDt );
  WellConstraintBase::setNextDtFromTable( m_injectionTemperatureTable, currentTime, nextDt );
}
template< typename ConstraintRateType >
void InjectionConstraint< ConstraintRateType >::validateInjectionStream( )
{
  bool const hasTemperature = m_injectionTemperature >= 0 || m_injectionTemperatureTable != nullptr;
  GEOS_THROW_IF( (m_injectionStream.empty()  && hasTemperature) ||
                 (!m_injectionStream.empty() && !hasTemperature),
                 this->getName() << " "  <<  this->getDataContext() << ": Both "
                                 <<  injectionStreamKey::injectionStreamString() << " and " <<  injectionStreamKey::injectionTemperatureString()
                                 << " must be specified for multiphase simulations",
                 InputError );

  if( !m_injectionStream.empty())
  {
    real64 sum = 0.0;
    for( localIndex ic = 0; ic < m_injectionStream.size(); ++ic )
    {
      GEOS_ERROR_IF( m_injectionStream[ic] < 0.0 || m_injectionStream[ic] > 1.0,
                     classtype::getWrapperDataContext( injectionStreamKey::injectionStreamString() ) << ": Invalid injection stream" );
      sum += m_injectionStream[ic];
    }
    GEOS_THROW_IF( LvArray::math::abs( 1.0 - sum ) > std::numeric_limits< real64 >::epsilon(),
                   classtype::getWrapperDataContext( injectionStreamKey::injectionStreamString() ) << ": Invalid injection stream",
                   InputError );
  }
}

// Register concrete wrapper constraint types and instantiate templates.
//template class InjectionConstraint< LiquidRateConstraint >;
//using InjectionLiquidRateConstraint = InjectionConstraint< LiquidRateConstraint >;
//REGISTER_CATALOG_ENTRY( WellConstraintBase, InjectionLiquidRateConstraint, string const &, Group * const )

template class InjectionConstraint< MassRateConstraint >;
using InjectionMassRateConstraint = InjectionConstraint< MassRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, InjectionMassRateConstraint, string const &, Group * const )

template class InjectionConstraint< PhaseVolumeRateConstraint >;
using InjectionPhaseVolumeRateConstraint = InjectionConstraint< PhaseVolumeRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, InjectionPhaseVolumeRateConstraint, string const &, Group * const )

template class InjectionConstraint< VolumeRateConstraint >;
using InjectionVolumeRateConstraint = InjectionConstraint< VolumeRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, InjectionVolumeRateConstraint, string const &, Group * const )


}
