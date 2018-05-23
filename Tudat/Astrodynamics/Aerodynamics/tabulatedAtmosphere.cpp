/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/make_shared.hpp>
#include <stdlib.h>
#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"

namespace tudat
{
namespace aerodynamics
{

//! Initialize atmosphere table reader.
void TabulatedAtmosphere::initialize( const std::string& atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    Eigen::MatrixXd containerOfAtmosphereTableFileData
            = input_output::readMatrixFromFile( atmosphereTableFile_, " \t", "%" );

    // Check whether data is present in the file.
    if ( containerOfAtmosphereTableFileData.rows( ) < 1
         || containerOfAtmosphereTableFileData.cols( ) < 1 )
    {
        std::string errorMessage = "The atmosphere table file " + atmosphereTableFile_ + " is empty";
        throw std::runtime_error( errorMessage );
    }

    // Initialize vectors.
    altitudeData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    densityData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    pressureData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    temperatureData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    specificHeatRatioData_.resize( containerOfAtmosphereTableFileData.rows( ) );

    int densityIndex = 0;
    int pressureIndex = 0;
    int temperatureIndex = 0;
    int specificHeatRatioIndex = 0;
    int gasConstantIndex = 0;
    containsSpecificHeatRatio_ = false;
    containsGasConstant_ = false;

    for( unsigned int i = 0; i < dependentVariables_.size( ); i++ )
    {
        switch( dependentVariables_[ i ] )
        {
        case density_dependent_atmosphere:
            densityIndex = i + 1;
            break;
        case pressure_dependent_atmosphere:
            pressureIndex = i + 1;
            break;
        case temperature_dependent_atmosphere:
            temperatureIndex = i + 1;
            break;
        case specific_heat_ratio_dependent_atmosphere:
            containsSpecificHeatRatio_ = true;
            specificHeatRatioIndex = i + 1;
            break;
        case gas_constant_dependent_atmosphere:
            containsGasConstant_ = true;
            gasConstantIndex = i + 1;
            break;
        default:
            std::string errorMessage = "Error, dependent variable " +
                    std::to_string( dependentVariables_[i] ) +
                    " not found in tabulated atmosphere";
            throw std::runtime_error( errorMessage );
        }
    }

    if( densityIndex == 0 || pressureIndex == 0 || temperatureIndex == 0 )
    {
        throw std::runtime_error(
                    "Error, tabulated atmosphere must be initialized with at least temperature, pressure and density" );
    }

    // Loop through all the strings stored in the container and store the data
    // in the right Eigen::VectorXd.
    for ( int i = 0; i < containerOfAtmosphereTableFileData.rows( ); i++  )
    {
        altitudeData_[ i ] = containerOfAtmosphereTableFileData( i, 0 );
        densityData_[ i ] = containerOfAtmosphereTableFileData( i, 1 );
        pressureData_[ i ] = containerOfAtmosphereTableFileData( i, 2 );
        temperatureData_[ i ] = containerOfAtmosphereTableFileData( i, 3 );
        specificHeatRatioData_[ i ] = containerOfAtmosphereTableFileData( i, 4 );
    }


    using namespace interpolators;

    cubicSplineInterpolationForDensity_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, densityData_ );
    cubicSplineInterpolationForPressure_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, pressureData_ );
    cubicSplineInterpolationForTemperature_
            = boost::make_shared< CubicSplineInterpolatorDouble >(
                altitudeData_, temperatureData_ );
    cubicSplineInterpolationForHeatRatio_
            = boost::make_shared< CubicSplineInterpolatorDouble >(
                altitudeData_, specificHeatRatioData_ );
}

//! Initialize species table reader.
void TabulatedAtmosphere::initializeSpecies( const std::string& speciesTableFile )
{
    // Locally store the atmosphere table file name.
    speciesTableFile_ = speciesTableFile;

    Eigen::MatrixXd containerOfSpeciesTableFileData
            = input_output::readMatrixFromFile( speciesTableFile_, " \t", "%" );

    // Check whether data is present in the file.
    if ( containerOfSpeciesTableFileData.rows( ) < 1
         || containerOfSpeciesTableFileData.cols( ) < 1 )
    {
        std::string errorMessage = "The species table file " + speciesTableFile_ + " is empty";
        throw std::runtime_error( errorMessage );
    }

    // Initialize vectors.
    CO2Data_.resize( containerOfSpeciesTableFileData.rows( ) );
    N2Data_.resize( containerOfSpeciesTableFileData.rows( ) );
    ArData_.resize( containerOfSpeciesTableFileData.rows( ) );
    COData_.resize( containerOfSpeciesTableFileData.rows( ) );
    OData_.resize( containerOfSpeciesTableFileData.rows( ) );
    O2Data_.resize( containerOfSpeciesTableFileData.rows( ) );
    O3Data_.resize( containerOfSpeciesTableFileData.rows( ) );
    HData_.resize( containerOfSpeciesTableFileData.rows( ) );
    H2Data_.resize( containerOfSpeciesTableFileData.rows( ) );

    // Loop through all the strings stored in the container and store the data
    // in the right Eigen::VectorXd.
    for ( int i = 0; i < containerOfSpeciesTableFileData.rows( ); i++  )
    {
        CO2Data_[ i ] = containerOfSpeciesTableFileData( i, 0 );
        N2Data_[ i ] = containerOfSpeciesTableFileData( i, 1 );
        ArData_[ i ] = containerOfSpeciesTableFileData( i, 2 );
        COData_[ i ] = containerOfSpeciesTableFileData( i, 3 );
        OData_[ i ] = containerOfSpeciesTableFileData( i, 4 );
        O2Data_[ i ] = containerOfSpeciesTableFileData( i, 5 );
        O3Data_[ i ] = containerOfSpeciesTableFileData( i, 6 );
        HData_[ i ] = containerOfSpeciesTableFileData( i, 7 );
        H2Data_[ i ] = containerOfSpeciesTableFileData( i, 8 );
    }

    using namespace interpolators;

    cubicSplineInterpolationForCO2_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, CO2Data_ );
    cubicSplineInterpolationForN2_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, N2Data_ );
    cubicSplineInterpolationForAr_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, ArData_ );
    cubicSplineInterpolationForCO_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, COData_ );
    cubicSplineInterpolationForO_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, OData_ );
    cubicSplineInterpolationForO2_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, O2Data_ );
    cubicSplineInterpolationForO3_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, O3Data_ );
    cubicSplineInterpolationForH_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, HData_ );
    cubicSplineInterpolationForH2_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, H2Data_ );
}

void TabulatedAtmosphere::computeProperties(const double altitude)
{
    double molGasConst = tudat::physical_constants::MOLAR_GAS_CONSTANT;
    const double AVOGADRO_CONSTANT = 6.02214129E23;
    double temp = cubicSplineInterpolationForTemperature_->interpolate( altitude );
    double pres = cubicSplineInterpolationForPressure_->interpolate( altitude );
    double conversion = AVOGADRO_CONSTANT*pres/(molGasConst*temp);

    double CO2 = cubicSplineInterpolationForCO2_->interpolate( altitude );
    double N2 = cubicSplineInterpolationForN2_->interpolate( altitude );
    double Ar = cubicSplineInterpolationForAr_->interpolate( altitude );
    double CO = cubicSplineInterpolationForCO_->interpolate( altitude );
    double O = cubicSplineInterpolationForO_->interpolate( altitude );
    double O2 = cubicSplineInterpolationForO2_->interpolate( altitude );
    double O3 = cubicSplineInterpolationForO3_->interpolate( altitude );
    double H = cubicSplineInterpolationForH_->interpolate( altitude );
    double H2 = cubicSplineInterpolationForH2_->interpolate( altitude );

    // Get number densities
    numberDensities_.resize(9);
    numberDensities_[0] = conversion*CO2 ;
    numberDensities_[1] = conversion*N2 ;
    numberDensities_[2] = conversion*Ar ;
    numberDensities_[3] = conversion*CO ;
    numberDensities_[4] = conversion*O ;
    numberDensities_[5] = conversion*O2 ;
    numberDensities_[6] = conversion*O3 ;
    numberDensities_[7] = conversion*H ;
    numberDensities_[8] = conversion*H2 ;

    // Get average number density
    double sumOfNumberDensity = 0.0 ;
    for( unsigned int i = 0 ; i < numberDensities_.size( ) ; i++)
    {
        sumOfNumberDensity += numberDensities_[ i ];
    }
    averageNumberDensity_ = sumOfNumberDensity / double( numberDensities_.size( ) );

    // Mean molar mass (Thermodynamics an Engineering Approach, Michael A. Boles)
    meanMolarMass_ = numberDensities_[0] * gasComponentProperties_.molarMassCarbonDiOxide;
    meanMolarMass_ += numberDensities_[1] * gasComponentProperties_.molarMassNitrogen;
    meanMolarMass_ += numberDensities_[2] * gasComponentProperties_.molarMassArgon;
    meanMolarMass_ += numberDensities_[3] * gasComponentProperties_.molarMassCarbonMonoOxide;
    meanMolarMass_ += numberDensities_[4] * gasComponentProperties_.molarMassAtomicOxygen;
    meanMolarMass_ += numberDensities_[5] * gasComponentProperties_.molarMassOxygen;
    meanMolarMass_ += numberDensities_[6] * gasComponentProperties_.molarMassOzon;
    meanMolarMass_ += numberDensities_[7] * gasComponentProperties_.molarMassAtomicHydrogen;
    meanMolarMass_ += numberDensities_[8] * gasComponentProperties_.molarMassHydrogen;
    meanMolarMass_ = meanMolarMass_ / sumOfNumberDensity ;


    // Collision diameter
    weightedAverageCollisionDiameter_ = numberDensities_[0]* gasComponentProperties_.diameterCarbonDiOxide ;
    weightedAverageCollisionDiameter_ += numberDensities_[1]* gasComponentProperties_.diameterNitrogen ;
    weightedAverageCollisionDiameter_ += numberDensities_[2]* gasComponentProperties_.diameterArgon ;
    weightedAverageCollisionDiameter_ += numberDensities_[3]* gasComponentProperties_.diameterCarbonMonoOxide ;
    weightedAverageCollisionDiameter_ += numberDensities_[4]* gasComponentProperties_.diameterAtomicOxygen ;
    weightedAverageCollisionDiameter_ += numberDensities_[5]* gasComponentProperties_.diameterOxygen ;
    weightedAverageCollisionDiameter_ += numberDensities_[6]* gasComponentProperties_.diameterOzon ;
    weightedAverageCollisionDiameter_ += numberDensities_[7]* gasComponentProperties_.diameterAtomicHydrogen ;
    weightedAverageCollisionDiameter_ += numberDensities_[8]* gasComponentProperties_.diameterHydrogen ;
    weightedAverageCollisionDiameter_ = weightedAverageCollisionDiameter_ / sumOfNumberDensity;

    ratioOfSpecificHeats_ = cubicSplineInterpolationForHeatRatio_->interpolate(altitude);

}


} // namespace aerodynamics
} // namespace tudat
