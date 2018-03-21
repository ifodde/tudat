/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FLIGHTCONDITIONS_H
#define TUDAT_FLIGHTCONDITIONS_H

#include <vector>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <math.h>

#include "Tudat/Astrodynamics/Aerodynamics/trimOrientation.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h"
#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00InputFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/SystemModels/vehicleSystems.h"
#include <Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>


namespace tudat
{

namespace aerodynamics
{

//! Class for calculating aerodynamic flight characteristics of a vehicle during numerical
//! integration.
/*!
 *  Class for calculating aerodynamic flight characteristics of a vehicle during numerical
 *  integration. Class is used to ensure that dependent variables such as density, altitude, etc.
 *  are only calculated once during each numerical integration step. The get functions of this class
 *  are linked to the various models in the code that subsequently require these values.
 */
class FlightConditions
{
private:

    //! List of variables that can be computed by flight condition
    enum FlightConditionVariables
    {
        altitude_flight_condition,
        density_flight_condition,
        pressure_flight_condition,
        temperature_flight_condition,
        latitude_flight_condition,
        longitude_flight_condition,
        mach_number_flight_condition,
        speed_of_sound_flight_condition,
        airspeed_flight_condition,
        geodetic_latitude_condition,
        dynamic_pressure_condition
    };

public:

    //! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param atmosphereModel Atmosphere model of atmosphere through which vehicle is flying
     *  \param shapeModel Model describing the shape of the body w.r.t. which the flight is taking place.
     *  \param aerodynamicCoefficientInterface Class from which the aerodynamic (force and moment)
     *  coefficients are retrieved
     *  \param aerodynamicAngleCalculator Object from which the aerodynamic/trajectory angles
     *  of the vehicle are calculated.
     *  \param controlSurfaceDeflectionFunction Function returning control surface deflection, with input the control
     *  surface identifier.
     */
    FlightConditions(const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
                      const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
                      const boost::shared_ptr< AerodynamicCoefficientInterface >
                      aerodynamicCoefficientInterface,
                      const boost::shared_ptr< system_models::VehicleSystems > vehicleSystem,
                      const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
                      aerodynamicAngleCalculator =
            boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >( ),
                      const boost::function< double( const std::string& )> controlSurfaceDeflectionFunction =
            boost::function< double( const std::string& )>( ));

    //! Function to update all flight conditions.
    /*!
     *  Function to update all flight conditions (altitude, density, force coefficients) to
     *  current state of vehicle and central body.
     *  \param currentTime Time to which conditions are to be updated.
     */
    void updateConditions( const double currentTime );

    //! Function to retrieve (and compute if necessary) the current altitude
    /*!
     * Function to retrieve (and compute if necessary) the current altitude
     * \return Current altitude
     */
    double getCurrentAltitude( )
    {
        if( scalarFlightConditions_.count( altitude_flight_condition ) == 0 )
        {
            computeAltitude( );
        }
        return scalarFlightConditions_.at( altitude_flight_condition );    }

    //! Function to retrieve (and compute if necessary) the current freestream density
    /*!
     * Function to retrieve (and compute if necessary) the current freestream density
     * \return Current freestream density
     */
    double getCurrentDensity( )
    {
        if( scalarFlightConditions_.count( density_flight_condition ) == 0 )
        {
            computeDensity( );
        }
        return scalarFlightConditions_.at( density_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current freestream temperature
    /*!
     * Function to retrieve (and compute if necessary) the current freestream temperature
     * \return Current freestream temperature
     */
    double getCurrentFreestreamTemperature( )
    {
        if( scalarFlightConditions_.count( temperature_flight_condition ) == 0 )
        {
            computeTemperature( );
        }
        return scalarFlightConditions_.at( temperature_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current freestream dynamic pressure
    /*!
     * Function to retrieve (and compute if necessary) the current freestream dynamic pressure
     * \return Current freestream dynamic pressure
     */
    double getCurrentDynamicPressure( )
    {
        if( scalarFlightConditions_.count( dynamic_pressure_condition ) == 0 )
        {
            computeDynamicPressure( );
        }
        return scalarFlightConditions_.at( dynamic_pressure_condition );
    }

    //! Function to retrieve (and compute if necessary) the current freestream pressure
    /*!
     * Function to retrieve (and compute if necessary) the current freestream pressure
     * \return Current freestream dynamic pressure
     */
    double getCurrentPressure( )
    {
        if( scalarFlightConditions_.count( pressure_flight_condition ) == 0 )
        {
            computeFreestreamPressure( );
        }
        return scalarFlightConditions_.at( pressure_flight_condition );
    }

    /*!
     * Function to retrieve (and compute if necessary) the current airspeed
     * \return Current airspeed
     */
    double getCurrentAirspeed( )
    {
        if( scalarFlightConditions_.count( airspeed_flight_condition ) == 0 )
        {
            computeAirspeed( );
        }
        return scalarFlightConditions_.at( airspeed_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current speed of sound
    /*!
     * Function to retrieve (and compute if necessary) the current speed of sound
     * \return Current speed of sound
     */
    double getCurrentSpeedOfSound( )
    {
        if( scalarFlightConditions_.count( speed_of_sound_flight_condition ) == 0 )
        {
            computeSpeedOfSound( );
        }
        return scalarFlightConditions_.at( speed_of_sound_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current Mach number
    /*!
     * Function to retrieve (and compute if necessary) the current Mach number
     * \return Current Mach number
     */
    double getCurrentMachNumber( )
    {
        if( scalarFlightConditions_.count( mach_number_flight_condition ) == 0 )
        {
            computeMachNumber( );
        }
        return scalarFlightConditions_.at( mach_number_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current geodetic latitude
    /*!
     * Function to retrieve (and compute if necessary) the current geodetic latitude
     * \return Current geodetic latitude
     */
    double getCurrentGeodeticLatitude( )
    {
        if( scalarFlightConditions_.count( geodetic_latitude_condition ) == 0 )
        {
            computeGeodeticLatitude( );
        }
        return scalarFlightConditions_.at( geodetic_latitude_condition );
    }

    //! Function to return the current time of the FlightConditions
    /*!
     *  Function to return the current time of the FlightConditions.
     *  \return Current time of the FlightConditions
     */
    double getCurrentTime( )
    {
        return currentTime_;
    }

    double getCurrentKnudsenNumber( )
    {
        const double AVOGADRO_CONSTANT = 6.02214129E23;

        double characteristicLength = vehicleSystem_->getRefLength();


        boost::shared_ptr< aerodynamics::AtmosphereModel > Atmosphere = atmosphereModel_;

        double molMass = Atmosphere->getMeanMolarMass(scalarFlightConditions_.at( altitude_flight_condition ),
                                                            scalarFlightConditions_.at( longitude_flight_condition ),
                                                            scalarFlightConditions_.at( latitude_flight_condition ), currentTime_);
        double colDiameter = Atmosphere->getWeightedAverageCollisionDiameter(scalarFlightConditions_.at( altitude_flight_condition ),
                                                                                   scalarFlightConditions_.at( longitude_flight_condition ),
                                                                                   scalarFlightConditions_.at( latitude_flight_condition ), currentTime_);

        double Density = getCurrentDensity();
        double particleMass = (molMass / 1000.0) / AVOGADRO_CONSTANT;
        return particleMass / (std::sqrt(2.0) *mathematical_constants::PI *
                            std::pow(colDiameter,2.0) * Density) / characteristicLength;
    }

    double getRatioOfSpecificHeats( )
    {
        boost::shared_ptr< aerodynamics::AtmosphereModel > Atmosphere = atmosphereModel_;
        return Atmosphere->getRatioOfSpecificHeats(scalarFlightConditions_.at( altitude_flight_condition ),
                                                      scalarFlightConditions_.at( longitude_flight_condition ),
                                                      scalarFlightConditions_.at( latitude_flight_condition ), currentTime_);
    }

    double getSuttonGravesConvHeatFlux( )
    {
        double velocity = getCurrentAirspeed();
        double density = getCurrentDensity();
        double noseRadius = vehicleSystem_->getNoseRadius();
        double k = 1.83e-4;

        double convectiveHeatFlux = k*std::sqrt(density/noseRadius)*std::pow(velocity,3.0);

        return convectiveHeatFlux;

    }

    double getTauberRadHeatFlux( )
    {
        std::vector<double> velocityFunction(37);
        std::vector<double> heatingFunction(37);
        velocityFunction[0]  = 6000;   heatingFunction[0]  = 0.2;
        velocityFunction[1]  = 6150;   heatingFunction[1]  = 1.0;
        velocityFunction[2]  = 6300;   heatingFunction[2]  = 1.95;
        velocityFunction[3]  = 6500;   heatingFunction[3]  = 3.42;
        velocityFunction[4]  = 6700;   heatingFunction[4]  = 5.1;
        velocityFunction[5]  = 6900;   heatingFunction[5]  = 7.1;
        velocityFunction[6]  = 7000;   heatingFunction[6]  = 8.1;
        velocityFunction[7]  = 7200;   heatingFunction[7]  = 10.2;
        velocityFunction[8]  = 7400;   heatingFunction[8]  = 12.5;
        velocityFunction[9]  = 7600;   heatingFunction[9]  = 14.8;
        velocityFunction[10] = 7800;   heatingFunction[10] = 17.1;
        velocityFunction[11] = 8000;   heatingFunction[11] = 19.2;
        velocityFunction[12] = 8200;   heatingFunction[12] = 21.4;
        velocityFunction[13] = 8400;   heatingFunction[13] = 24.1;
        velocityFunction[14] = 8600;   heatingFunction[14] = 26.0;
        velocityFunction[15] = 8800;   heatingFunction[15] = 28.9;
        velocityFunction[16] = 9000;   heatingFunction[16] = 32.8;
        velocityFunction[17] = 9200; 	heatingFunction[17] = 35.138107;
        velocityFunction[18] = 9400; 	heatingFunction[18] = 38.156042;
        velocityFunction[19] = 9600; 	heatingFunction[19] = 41.269273;
        velocityFunction[20] = 9800; 	heatingFunction[20] = 44.477800;
        velocityFunction[21] = 10000; 	heatingFunction[21] = 47.781623;
        velocityFunction[22] = 10200; 	heatingFunction[22] = 51.180743;
        velocityFunction[23] = 10400; 	heatingFunction[23] = 54.675159;
        velocityFunction[24] = 10600; 	heatingFunction[24] = 58.264871;
        velocityFunction[25] = 10800; 	heatingFunction[25] = 61.949880;
        velocityFunction[26] = 11000; 	heatingFunction[26] = 65.730185;
        velocityFunction[27] = 11200; 	heatingFunction[27] = 69.605786;
        velocityFunction[28] = 11400; 	heatingFunction[28] = 73.576683;
        velocityFunction[29] = 11600; 	heatingFunction[29] = 77.642877;
        velocityFunction[30] = 11800; 	heatingFunction[30] = 81.804367;
        velocityFunction[31] = 12000; 	heatingFunction[31] = 86.061153;
        velocityFunction[32] = 12200; 	heatingFunction[32] = 90.413235;
        velocityFunction[33] = 12400; 	heatingFunction[33] = 94.860614;
        velocityFunction[34] = 12600; 	heatingFunction[34] = 99.403289;
        velocityFunction[35] = 12800; 	heatingFunction[35] = 104.041260;
        velocityFunction[36] = 13000; 	heatingFunction[36] = 108.774528;

        std::map< double, double > heatingMap;
        std::transform( velocityFunction.begin(), velocityFunction.end(), heatingFunction.begin(),
               std::inserter(heatingMap, heatingMap.end() ), std::make_pair<double const&,double const&> );

        boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            boost::make_shared< interpolators::InterpolatorSettings >( interpolators::cubic_spline_interpolator );

        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
                interpolators::createOneDimensionalInterpolator(
                    heatingMap, interpolatorSettings );

        double radiativeConstant_a = 0.526;
        double radiativeConstant_b = 1.19;
        double radiativeConstant_C = 2.35e8;

        double airspeed = getCurrentAirspeed();
        double density = getCurrentDensity();
        double noseRadius = vehicleSystem_->getNoseRadius();

        double radiativeHeatFlux = radiativeConstant_C * std::pow(noseRadius,radiativeConstant_a) *
                                std::pow(density,radiativeConstant_b) * interpolator->interpolate(airspeed);

        return radiativeHeatFlux;

    }

    //! Function to return atmosphere model object
    /*!
     *  Function to return atmosphere model object
     *  \return Atmosphere model object
     */
    boost::shared_ptr< aerodynamics::AtmosphereModel > getAtmosphereModel( ) const
    {
        return atmosphereModel_;
    }

    //! Function to (re)set aerodynamic angle calculator object
    /*!
     *  Function to (re)set aerodynamic angle calculator object
     *  \param aerodynamicAngleCalculator Aerodynamic angle calculator object to set.
     */
    void setAerodynamicAngleCalculator(
            const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
            aerodynamicAngleCalculator )
    {
        aerodynamicAngleCalculator_ = aerodynamicAngleCalculator;
        bodyCenteredPseudoBodyFixedStateFunction_ = boost::bind(
                    &reference_frames::AerodynamicAngleCalculator::getCurrentAirspeedBasedBodyFixedState, aerodynamicAngleCalculator_ );
    }

    //! Function to set custom dependency of aerodynamic coefficients
    /*!
     * Function to set custom dependency of aerodynamic coefficients. If needed, the
     * AerodynamicCoefficientsIndependentVariables enum may be expanded to include e.g. control surface deflections, the
     * values of which will then be retrieved from the function set here
     * \param independentVariable Identifier of independent variable
     * \param coefficientDependency Function returning the current value of the independent variable.
     */
    void setAerodynamicCoefficientsIndependentVariableFunction(
            const AerodynamicCoefficientsIndependentVariables independentVariable,
            const boost::function< double( ) > coefficientDependency );

    //! Function to return current central body-fixed state of vehicle.
    /*!
     *  Function to return central body-fixed state of vehicle.
     *  \return Current central body-fixed state of vehicle.
     */
    Eigen::Vector6d getCurrentBodyCenteredBodyFixedState( )
    {
        return currentBodyCenteredAirspeedBasedBodyFixedState_;
    }

    //! Function to return current central body-fixed velocity of vehicle.
    /*!
     *  Function to return central body-fixed velocity of vehicle.
     *  \return Current central body-fixed velocity of vehicle.
     */
    Eigen::Vector3d getCurrentAirspeedBasedVelocity( )
    {
        return currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 3, 3 );
    }


    //! Function to return aerodynamic angle calculator object
    /*!
     *  Function to return aerodynamic angle calculator object
     *  \return Aerodynamic angle calculator object
     */
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
    getAerodynamicAngleCalculator( )
    {
        return aerodynamicAngleCalculator_;
    }

    //! Function to return object from which the aerodynamic coefficients are obtained.
    /*!
     *  Function to return object from which the aerodynamic coefficients are obtained.
     *  \return Object from which the aerodynamic coefficients are obtained.
     */
    boost::shared_ptr< AerodynamicCoefficientInterface > getAerodynamicCoefficientInterface( )
    {
        return aerodynamicCoefficientInterface_;
    }

    //! Function to return list of independent variables of the aerodynamic coefficient interface
    /*!
     *  Function to return list of independent variables of the aerodynamic coefficient interface
     *  \return List of independent variables of the aerodynamic coefficient interface
     */
    std::vector< double > getAerodynamicCoefficientIndependentVariables( )
    {
        if( aerodynamicCoefficientIndependentVariables_.size( ) !=
                aerodynamicCoefficientInterface_->getNumberOfIndependentVariables( ) )
        {
            updateAerodynamicCoefficientInput( );
        }

        return aerodynamicCoefficientIndependentVariables_;
    }

    //! Function to return list of independent variables of the control surface aerodynamic coefficient interface
    /*!
     *  Function to return list of independent variables of the control surface aerodynamic coefficient interface
     *  \return List of independent variables of the control surface aerodynamic coefficient interface, with map key
     *  the control surface identifiers.
     */
    std::map< std::string, std::vector< double > > getControlSurfaceAerodynamicCoefficientIndependentVariables( )
    {
        if( controlSurfaceAerodynamicCoefficientIndependentVariables_.size( ) !=
                aerodynamicCoefficientInterface_->getNumberOfControlSurfaces( ) )
        {
            updateAerodynamicCoefficientInput( );
        }

        return controlSurfaceAerodynamicCoefficientIndependentVariables_;
    }


    //! Function to reset the current time of the flight conditions.
    /*!
     *  Function to reset the current time of the flight conditions. This function is typically sused to set the current time
     *  to NaN, indicating the need to recompute all quantities for the next time computation.
     * \param currentTime
     */
    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;

        scalarFlightConditions_.clear( );
        isLatitudeAndLongitudeSet_ = 0;

        aerodynamicAngleCalculator_->resetCurrentTime( currentTime_ );
        aerodynamicCoefficientIndependentVariables_.clear( );
        controlSurfaceAerodynamicCoefficientIndependentVariables_.clear( );
    }

private:

    //! Function to (compute and) retrieve the value of an independent variable of aerodynamic coefficients
    /*!
     * Function to (compute and) retrieve the value of an independent variable of aerodynamic coefficients
     * \param independentVariableType Identifier for type of independent variable
     * \param secondaryIdentifier String used as secondary identifier of independent variable (e.g. control surface name).
     * \return Current value of requested independent variable.
     */
    double getAerodynamicCoefficientIndependentVariable(
            const AerodynamicCoefficientsIndependentVariables independentVariableType,
            const std::string& secondaryIdentifier = "" );

    //! Function to compute and set the current latitude and longitude
    void computeLatitudeAndLongitude( )
    {
        scalarFlightConditions_[ latitude_flight_condition ] = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::latitude_angle );
        scalarFlightConditions_[ longitude_flight_condition ] = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::longitude_angle );
        isLatitudeAndLongitudeSet_ = 1;
    }

    //! Function to compute and set the current altitude
    void computeAltitude( )
    {
        scalarFlightConditions_[ altitude_flight_condition ] =
                shapeModel_->getAltitude( currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 0, 3 ) );
    }

    //! Function to update input to atmosphere model (altitude, as well as latitude and longitude if needed).
    void updateAtmosphereInput( )
    {
        if( ( scalarFlightConditions_.count( latitude_flight_condition ) == 0 ||
              scalarFlightConditions_.count( longitude_flight_condition ) == 0 ) )
        {
           if( updateLatitudeAndLongitudeForAtmosphere_ )
            {
                computeLatitudeAndLongitude( );
            }
            else
            {
                scalarFlightConditions_[ latitude_flight_condition ] = 0.0;
                scalarFlightConditions_[ longitude_flight_condition ] = 0.0;
            }
        }

        if( scalarFlightConditions_.count( altitude_flight_condition ) == 0 )
        {
            computeAltitude( );
        }
    }

    //! Function to compute and set the current freestream density
    void computeDensity( )
    {
        updateAtmosphereInput( );
        scalarFlightConditions_[ density_flight_condition ] =
                atmosphereModel_->getDensity(
                    scalarFlightConditions_.at( altitude_flight_condition ),
                    scalarFlightConditions_.at( longitude_flight_condition ),
                    scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
    }

    //! Function to compute and set the current freestream temperature
    void computeTemperature( )
    {
        updateAtmosphereInput( );

             scalarFlightConditions_[ temperature_flight_condition ] =
                     atmosphereModel_->getTemperature(
                         scalarFlightConditions_.at( altitude_flight_condition ),
                         scalarFlightConditions_.at( longitude_flight_condition ),
                         scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
    }

    //! Function to compute and set the current freestream pressure.
    void computeFreestreamPressure( )
    {
        updateAtmosphereInput( );

             scalarFlightConditions_[ pressure_flight_condition ] =
                     atmosphereModel_->getPressure(
                         scalarFlightConditions_.at( altitude_flight_condition ),
                         scalarFlightConditions_.at( longitude_flight_condition ),
                         scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
    }


    //! Function to compute and set the current speed of sound
    void computeSpeedOfSound( )
    {
        updateAtmosphereInput( );
        scalarFlightConditions_[ speed_of_sound_flight_condition ]  =
                atmosphereModel_->getSpeedOfSound(
                    scalarFlightConditions_.at( altitude_flight_condition ),
                    scalarFlightConditions_.at( longitude_flight_condition ),
                    scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
    }

    //! Function to compute and set the current airspeed
    void computeAirspeed( )
    {
        scalarFlightConditions_[ airspeed_flight_condition ] = currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 3, 3 ).norm( );
    }

    //! Function to compute and set the current freestream dynamic pressure.
    void computeDynamicPressure( )
    {
        double currentAirspeed = getCurrentAirspeed( );
        scalarFlightConditions_[ dynamic_pressure_condition ] = 0.5 *
                getCurrentDensity( ) * currentAirspeed * currentAirspeed;
    }

    //! Function to compute and set the current Mach number
    void computeMachNumber( )
    {
        scalarFlightConditions_[ mach_number_flight_condition ] =
                getCurrentAirspeed( ) / getCurrentSpeedOfSound( );
    }

    //! Function to compute and set the current geodetic latitude.
    void computeGeodeticLatitude( )
    {
        if( !geodeticLatitudeFunction_.empty( ) )
        {
            scalarFlightConditions_[ geodetic_latitude_condition ] = geodeticLatitudeFunction_(
                        currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 0, 3 ) );
        }
        else
        {
            if( scalarFlightConditions_.count( latitude_flight_condition ) == 0 || !isLatitudeAndLongitudeSet_ )
            {
                computeLatitudeAndLongitude( );
            }
            scalarFlightConditions_[ geodetic_latitude_condition ] = scalarFlightConditions_[ latitude_flight_condition ] ;
        }
    }

    //! Function to update the independent variables of the aerodynamic coefficient interface
    void updateAerodynamicCoefficientInput( );

    //! Name of central body (i.e. body with the atmosphere)
    std::string centralBody_;

    //! Atmosphere model of atmosphere through which vehicle is flying
    boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel_;

    //! Model describing the shape of the body w.r.t. which the flight is taking place.
    const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel_;

    //! Function to return the current state of the vehicle in a body-fixed frame.
    boost::function< Eigen::Vector6d( ) > bodyCenteredPseudoBodyFixedStateFunction_;

    //! Object from which the aerodynamic coefficients are obtained.
    boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface_;

    //! Object from which the aerodynamic/trajectory angles of the vehicle are calculated.
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    //! Function returning control surface deflection, with input the control surface identifier.
    boost::function< double( const std::string& ) > controlSurfaceDeflectionFunction_;

    //! Current state of vehicle in base frame for Body objects.
    Eigen::Vector6d currentBodyCenteredState_;

    //! Current state of vehicle in body-fixed frame.
    Eigen::Vector6d currentBodyCenteredAirspeedBasedBodyFixedState_;

    //! List of atmospheric/flight properties computed at current time step.
    std::map< FlightConditionVariables, double > scalarFlightConditions_;

    //! Current time of propagation.
    double currentTime_;

    //! List of custom functions for aerodynamic coefficient dependencies.
    std::map< AerodynamicCoefficientsIndependentVariables, boost::function< double( ) > > customCoefficientDependencies_;

    //! Boolean setting whether latitude and longitude are to be updated by updateConditions().
    bool updateLatitudeAndLongitudeForAtmosphere_;

    //! Boolean denoting whether the current latitude and longitude have been computed at current time step
    bool isLatitudeAndLongitudeSet_;

    //! Function from which to compute the geodetic latitude as function of body-fixed position (empty if equal to
    //! geographic latitude).
    boost::function< double( const Eigen::Vector3d& ) > geodeticLatitudeFunction_;

    //! Current list of independent variables of the aerodynamic coefficient interface
    std::vector< double > aerodynamicCoefficientIndependentVariables_;

    //! List of independent variables of the control surface aerodynamic coefficient interface, with map key
    //! control surface identifiers.
    std::map< std::string, std::vector< double > > controlSurfaceAerodynamicCoefficientIndependentVariables_;

    boost::shared_ptr< system_models::VehicleSystems > vehicleSystem_;
};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_FLIGHTCONDITIONS_H
