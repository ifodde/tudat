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

#ifndef TUDAT_ATMOSPHERE_MODEL_H
#define TUDAT_ATMOSPHERE_MODEL_H

#include <boost/shared_ptr.hpp>

#include <iostream>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/windModel.h"

namespace tudat
{
namespace aerodynamics
{

//! Atmosphere model class.
/*!
 * Base class for all atmosphere models.
 * To keep the function generic for both reference atmospheres and standard
 * atmospheres, altitude, longitude, latitude, and time are inputs for all functions.
 */
class AtmosphereModel
{
public:

    //! Default destructor.
    /*!
    * Default destructor.
    */
    virtual ~AtmosphereModel( ) { }

    //! Get local density.
    /*!
    * Returns the local density parameter of the atmosphere in kg per meter^3.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric density.
    */
    virtual double getDensity( const double altitude, const double longitude,
                               const double latitude, const double time ) = 0;

    //! Get local pressure.
    /*!
    * Returns the local pressure of the atmosphere parameter in Newton per meter^2.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric pressure.
    */
    virtual double getPressure( const double altitude, const double longitude,
                                const double latitude, const double time ) = 0;

    //! Get local temperature.
    /*!
    * Returns the local temperature of the atmosphere parameter in Kelvin.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric temperature.
    */
    virtual double getTemperature( const double altitude, const double longitude,
                                   const double latitude, const double time ) = 0;

    //! Get local speed of sound.
    /*!
    * Returns the local speed of sound of the atmosphere in m/s.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric speed of sound.
    */
    virtual double getSpeedOfSound( const double altitude, const double longitude,
                                    const double latitude, const double time ) = 0;

    //! Function to retrieve the model describing the wind velocity vector of the atmosphere
    /*!
     * Function to retrieve the model describing the wind velocity vector of the atmosphere
     * \return Model describing the wind velocity vector of the atmosphere
     */
    boost::shared_ptr< WindModel > getWindModel( )
    {
        return windModel_;
    }

    //! Function to set the model describing the wind velocity vector of the atmosphere
    /*!
     * Function to set the model describing the wind velocity vector of the atmosphere
     * \param windModel New model describing the wind velocity vector of the atmosphere
     */
    void setWindModel( const boost::shared_ptr< WindModel > windModel )
    {
        windModel_ = windModel;
    }
    //! Get local mean molar mass.
    /*!
    * Returns the local mean molar mass in kg/mol.
    * \param altitude Altitude at which mean molar mass is to be computed [m].
    * \param longitude Longitude at which mean molar mass  is to be computed [rad].
    * \param latitude Latitude at which mean molar mass  is to be computed [rad].
    * \param time Time at which mean molar mass  is to be computed (seconds since J2000).
    * \return mean molar mass.
    */
    virtual double getMeanMolarMass( const double altitude, const double longitude,
                          const double latitude, const double time )
    {
        std::cout<<"here in flight conditions 1"<<std::endl;
        throw std::runtime_error("ERROR: no getMeanMolarMass for this atmosphere model" );
        return 0;
    }
    //! Get local weighted average collision diameter.
    /*!
    * Returns the local weighted average collision diameter using the number densities as weights.
    * \param altitude Altitude at which average collision diameter is to be computed [m].
    * \param longitude Longitude at which average collision diameter is to be computed [rad].
    * \param latitude Latitude at which average collision diameter is to be computed [rad].
    * \param time Time at which average collision diameter is to be computed (seconds since J2000).
    * \return weighted average collision diameter.
    */
    virtual double getWeightedAverageCollisionDiameter( const double altitude, const double longitude,
                          const double latitude, const double time )
    {
        std::cout<<"here in flight conditions 2"<<std::endl;
        throw std::runtime_error("ERROR: no getCollisionDiameter for this atmosphere model" );
        return 0;
    }

protected:

    //! Model describing the wind velocity vector of the atmosphere
    boost::shared_ptr< WindModel > windModel_;
private:
};

//! Typedef for shared-pointer to AtmosphereModel object.
typedef boost::shared_ptr< AtmosphereModel > AtmosphereModelPointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_ATMOSPHERE_MODEL_H
