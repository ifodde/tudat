/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 *      By default the eccentricity is used as the iteration procedure. This is because in
 *      optimizing a Cassini-like trajectory, the pericenter radius had about 2-4 NaN values in
 *      100000 times the gravity assist calculation. The eccentricity iteration had no NaN values
 *      for a similar run in which 100000 gravity assist calculations were done. Also the
 *      eccentricity seemed to require slightly less iterations (does not necessarily mean it is
 *      faster or more accurate).
 *
 */

#ifndef TUDAT_GRAVITY_ASSIST_H
#define TUDAT_GRAVITY_ASSIST_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

namespace tudat
{
namespace mission_segments
{

//!
/*!
 * WARNING: these functions only work for a specific waverider vehicle, if a different waverider is used
 * then the one in (Hess 2017), the functions are not accurate.
 */

//! Calculate limiting surface of an aero gravity assist at Mars.
/*!
 * Calculates if the input values are feasible for an aerogravity assist, and if not, calculates the
 * extra Delta V needed to correct for it.
 * \param incomingHyperbolicVelocity Hyperbolic excess velocity of the spacecraft before the swing-by.     [m s^-1]
 * \param velocityBendingAngle Angle between the incoming and outgoing hyperbolic excess velocity vectors.    [rad]
 * \return outgoingHyperbolicVelocity the maximum outgoing hyperbolic velocity possible.                   [m s^-1]
 */
double MarsAtmosphericTrajectoryFunction( const double incomingHyperbolicVelocity,
                                          const double velocityBendingAngle );

//! Calculate limiting surface of an aero gravity assist at Earth.
/*!
 * Calculates if the input values are feasible for an aerogravity assist, and if not, calculates the
 * extra Delta V needed to correct for it.
 * \param incomingHyperbolicVelocity Hyperbolic excess velocity of the spacecraft before the swing-by.     [m s^-1]
 * \param velocityBendingAngle Angle between the incoming and outgoing hyperbolic excess velocity vectors.    [rad]
 * \return outgoingHyperbolicVelocity the maximum outgoing hyperbolic velocity possible.                   [m s^-1]
 */
double EarthAtmosphericTrajectoryFunction( const double incomingHyperbolicVelocity,
                                          const double velocityBendingAngle );

//! Calculate limiting surface of an aero gravity assist at Venus.
/*!
 * Calculates if the input values are feasible for an aerogravity assist, and if not, calculates the
 * extra Delta V needed to correct for it.
 * \param incomingHyperbolicVelocity Hyperbolic excess velocity of the spacecraft before the swing-by.     [m s^-1]
 * \param velocityBendingAngle Angle between the incoming and outgoing hyperbolic excess velocity vectors.    [rad]
 * \return outgoingHyperbolicVelocity the maximum outgoing hyperbolic velocity possible.                   [m s^-1]
 */
double VenusAtmosphericTrajectoryFunction( const double incomingHyperbolicVelocity,
                                          const double velocityBendingAngle );

//! Calculate deltaV of an aero gravity assist.
/*!
 * Calculates if the input values are feasible for an aerogravity assist, and if not, calculates the
 * extra Delta V needed to correct for it.
 * \param centralBodyName Name of the swingby body.                                             [-]
 * \param centralBodyVelocity Heliocentric velocity of the swing-by body.                  [m s^-1]
 * \param incomingVelocity Heliocentric velocity of the spacecraft before the swing-by.    [m s^-1]
 * \param outgoingVelocity Heliocentric velocity of the spacecraft after the swing-by.     [m s^-1]
 * \return deltaV The deltaV required for the gravity assist maneuver.                     [m s^-1]
 */
double AeroGravityAssist( const std::string centralBodyName,
                          const Eigen::Vector3d& centralBodyVelocity,
                          const Eigen::Vector3d& incomingVelocity,
                          const Eigen::Vector3d& outgoingVelocity );

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_GRAVITY_ASSIST_H
