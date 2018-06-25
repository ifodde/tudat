/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      References for deltaV computation function:
 *          Melman J. Trajectory optimization for a mission to Neptune and Triton, MSc thesis
 *              report, Delft University of Technology, 2007.
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far].
 *      Reference for unpowered gravity assist propagation function:
 *          Conway, B.A., Spacecraft Trajectory Optimization, Chapter 7, Cambridge University
 *              Press, 2010.
 *      Reference for powered gravity assist propagation function:
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far].
 *
 *    Notes
 *      Gravity assist and swing-by are two different words for the same thing. The delta-V that is
 *      computed for a powered swing-by has not been proven to be the optimum (lowest) to achieve
 *      the desired geometry of incoming and outgoing hyperbolic legs. Some literature research will
 *      have to be done to look at the alternatives.
 *
 *      Note that the exact implementation of Newton Raphson as root finder should be updated if
 *      someone would want to use a different root finding technique.
 *
 *      Note that by default a velocity effect deltaV of less than 1 micrometer/second is deemed
 *      negligable in this code. This value can be set though.
 *
 */

#include <cmath>

#include <boost/bind.hpp>

#include <Eigen/Dense>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/MissionSegments/gravityAssist.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"

namespace tudat
{
namespace mission_segments
{

//! Calculates Mars' limiting function
double MarsAtmosphericTrajectoryFunction( const double incomingHyperbolicVelocity,
                                          const double velocityBendingAngle )
{
    if (incomingHyperbolicVelocity > 7000 || incomingHyperbolicVelocity < 2000)
    {
        return -10000.0;
    }
    double p00 = 2380.0;
    double p10 = 0.85;
    double p01 = -1048;
    double p11 = -0.3385;
    double p02 = 641.6;
    double p12 = 0.03174;
    double p03 = -160.4;

    double x = incomingHyperbolicVelocity;
    double y = velocityBendingAngle;

    double outgoingHyperbolicVelocity = p00 + p10*x + p01*y + p11*x*y + p02*std::pow(y,2) + p12*x*std::pow(y,2) + p03*pow(y,3);

    if (outgoingHyperbolicVelocity >= 0 && outgoingHyperbolicVelocity < incomingHyperbolicVelocity){
        return outgoingHyperbolicVelocity;
    } else {
        return -10000.0;
    }


}

//! Calculates Earth's limiting function
double EarthAtmosphericTrajectoryFunction( const double incomingHyperbolicVelocity,
                                          const double velocityBendingAngle )
{
    if ( (incomingHyperbolicVelocity > 6000 || incomingHyperbolicVelocity < 2000) ||
         (incomingHyperbolicVelocity > 5000 && velocityBendingAngle > 100*mathematical_constants::PI / 180.0) )
    {
        return -10000.0;
    }
    double p00 = -2185.0;
    double p10 = 1.801;
    double p01 = 9602;
    double p11 = -1.495;
    double p02 = -4120.0;
    double p12 = 0.2912;
    double p03 = 368.0;

    double x = incomingHyperbolicVelocity;
    double y = velocityBendingAngle;

    double outgoingHyperbolicVelocity = p00 + p10*x + p01*y + p11*x*y + p02*std::pow(y,2) + p12*x*std::pow(y,2) + p03*pow(y,3);

    if (outgoingHyperbolicVelocity >= 0 && outgoingHyperbolicVelocity < incomingHyperbolicVelocity){
        return outgoingHyperbolicVelocity;
    } else {
        return -10000.0;
    }


}

//! Calculates Venus' limiting function
double VenusAtmosphericTrajectoryFunction( const double incomingHyperbolicVelocity,
                                          const double velocityBendingAngle )
{
    if (incomingHyperbolicVelocity > 6000 || incomingHyperbolicVelocity < 4000)
    {
        return -10000.0;
    }
    double p00 = 27010.0;
    double p10 = -1.36;
    double p01 = -13660.0;
    double p11 = 0.5231;
    double p02 = 1470;

    double x = incomingHyperbolicVelocity;
    double y = velocityBendingAngle;

    double outgoingHyperbolicVelocity = p00 + p10*x + p01*y + p11*x*y + p02*std::pow(y,2);

    if (outgoingHyperbolicVelocity >= 0 && outgoingHyperbolicVelocity < incomingHyperbolicVelocity && outgoingHyperbolicVelocity <= 3.0E3){
        return outgoingHyperbolicVelocity;
    } else {
        return -10000.0;
    }


}

//! Calculate deltaV of an aero gravity assist.
double AeroGravityAssist( const std::string centralBodyName,
                          const Eigen::Vector3d& centralBodyVelocity,
                          const Eigen::Vector3d& incomingVelocity,
                          const Eigen::Vector3d& outgoingVelocity)
{

    // Compute incoming and outgoing hyperbolic excess velocity.
    const Eigen::Vector3d incomingHyperbolicExcessVelocity
            = incomingVelocity - centralBodyVelocity;
    const Eigen::Vector3d outgoingHyperbolicExcessVelocity
            = outgoingVelocity - centralBodyVelocity;

    // Compute absolute values of the hyperbolic excess velocities.
    const double absoluteIncomingExcessVelocity = incomingHyperbolicExcessVelocity.norm( );
    const double absoluteOutgoingExcessVelocity = outgoingHyperbolicExcessVelocity.norm( );

    // Compute bending angle.
    double bendingAngle = linear_algebra::computeAngleBetweenVectors(
                            incomingHyperbolicExcessVelocity, outgoingHyperbolicExcessVelocity );

    // Calculate limit of outgoing hyperbolic exces velocity
    double maxAGAOutgoingExcessVelocity = 0.0;
    if ( centralBodyName == "Mars"  )
    {
        maxAGAOutgoingExcessVelocity = MarsAtmosphericTrajectoryFunction( absoluteIncomingExcessVelocity, bendingAngle );
    } else if ( centralBodyName == "Earth"  )
    {
        maxAGAOutgoingExcessVelocity = EarthAtmosphericTrajectoryFunction( absoluteIncomingExcessVelocity, bendingAngle );
    } else if ( centralBodyName == "Venus"  )
    {
        maxAGAOutgoingExcessVelocity = VenusAtmosphericTrajectoryFunction( absoluteIncomingExcessVelocity, bendingAngle );
    } else
    {
        std::cerr<<"Planet does not have a limiting function.";
    }

    double deltaV = 0;

    if (maxAGAOutgoingExcessVelocity < absoluteOutgoingExcessVelocity)
    {
        deltaV += absoluteOutgoingExcessVelocity - maxAGAOutgoingExcessVelocity;
    }

    return deltaV;

}

} // namespace mission_segments
} // namespace tudat
