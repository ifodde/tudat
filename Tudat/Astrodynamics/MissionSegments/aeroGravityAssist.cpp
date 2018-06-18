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
    if (incomingHyperbolicVelocity > 9000 || incomingHyperbolicVelocity < 3000)
    {
        return 0;
    }
    double p00 = 1022.0;//  (524.1, 1520)
    double p10 = 1.008;//  (0.9565, 1.059)
    double p01 = 913.9;//  (324.6, 1503)
    double p11 = -0.5519;//  (-0.597, -0.5068)
    double p02 = -7.108;//  (-234.6, 220.4)
    double p12 = 0.07504;//  (0.06548, 0.0846)
    double p03 = -92.47;//  (-123.6, -61.38)

    double x = incomingHyperbolicVelocity;
    double y = velocityBendingAngle;

    if (p00 + p10*x + p01*y + p11*x*y + p02*std::pow(y,2) + p12*x*std::pow(y,2) + p03*pow(y,3) >= 0){
        return p00 + p10*x + p01*y + p11*x*y + p02*std::pow(y,2) + p12*x*std::pow(y,2) + p03*pow(y,3);
    } else {
        return 0.0;
    }


}

//! Calculate deltaV of an aero gravity assist.
double AeroGravityAssist( const std::string centralBodyName,
                          const Eigen::Vector3d& centralBodyVelocity,
                          const Eigen::Vector3d& incomingVelocity,
                          const Eigen::Vector3d& outgoingVelocity)
{
    // Check if Body is available for aerogravity assist.
    if (centralBodyName != "Mars")
    {
        std::cerr<<"Planet does not have a limiting function.";
    }

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
    double maxAGAOutgoingExcessVelocity = MarsAtmosphericTrajectoryFunction( absoluteIncomingExcessVelocity, bendingAngle );

    double deltaV = 0;

    // ADD CHECK FOR NEGATIVE maxAGAOutgoingExcessVelocity

    if (maxAGAOutgoingExcessVelocity < absoluteOutgoingExcessVelocity)
    {
        deltaV += absoluteOutgoingExcessVelocity - maxAGAOutgoingExcessVelocity;
    }

    return deltaV;

}

} // namespace mission_segments
} // namespace tudat
