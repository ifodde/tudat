/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

namespace tudat
{
namespace ephemerides
{

//! Get rotation quaternion to target frame from base frame.
Eigen::Quaterniond SimpleRotationalEphemeris::getRotationToTargetFrame(
        const double secondsSinceEpoch )
{
    // Determine rotation angle compared to initial rotational state.
    double rotationAngle = basic_mathematics::computeModulo(
                ( secondsSinceEpoch - initialSecondsSinceEpoch_ ) * rotationRate_,
                2.0 * mathematical_constants::PI );

    // Calculate and return rotation to base frame.
    return reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                rotationAngle ) * initialRotationToTargetFrame_;
}

Eigen::Quaterniond  SimpleRotationalEphemeris::getRotationToTargetFrameFromExtendedTime(
        const Time timeSinceEpoch )
{
    long double rotationPeriod = 2.0L * mathematical_constants::LONG_PI /
            static_cast< long double >( rotationRate_ );
    long double timeIntoCurrentRotation = basic_mathematics::computeModuloWithLongDoubleRatio(
                ( timeSinceEpoch -  Time( initialSecondsSinceEpoch_ ) ), Time( rotationPeriod ) );
    long double rotationAngle =
                2.0L * mathematical_constants::LONG_PI * timeIntoCurrentRotation / rotationPeriod;

    // Calculate and return rotation to base frame.
    return reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                rotationAngle ) * initialRotationToTargetFrame_;
}

//! Function to calculate the derivative of the rotation matrix from target frame to original frame.
Eigen::Matrix3d SimpleRotationalEphemeris::getDerivativeOfRotationToTargetFrame(
        const double secondsSinceEpoch )
{
    // Determine rotation angle compared to initial rotational state.
    double rotationAngle = basic_mathematics::computeModulo(
                ( secondsSinceEpoch - initialSecondsSinceEpoch_ ) * rotationRate_,
                2.0 * mathematical_constants::PI );

    // Calculate derivative of rotation matrix.
    return rotationRate_ * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER * tudat::reference_frames::
            getInertialToPlanetocentricFrameTransformationQuaternion( rotationAngle )
            * Eigen::Matrix3d( initialRotationToTargetFrame_ );
}

Eigen::Matrix3d SimpleRotationalEphemeris::getDerivativeOfRotationToTargetFrameFromExtendedTime(
        const Time timeSinceEpoch )
{
    long double rotationPeriod = 2.0L * mathematical_constants::LONG_PI /
            static_cast< long double >( rotationRate_ );
    long double timeIntoCurrentRotation = basic_mathematics::computeModuloWithLongDoubleRatio(
                ( timeSinceEpoch -  Time( initialSecondsSinceEpoch_ ) ),
                Time( rotationPeriod ) );
    double rotationAngle =  2.0L * mathematical_constants::LONG_PI * timeIntoCurrentRotation / rotationPeriod;

    // Calculate derivative of rotation matrix.
    return rotationRate_ * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER * tudat::reference_frames::
            getInertialToPlanetocentricFrameTransformationQuaternion( rotationAngle )
            * Eigen::Matrix3d( initialRotationToTargetFrame_ );
}

//! Function to reset the right ascension and declination of body's north pole.
void SimpleRotationalEphemeris::resetInitialPoleRightAscensionAndDeclination(
        const double rightAscension, const double declination )
{
    // Recalculate initial rotation quaternion
    initialRotationToTargetFrame_ =
        reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
            declination, rightAscension, initialEulerAngles_.z( ) );

    // Reset angles in vector of Euler angles.
    initialEulerAngles_.x( ) = rightAscension;
    initialEulerAngles_.y( ) = declination;
}

} // namespace tudat
} // namespace ephemerides
