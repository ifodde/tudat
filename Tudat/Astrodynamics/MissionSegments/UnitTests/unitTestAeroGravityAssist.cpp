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
 *      Izzo, D. and Vinko, T. ACT - Informatics - GTOP Database, ESA Advanced Concept Team, last
 *          accessed on 2012-01-12. http://www.esa.int/gsp/ACT/inf/op/globopt.htm.
 *      Musegaas, P. Gravity Assist calculation Verification.xlsx, last accessed: 3 December 2012,
 *          http://tudat.tudelft.nl/projects/tudat/wiki/Unit_tests, 2012.
 *
 *    Notes
 *      Three main functions are tested in these unit tests.
 *        Regarding the deltaV calculation gravity assist method:
 *          There is a complicated if-statement in this method. Hence many unit test are performed
 *          to test the functionality. Also various limit cases failed previously, hence many tests
 *          for this are also included:
 *              Case 1: required bending angle > maximum bending angle:
 *                  Two tests were written. In the first one no velocity effect is needed. This
 *                  test has a low accuracy, which should be replaced one day (it still relies on
 *                  hand calculator calculations done in 2011). In the second one a combination of
 *                  bending-effect deltaV and velocity-effect deltaV is calculated. This test has
 *                  been calculated using Tudat, and was verified using Excel.
 *                  Could definitely be improved.
 *              Case 2: no assist is required:
 *                  One test was written.
 *              Case 3: velocity effect deltaV only, using eccentricity iteration scheme:
 *                  Four tests were written. The first one calculates a case from Cassini-1 of GTOP
 *                  with high precision. The other three test limit cases: low incoming, high
 *                  outgoing velocity; high incoming, low outgoing velocity; low incoming, low
 *                  outgoing velocity. These tests were calculated using Tudat, but verified in
 *                  Excel to be exactly correct.
 *              Case 4: velocity effect deltaV only, using pericenter radius iteration scheme:
 *                  The same four tests as for case 3 were used.
 *        Regarding the unpowered gravity assist propagator:
 *          One test was written, based on GTOP. This should be a satisfactory test.
 *        Regarding the powered gravity assist propagator:
 *          Two tests were written. The first one is similar to the unpowered gravity assist
 *          propagator. The second one is reverse engineered from the Cassini-1 test, similar to
 *          the one in the deltaV calculation test.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <iostream>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/MissionSegments/aeroGravityAssist.h"

namespace tudat
{
namespace unit_tests
{

//! Test of aerogravity assist code.
BOOST_AUTO_TEST_SUITE( test_aerogravity_assist )

//! Test bending angle Delta-V computation.
BOOST_AUTO_TEST_CASE( testExtraDeltaV )
{
    // Tolerance, determined primarily by the accuracy of the hand calculations for this test case.
    const double velocityTolerance = 0.02;

    // Expected delta-V for a powered aga around Mars.
    const double expectedDeltaV = 170.9209;

    // Define Sun gravitational parameter.
    const double gravitationalParameterSun = 1.32712440018e20;

    // Define planet-Sun distance.
    const double distanceMarsToSun = unit_conversions::
            convertAstronomicalUnitsToMeters( 1.5 );

    // Define planet heliocentric velocity vector. The orbit is considered to be circular.
    const Eigen::Vector3d marsVelocity( 0.0,
                                        std::sqrt( gravitationalParameterSun / distanceMarsToSun ),
                                        0.0 );

    // Define satellite incoming vector.
    using mathematical_constants::PI;
    const Eigen::Vector3d incomingVelocity( 4470.0,
                                            4000.0 + marsVelocity( 1 ),
                                            0.0 );

    // Define satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( -2969.50444,
                                            426.645399 + marsVelocity( 1 ),
                                            0.0 );

    // Perform the gravity assist.
    const double deltaV = mission_segments::AeroGravityAssist( "Mars",
                                                           marsVelocity, incomingVelocity,
                                                           outgoingVelocity);


    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, velocityTolerance );
}

//! Test bending angle Delta-V computation.
BOOST_AUTO_TEST_CASE( testNoDeltaV )
{
    // Tolerance, determined primarily by the accuracy of the hand calculations for this test case.
    const double velocityTolerance = 0.02;

    // Expected delta-V for a aga around Mars.
    const double expectedDeltaV = 0.0;

    // Define Sun gravitational parameter.
    const double gravitationalParameterSun = 1.32712440018e20;

    // Define planet-Sun distance.
    const double distanceMarsToSun = unit_conversions::
            convertAstronomicalUnitsToMeters( 1.5 );

    // Define planet heliocentric velocity vector. The orbit is considered to be circular.
    const Eigen::Vector3d marsVelocity( 0.0,
                                        std::sqrt( gravitationalParameterSun / distanceMarsToSun ),
                                        0.0 );

    // Define satellite incoming vector.
    using mathematical_constants::PI;
    const Eigen::Vector3d incomingVelocity( 4470.0,
                                            4000.0 + marsVelocity( 1 ),
                                            0.0 );

    // Define satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( -1581.93,
                                            2549.0147 + marsVelocity( 1 ),
                                            0.0 );

    // Perform the gravity assist.
    const double deltaV = mission_segments::AeroGravityAssist( "Mars",
                                                           marsVelocity, incomingVelocity,
                                                           outgoingVelocity);

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, velocityTolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
