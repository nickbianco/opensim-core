/* -------------------------------------------------------------------------- *
 *                  OpenSim:  testScholz2015GeometryPath.cpp                  *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2025 Stanford University and the Authors                *
 * Author(s): Nicholas Bianco                                                 *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include <OpenSim/Actuators/ModelFactory.h>
#include <OpenSim/Simulation/Manager/Manager.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PathSpring.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>
#include <OpenSim/Simulation/SimbodyEngine/SliderJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/FreeJoint.h>

#include <catch2/catch_all.hpp>

using namespace OpenSim;

TEST_CASE("Interface") {
    Model model = ModelFactory::createDoublePendulum();

    Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
    path->setName("path");
    path->setOrigin(model.getGround(), SimTK::Vec3(0.05, 0.05, 0.));
    path->setInsertion(model.getComponent<Body>("/bodyset/b1"),
            SimTK::Vec3(-0.25, 0.1, 0.));
    model.addComponent(path);

    SECTION("No obstacles or via points") {
        CHECK(path->getNumViaPoints() == 0);
        CHECK(path->getNumObstacles() == 0);
        CHECK(path->getProperty_segments().size() == 1);
    }

    // Add a ContactCylinder wrapping obstacle to the path.
    auto* obstacle = new ContactCylinder(0.1,
        SimTK::Vec3(-0.5, 0.1, 0.), SimTK::Vec3(0),
        model.getComponent<Body>("/bodyset/b0"));
    model.addComponent(obstacle);
    path->addObstacle(*obstacle, SimTK::Vec3(0., 0.1, 0.));

    SECTION("One obstacle, no via points") {
        CHECK(path->getNumViaPoints() == 0);
        CHECK(path->getNumObstacles() == 1);
        CHECK(path->getProperty_segments().size() == 1);
    }

    // Add a via point to the path.
    path->addViaPoint(model.getComponent<Body>("/bodyset/b1"),
            SimTK::Vec3(-0.75, 0.1, 0.));

    SECTION("One obstacle, one via point") {
        CHECK(path->getNumViaPoints() == 1);
        CHECK(path->getNumObstacles() == 1);
        CHECK(path->getProperty_segments().size() == 2);
    }

    SECTION("Invalid algorithm") {
        path->setAlgorithm("invalid");
        CHECK_THROWS_WITH(model.initSystem(),
            Catch::Matchers::ContainsSubstring(
                    "Property 'algorithm' has invalid value invalid; expected "
                    "one of the following: MinimumLength, Scholz2015"));
    }

    SECTION("Invalid curve segment accuracy") {
        path->setCurveSegmentAccuracy(-1.0);
        CHECK_THROWS_WITH(model.initSystem(),
            Catch::Matchers::ContainsSubstring(
                    "Scholz2015GeometryPath::extendFinalizeFromProperties: "
                    "path: curve_segment_accuracy must be greater than zero"));
    }

    SECTION("Invalid smoothness tolerance") {
        path->setSmoothnessTolerance(-1.0);
        CHECK_THROWS_WITH(model.initSystem(),
            Catch::Matchers::ContainsSubstring(
                    "Scholz2015GeometryPath::extendFinalizeFromProperties: "
                    "path: smoothness_tolerance must be greater than zero"));
    }

    SECTION("Invalid solver max iterations") {
        path->setSolverMaxIterations(-1);
        CHECK_THROWS_WITH(model.initSystem(),
            Catch::Matchers::ContainsSubstring(
                    "Scholz2015GeometryPath::extendFinalizeFromProperties: "
                    "path: solver_max_iterations must be greater than zero"));
    }
}

TEST_CASE("Suspended pendulum") {

    // Create a single pendulum model.
    Model model = ModelFactory::createPendulum();

    // Add a path spring that will wrap over each ContactGeometry, suspending
    // the pendulum. We use a relatively large dissipation coefficient to ensure
    // that the pendulum comes to rest quickly.
    auto* actu = new PathSpring();
    actu->setName("path_spring");
    actu->setRestingLength(0.0);
    actu->setDissipation(5.0);
    actu->setStiffness(25.0);
    actu->set_path(Scholz2015GeometryPath());
    model.addComponent(actu);

    // Set the path's origin and insertion.
    Scholz2015GeometryPath& path = actu->updPath<Scholz2015GeometryPath>();
    path.setOrigin(model.getGround(), SimTK::Vec3(-0.1, 0, 0));
    path.setInsertion(model.getComponent<Body>("/bodyset/b0"),
            SimTK::Vec3(-0.5, 0.1, 0));

    // Check that pendulum is at rest with the expected length.
    const double expectedLength = 0.81268778;
    auto simulateHangingPendulum = [&](Model& model) {
        SimTK::State state = model.initSystem();
        Manager manager(model);
        manager.initialize(state);
        state = manager.integrate(20.0);

        model.realizeVelocity(state);
        CHECK_THAT(path.getLength(state),
            Catch::Matchers::WithinRel(expectedLength, 1e-3));
        CHECK_THAT(path.getLengtheningSpeed(state),
            Catch::Matchers::WithinAbs(0.0, 1e-3));
        CHECK_THAT(state.getU()[0],
            Catch::Matchers::WithinAbs(0.0, 1e-3));
    };

    // Simulate the "suspended" pendulum over each contact geometry. The
    // geometries are configured such that each test will result in the same
    // final path length.

    SECTION("ContactCylinder") {
        auto* obstacle = new ContactCylinder(0.1,
            SimTK::Vec3(0.25, 0, 0), SimTK::Vec3(0), model.getGround());
        model.addComponent(obstacle);
        path.addObstacle(*obstacle, SimTK::Vec3(0.0, 0.1, 0.0));
        simulateHangingPendulum(model);
    }

    SECTION("ContactEllipsoid") {
        auto* obstacle = new ContactEllipsoid(SimTK::Vec3(0.1, 0.1, 0.3),
            SimTK::Vec3(0.25, 0, 0), SimTK::Vec3(0), model.getGround());
        model.addComponent(obstacle);
        path.addObstacle(*obstacle, SimTK::Vec3(0.0, 0.1, 0.0));
        simulateHangingPendulum(model);
    }

    SECTION("ContactSphere") {
        auto* obstacle = new ContactSphere(0.1,
            SimTK::Vec3(0.25, 0, 0), model.getGround());
        model.addComponent(obstacle);
        path.addObstacle(*obstacle, SimTK::Vec3(0.0, 0.1, 0.0));
        simulateHangingPendulum(model);
    }

    SECTION("ContactTorus") {
        auto* obstacle = new ContactTorus(0.4, 0.1,
            SimTK::Vec3(0.25, 0.4, 0),
            SimTK::Vec3(0, 0.5*SimTK::Pi, 0),
            model.getGround());
        model.addComponent(obstacle);
        path.addObstacle(*obstacle, SimTK::Vec3(0.0, -0.3, 0.0));
        simulateHangingPendulum(model);
    }
}

TEST_CASE("Moment arms") {

    // Randomly choose the radius of the cylinder within a valid range.
    SimTK::Random::Uniform random(0.1, 0.2);
    const SimTK::Real radius = random.getValue();
    const SimTK::Real offset = 0.5;

    // Construct a double pendulum model.
    Model model = ModelFactory::createDoublePendulum();

    // Configure the Scholz2015GeometryPath.
    Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
    path->setName("path");
    path->setOrigin(model.getComponent<Body>("/bodyset/b0"),
            SimTK::Vec3(-offset, 0., 0.));
    path->setInsertion(model.getComponent<Body>("/bodyset/b1"),
            SimTK::Vec3(-offset, 0., 0.));

    // Add a ContactCylinder along the axis definiing the PinJoint between
    // body 'b0' and body 'b1'.
    auto* obstacle = new ContactCylinder(radius,
        SimTK::Vec3(0., 0., 0.), SimTK::Vec3(0),
        model.getComponent<Body>("/bodyset/b0"));
    model.addComponent(obstacle);
    path->addObstacle(*obstacle, SimTK::Vec3(0., radius, 0.));

    // Add the path to the model.
    model.addComponent(path);

    // The model is in the default configuration with both bodies aligned along
    // their x-axes. Check that the computed moment arm is equal to the cylinder
    // radius.
    SECTION("Zero degrees") {
        SimTK::State state = model.initSystem();
        model.realizePosition(state);
        const Coordinate& coord = model.getCoordinateSet().get("q1");
        CHECK_THAT(path->computeMomentArm(state, coord),
                Catch::Matchers::WithinRel(radius, 1e-6));
    }

    // The model is now posed with the second body rotated 90 degrees relative
    // to the first body. The path should be lifted off the cylinder, no matter
    // what cylinder radius is used in the range defined above.
    SECTION("90 degrees") {
        model.updCoordinateSet().get("q1").setDefaultValue(SimTK::Pi/2.);
        SimTK::State state = model.initSystem();
        model.realizePosition(state);
        const Coordinate& coord = model.getCoordinateSet().get("q1");

        SimTK::Real pathLength = std::sqrt(2*offset*offset);
        CHECK_THAT(path->getLength(state),
                Catch::Matchers::WithinRel(pathLength, 1e-6));

        SimTK::Real momentArm =
            std::sqrt(offset*offset - 0.25*pathLength*pathLength);
        CHECK_THAT(path->computeMomentArm(state, coord),
                Catch::Matchers::WithinRel(momentArm, 1e-6));
    }

    // The model is now posed with the second body rotated 180 degrees relative
    // to the first body, such that the two bodies are fully overlapping. In
    // this configuration, the path length should be zero and the moment arm
    // should be equal to the offset used for the origin and insertion points.
    SECTION("180 degrees") {
        model.updCoordinateSet().get("q1").setDefaultValue(SimTK::Pi);
        SimTK::State state = model.initSystem();
        model.realizePosition(state);
        const Coordinate& coord = model.getCoordinateSet().get("q1");

        CHECK_THAT(path->getLength(state),
                Catch::Matchers::WithinAbs(0.0, 1e-6));
        CHECK_THAT(path->computeMomentArm(state, coord),
                Catch::Matchers::WithinRel(offset, 1e-6));
    }
}

// Test a simple Scholz2015GeometryPath with a known solution. The path wraps
// over the following obstacles, in order: torus, ellipsoid, torus, cylinder.
// We wrap the cable conveniently over the obstacles such that each curve
// segment becomes a circular-arc shape. This allows us to check the results by
// hand.
//
// This test is based on "testSimpleCable" in Simbody's TestCableSpan.cpp.
TEST_CASE("Simple path") {

    // Create an empty model and path object.
    Model model;
    Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
    path->setOrigin(model.getGround(), SimTK::Vec3(0.));
    path->setInsertion(model.getGround(), SimTK::Vec3(-4., 0., 0.));
    path->setCurveSegmentAccuracy(1e-12);
    path->setSmoothnessTolerance(1e-7);

    // Add the first torus obstacle.
    auto* torus1 = new ContactTorus(10., 1.,
        SimTK::Vec3(1., 11., 0.), SimTK::Vec3(0., 0.5*SimTK::Pi, 0.),
        model.getGround());
    torus1->setName("torus1");
    model.addComponent(torus1);
    path->addObstacle(*torus1, SimTK::Vec3(0., -9., 0.));

    // Add the ellipsoid obstacle.
    auto* ellipsoid = new ContactEllipsoid(SimTK::Vec3(1., 1., 6.),
        SimTK::Vec3(4.5, 1., 0.), SimTK::Vec3(0), model.getGround());
    ellipsoid->setName("ellipsoid");
    model.addComponent(ellipsoid);
    path->addObstacle(*ellipsoid, SimTK::Vec3(1., 1., 0.5));

    // Add the second torus obstacle.
    auto* torus2 = new ContactTorus(10., 1.5,
        SimTK::Vec3(4., -12., 0.), SimTK::Vec3(0., 0.5*SimTK::Pi, 0.),
        model.getGround());
    torus2->setName("torus2");
    model.addComponent(torus2);
    path->addObstacle(*torus2, SimTK::Vec3(0., 1., 0.));

    // Add the cylinder obstacle.
    auto* cylinder = new ContactCylinder(2.,
        SimTK::Vec3(-2., -1.5, 0.), SimTK::Vec3(0), model.getGround());
    cylinder->setName("cylinder");
    model.addComponent(cylinder);
    path->addObstacle(*cylinder, SimTK::Vec3(0., -1., 0.));

    // Add the path to the model.
    model.addComponent(path);

    // Initialize the system and state.
    SimTK::State state = model.initSystem();
    model.realizePosition(state);

    SECTION("Configuration and smoothness") {
        CHECK(path->getNumObstacles() == 4);
        CHECK(path->getNumViaPoints() == 0);
        CHECK(path->getSmoothness(state) <= path->getSmoothnessTolerance());
    }

    SECTION("Path length") {
        // Note that the length deviates because the path is solved up to angle
        // tolerance, not length tolerance.
        const SimTK::Real lengthTolerance = 1e-5;

        auto assertCurveSegmentLength =
            [&](SimTK::CableSpanObstacleIndex obsIx, SimTK::Real obsRadius)
        {
            REQUIRE(path->isInContactWithObstacle(state, obsIx));

            const SimTK::Real angle = 0.5 * SimTK::Pi;
            const SimTK::Real expectedLength = angle * obsRadius;
            const SimTK::Real gotLength =
                    path->calcCurveSegmentArcLength(state, obsIx);

            CHECK_THAT(gotLength, Catch::Matchers::WithinRel(
                    expectedLength, lengthTolerance));
        };
        assertCurveSegmentLength(SimTK::CableSpanObstacleIndex(0), 1.);
        assertCurveSegmentLength(SimTK::CableSpanObstacleIndex(1), 1.);
        assertCurveSegmentLength(SimTK::CableSpanObstacleIndex(2), 1.5);
        assertCurveSegmentLength(SimTK::CableSpanObstacleIndex(3), 2.);

        // Sum all straight line segment lengths + curve line segment lengths.
        const SimTK::Real sumStraightLineSegmentLengths =
                1. + 3.5 + 3. + 6. + 1.5;
        const SimTK::Real sumCurveLineSegmentLengths =
                0.5 * SimTK::Pi * (1. + 1. + 1.5 + 2.);
        const SimTK::Real expectedTotalPathLength =
                sumStraightLineSegmentLengths + sumCurveLineSegmentLengths;
        const SimTK::Real gotTotalPathLength = path->getLength(state);
        CHECK_THAT(gotTotalPathLength, Catch::Matchers::WithinRel(
                expectedTotalPathLength, lengthTolerance));
    }

    SECTION("Curve segment initial and final (Frenet) frames") {
        const SimTK::Real frenetFrameTolerance = 1e-6;
        auto assertCurveSegmentFrenetFrames = [&](
                SimTK::CableSpanObstacleIndex obsIx,
                const SimTK::Transform& expected_X_GP,
                const SimTK::Transform& expected_X_GQ)
        {
            const SimTK::Transform& got_X_GP =
                path->calcCurveSegmentInitialFrenetFrame(state, obsIx);
            const SimTK::Transform& got_X_GQ =
                path->calcCurveSegmentFinalFrenetFrame(state, obsIx);

            CHECK_THAT((expected_X_GP.p() - got_X_GP.p()).norm(),
                Catch::Matchers::WithinAbs(0., frenetFrameTolerance));
            CHECK_THAT(
                (expected_X_GP.R().asMat33() - got_X_GP.R().asMat33()).norm(),
                Catch::Matchers::WithinAbs(0., frenetFrameTolerance));
            CHECK_THAT(
                (expected_X_GQ.p() - got_X_GQ.p()).norm(),
                Catch::Matchers::WithinAbs(0., frenetFrameTolerance));
            CHECK_THAT(
                (expected_X_GQ.R().asMat33() - got_X_GQ.R().asMat33()).norm(),
                Catch::Matchers::WithinAbs(0., frenetFrameTolerance));
        };

        assertCurveSegmentFrenetFrames(
            SimTK::CableSpanObstacleIndex(0),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(0., 1., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(-1., 0., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(0., 1., 0.)),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(1., 0., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(0., 1., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(1., 2., 0.)));
        assertCurveSegmentFrenetFrames(
            SimTK::CableSpanObstacleIndex(1),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(1., 0., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(0., 1., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(4.5, 2., 0.)),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(0., -1., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(1., 0., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(5.5, 1., 0.)));
        assertCurveSegmentFrenetFrames(
            SimTK::CableSpanObstacleIndex(2),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(0., -1., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(1., 0., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(5.5, -2., 0.)),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(-1., 0., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(0., -1., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(4., -3.5, 0.)));
        assertCurveSegmentFrenetFrames(
            SimTK::CableSpanObstacleIndex(3),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(-1., 0., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(0., -1., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(-2., -3.5, 0.)),
            SimTK::Transform(
                SimTK::Rotation().setRotationFromTwoAxes(
                    SimTK::UnitVec3(SimTK::Vec3(0., 1., 0.)),
                    SimTK::XAxis,
                    SimTK::UnitVec3(SimTK::Vec3(-1., 0., 0.)),
                    SimTK::YAxis),
                SimTK::Vec3(-4., -1.5, 0.)));
    }
}

// Simple Scholz2015GeometryPath comprised of via points with a known solution.
//
// This test is based on "testViaPoints" in Simbody's TestCableSpan.cpp.
TEST_CASE("Via points") {

    // Create an empty model and path object.
    Model model;
    Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
    path->setOrigin(model.getGround(), SimTK::Vec3(0.));
    path->setInsertion(model.getGround(), SimTK::Vec3(-4., 0., 0.));

    // Add a sequence of via points to the path.
    path->addViaPoint(model.getGround(), SimTK::Vec3(0., 1., 0.));
    path->addViaPoint(model.getGround(), SimTK::Vec3(0., 1., 2.));
    path->addViaPoint(model.getGround(), SimTK::Vec3(-4., 1., 2.));

    // Add the path to the model.
    model.addComponent(path);

    // Initialize the system and state.
    SimTK::State state = model.initSystem();
    model.realizePosition(state);

    CHECK(path->getNumObstacles() == 0);
    CHECK(path->getNumViaPoints() == 3);

    // With no obstacles, the smoothness should always be zero.
    CHECK(path->getSmoothness(state) == 0.);

    // Sum all straight line segment lengths.
    const SimTK::Real lengthTolerance = 1e-5;
    const SimTK::Real expectedTotalCableLength = 1. + 2. + 4. + std::sqrt(5.);
    const SimTK::Real gotTotalCableLength = path->getLength(state);
    CHECK_THAT(gotTotalCableLength, Catch::Matchers::WithinRel(
            expectedTotalCableLength, lengthTolerance));
}


// Simulate a Scholz2015GeometryPath over all supported surfaces and a via
// point, testing that we get the correct path kinematics. The path wraps over
// the following obstacles, in order: torus, ellipsoid, via point, sphere,
// cylinder, torus.
//
// This test is based on "testAllSurfaceKinds" in Simbody's TestCableSpan.cpp.
TEST_CASE("Path with all surfaces and a via point") {
    // Create an empty model.
    Model model;

    // Origin body and joint.
    auto* originBody = new Body("origin_body", 1.0, SimTK::Vec3(0),
        SimTK::Inertia(1.0));
    model.addComponent(originBody);
    auto* originJoint = new FreeJoint("origin_joint",
        model.getGround(), SimTK::Vec3(-8., 0.1, 0.), SimTK::Vec3(0),
        *originBody, SimTK::Vec3(0), SimTK::Vec3(0));
    model.addComponent(originJoint);

    // Via point body and joint.
    auto* viaPointBody = new Body("via_point_body", 1.0, SimTK::Vec3(0),
        SimTK::Inertia(1.0));
    model.addComponent(viaPointBody);
    auto* viaPointJoint = new FreeJoint("via_point_joint",
        model.getGround(), SimTK::Vec3(0., 0.9, 0.5), SimTK::Vec3(0),
        *viaPointBody, SimTK::Vec3(0), SimTK::Vec3(0));
    model.addComponent(viaPointJoint);

    // Termination body and joint.
    auto* terminationBody = new Body("termination_body", 1.0, SimTK::Vec3(0),
        SimTK::Inertia(1.0));
    model.addComponent(terminationBody);
    auto* terminationJoint = new FreeJoint("termination_joint",
        model.getGround(), SimTK::Vec3(20., 1.0, -1.), SimTK::Vec3(0),
        *terminationBody, SimTK::Vec3(0), SimTK::Vec3(0));
    model.addComponent(terminationJoint);

    // Create the path object.
    Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
    path->setOrigin(*originBody, SimTK::Vec3(0.));
    path->setInsertion(*terminationBody, SimTK::Vec3(0.));
    path->setCurveSegmentAccuracy(1e-9);
    path->setSmoothnessTolerance(1e-4);

    // Add the first torus obstacle.
    auto* torus1 = new ContactTorus(1., 0.2,
        SimTK::Vec3(-4., 0., 0.), SimTK::Vec3(0., 0.5*SimTK::Pi, 0.),
        model.getGround());
    torus1->setName("torus1");
    model.addComponent(torus1);
    path->addObstacle(*torus1, SimTK::Vec3(0.1, 0.2, 0.));

    // Add the ellipsoid obstacle.
    auto* ellipsoid = new ContactEllipsoid(SimTK::Vec3(1.5, 2.6, 1.),
        SimTK::Vec3(-2., 0., 0.), SimTK::Vec3(0), model.getGround());
    ellipsoid->setName("ellipsoid");
    model.addComponent(ellipsoid);
    path->addObstacle(*ellipsoid, SimTK::Vec3(0., 0., 1.1));

    // Add the via point.
    path->addViaPoint(*viaPointBody, SimTK::Vec3(0.));

    // Add the sphere obstacle.
    auto* sphere = new ContactSphere(1.,
        SimTK::Vec3(2., 0., 0.), model.getGround());
    sphere->setName("sphere");
    model.addComponent(sphere);
    path->addObstacle(*sphere, SimTK::Vec3(0.1, 1.1, 0.));

    // Add the cylinder obstacle.
    auto* cylinder = new ContactCylinder(1.,
        SimTK::Vec3(5., 0., 0.), SimTK::Vec3(0.5*SimTK::Pi, 0., 0.),
        model.getGround());
    cylinder->setName("cylinder");
    model.addComponent(cylinder);
    path->addObstacle(*cylinder, SimTK::Vec3(0., -1., 0.));

    // Add the second torus obstacle.
    auto* torus2 = new ContactTorus(1., 0.2,
        SimTK::Vec3(14., 0., 0.), SimTK::Vec3(0., 0.5*SimTK::Pi, 0.),
        model.getGround());
    torus2->setName("torus2");
    model.addComponent(torus2);
    path->addObstacle(*torus2, SimTK::Vec3(0.1, 0.2, 0.));

    // Add the path to the model.
    model.addComponent(path);

    // Initialize the system and state.
    SimTK::State state = model.initSystem();
    model.realizePosition(state);

    // Use this to assert cable length time derivative.
    SimTK::Real prevCableLength = SimTK::NaN;

    auto setJointKinematics =
            [&](FreeJoint& joint, SimTK::Vec3 q, SimTK::Vec3 u)
    {
        auto& x_coord = joint.updCoordinate(FreeJoint::Coord::TranslationX);
        x_coord.setValue(state, q[0]);
        x_coord.setSpeedValue(state, u[0]);

        auto& y_coord = joint.updCoordinate(FreeJoint::Coord::TranslationY);
        y_coord.setValue(state, q[1]);
        y_coord.setSpeedValue(state, u[1]);

        auto& z_coord = joint.updCoordinate(FreeJoint::Coord::TranslationZ);
        z_coord.setValue(state, q[2]);
        z_coord.setSpeedValue(state, u[2]);
    };

    // Let the path end and via points be parameterized by an angle, and draw
    // the path for different angles.
    const SimTK::Real dAngle = 1e-4;
    const SimTK::Real finalAngle = 0.5 * SimTK::Pi;
    for (SimTK::Real angle = 0.; angle < finalAngle; angle += dAngle) {

        // Move the cable origin.
        setJointKinematics(*originJoint,
                SimTK::Vec3(1.1 * sin(angle),
                            5. * sin(angle * 1.5),
                            5. * sin(angle * 2.)),
                SimTK::Vec3(1.1 * cos(angle),
                            5. * 1.5 * cos(angle * 1.5),
                            5. * 2. * cos(angle * 2.)));

        // Move the via point.
        setJointKinematics(*viaPointJoint,
                SimTK::Vec3(0., 0.5 * cos(angle), 0.),
                SimTK::Vec3(0., 0.5 * -sin(angle), 0.));

        // Move the cable termination.
        setJointKinematics(*terminationJoint,
                SimTK::Vec3(0.1 * sin(angle),
                            4. * sin(angle * 0.7),
                            10. * sin(angle * 1.3)),
                SimTK::Vec3(0.1 * cos(angle),
                            4. * 0.7 * cos(angle * 0.7),
                            10. * 1.3 * cos(angle * 1.3)));

        // Compute the path.
        model.realizeVelocity(state);
        const SimTK::Real cableLength = path->getLength(state);

        // Assert length derivative using the change in length.
        if (!SimTK::isNaN(prevCableLength)) {
            const SimTK::Real tolerance = 5e-3;
            const SimTK::Real expectedCableLengtheningSpeed =
                    (cableLength - prevCableLength) / dAngle;
            const SimTK::Real gotCableLengtheningSpeed =
                    path->getLengtheningSpeed(state);
            CHECK_THAT(gotCableLengtheningSpeed,
                Catch::Matchers::WithinAbs(expectedCableLengtheningSpeed,
                    tolerance));
        }

        // Total path length should be longer than direct distance between the
        // path endpoints.
        const SimTK::Real distanceBetweenEndPoints =
            (terminationBody->getPositionInGround(state) -
             originBody->getPositionInGround(state))
                .norm();
        CHECK(cableLength > distanceBetweenEndPoints);

        // Make sure that we actually solved the path up to tolerance.
        CHECK(path->getSmoothness(state) <= path->getSmoothnessTolerance());

        prevCableLength = cableLength;
    }
}

// Simulate touchdown and liftoff events in a Scholz2015GeometryPath. A simple
// touchdown and liftoff case on torus, ellipsoid, sphere, and cylinder wrapping
// obstacles. A path is spanned over each obstacle individually, so there are 4
// paths, each with one obstacle.
//
// This test is based on "testTouchdownAndLiftoff" in Simbody's
// TestCableSpan.cpp.
TEST_CASE("Touchdown and liftoff") {
    // Create an empty model.
    Model model;

    // Helper for creating a path with a single obstacle at a certain offset
    // location.
    SimTK::Vec3 originShift(-2., 0., 0.);
    SimTK::Vec3 terminationShift(2., 0., 0.);
    auto createPath = [&](
        const std::string& name, const ContactGeometry& obstacle,
        const SimTK::Vec3& sceneOffset) {

        auto* originBody = new Body(
            fmt::format("{}_origin_body", name), 1.0, SimTK::Vec3(0),
            SimTK::Inertia(1.0));
        model.addComponent(originBody);
        auto* originJoint = new FreeJoint(fmt::format("{}_origin_joint", name),
            model.getGround(), sceneOffset + originShift, SimTK::Vec3(0),
            *originBody, SimTK::Vec3(0), SimTK::Vec3(0));
        model.addComponent(originJoint);

        // Termination body and joint.
        auto* terminationBody = new Body(
            fmt::format("{}_termination_body", name), 1.0, SimTK::Vec3(0),
            SimTK::Inertia(1.0));
        model.addComponent(terminationBody);
        auto* terminationJoint = new FreeJoint(
            fmt::format("{}_termination_joint", name),
            model.getGround(), sceneOffset + terminationShift, SimTK::Vec3(0),
            *terminationBody, SimTK::Vec3(0), SimTK::Vec3(0));
        model.addComponent(terminationJoint);

        // Create the path object.
        Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
        path->setName(fmt::format("{}_path", name));
        path->setOrigin(*originBody, SimTK::Vec3(0.));
        path->setInsertion(*terminationBody, SimTK::Vec3(0.));
        path->setCurveSegmentAccuracy(1e-12);
        path->setSmoothnessTolerance(1e-6);

        // Add the obstacle.
        path->addObstacle(obstacle, SimTK::Vec3(SimTK::NaN));

        // Add the path to the model.
        model.addComponent(path);
    };

    // Create a path with a torus obstacle.
    SimTK::Vec3 sceneOffset(0., 2., 0.);
    auto* torus = new ContactTorus(2., 0.25,
        SimTK::Vec3(0., 1.75, 0.) + sceneOffset,
        SimTK::Vec3(0., 0.5*SimTK::Pi, 0.),
        model.getGround());
    torus->setName("torus");
    model.addComponent(torus);
    createPath("torus", *torus, sceneOffset);

    // Create a path with an ellipsoid obstacle.
    sceneOffset = SimTK::Vec3(0., 0., 0.);
    auto* ellipsoid = new ContactEllipsoid(SimTK::Vec3(1., 0.5, 0.75),
        SimTK::Vec3(0., -0.5, 0.) + sceneOffset,
        SimTK::Vec3(0., 0., 0.),
        model.getGround());
    ellipsoid->setName("ellipsoid");
    model.addComponent(ellipsoid);
    createPath("ellipsoid", *ellipsoid, sceneOffset);

    // Create a path with a sphere obstacle.
    sceneOffset = SimTK::Vec3(0., -2., 0.);
    auto* sphere = new ContactSphere(1.5,
        SimTK::Vec3(0., -1.5, 0.) + sceneOffset,
        model.getGround());
    sphere->setName("sphere");
    model.addComponent(sphere);
    createPath("sphere", *sphere, sceneOffset);

    // Create a path with a cylinder obstacle.
    sceneOffset = SimTK::Vec3(0., -6., 0.);
    auto* cylinder = new ContactCylinder(2.,
        SimTK::Vec3(0., -2., 0.) + sceneOffset,
        SimTK::Vec3(0),
        model.getGround());
    cylinder->setName("cylinder");
    model.addComponent(cylinder);
    createPath("cylinder", *cylinder, sceneOffset);

    // Initialize the system and state.
    SimTK::State state = model.initSystem();
    const SimTK::SimbodyMatterSubsystem& matter = model.getMatterSubsystem();
    const SimTK::CableSubsystem& cables = model.getCableSubsystem();
    CHECK(cables.getNumCables() == 4);

    // Simulate touchdown and liftoff events.
    for (SimTK::Real angle = 1e-2; angle < 4. * SimTK::Pi; angle += 0.02) {
        // Move the cable end points.
        const SimTK::Real yCoord = 0.1 * sin(angle);
        for (SimTK::CableSpanIndex cableIx(0); cableIx < cables.getNumCables();
             ++cableIx) {
            matter
                .getMobilizedBody(cables.getCable(cableIx).getOriginBodyIndex())
                .setQToFitTranslation(state, SimTK::Vec3(0., yCoord, 0.));
            matter
                .getMobilizedBody(
                    cables.getCable(cableIx).getTerminationBodyIndex())
                .setQToFitTranslation(state, SimTK::Vec3(0., yCoord, 0.));
        }

        model.realizePosition(state);

        // All obstacles are positioned such that for negative yCoord of the
        // endpoints, the cable touches down on the obstacle.
        for (SimTK::CableSpanIndex cableIx(0); cableIx < cables.getNumCables();
             ++cableIx) {
            cables.getCable(cableIx).calcLength(state);
            const bool gotContactStatus =
                cables.getCable(cableIx).isInContactWithObstacle(
                    state,
                    SimTK::CableSpanObstacleIndex(0));
            const bool expectedContactStatus = yCoord < 0.;
            CHECK(gotContactStatus == expectedContactStatus);
        }
    }
}


// In this test the optimal path is far from the initial path. This tests the
// robustness of the algorithm to the initial conditions.
//
// This test is based on "testRobustInitialPath" in Simbody's TestCableSpan.cpp.
TEST_CASE("Robust initial path") {
    // Create an empty model.
    Model model;

    // Create the path object.
    Scholz2015GeometryPath* path = new Scholz2015GeometryPath();
    path->setOrigin(model.getGround(), SimTK::Vec3(-0.1, 0., 0.));
    path->setInsertion(model.getGround(), SimTK::Vec3(0.1, 0., 0.));
    path->setCurveSegmentAccuracy(1e-12);
    path->setSmoothnessTolerance(1e-8);
    path->setSolverMaxIterations(100);

    // Set the algorithm to 'MinimumLength'. This test will fail if set to the
    // default algorithm, 'Scholz2015'.
    path->setAlgorithm("MinimumLength");

    // Add a torus obstacle
    auto* torus = new ContactTorus(2., 0.5,
        SimTK::Vec3(0., 4., 0.), SimTK::Vec3(0., 0.5*SimTK::Pi, 0.),
        model.getGround());
    torus->setName("torus");
    model.addComponent(torus);
    path->addObstacle(*torus, SimTK::Vec3(0.0, -0.1, 0.));

    // Add the path to the model.
    model.addComponent(path);

    // Initialize the system and state.
    SimTK::State state = model.initSystem();

    // Realize to the position stage.
    model.realizePosition(state);

    // Test that the solution was found.
    CHECK(path->getSmoothness(state) < path->getSmoothnessTolerance());

    // Test the cable length.
    CHECK_THAT(path->getLength(state),
        Catch::Matchers::WithinRel(5.6513, 1e-3));
}
