/* -------------------------------------------------------------------------- *
 *          OpenSim:  testScholz2015GeometryPathMovingPoints.cpp              *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2026 Stanford University and the Authors                *
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

#include <OpenSim/Common/ExpressionBasedFunction.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PathSpring.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>
#include <OpenSim/Simulation/SimbodyEngine/CantileverFreeBeamJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>

#include <catch2/catch_all.hpp>

using namespace OpenSim;

namespace {

const std::vector<std::string> c_lumbarCoordinates = {
        "/jointset/lumbar/lumbar_coord_0",
        "/jointset/lumbar/lumbar_coord_1",
        "/jointset/lumbar/lumbar_coord_2"};

// A pelvis welded to ground, a spine attached to the pelvis with a
// CantileverFreeBeamJoint (whose parent frame is offset and rotated relative
// to the pelvis), and a thigh attached to the pelvis with a pin joint.
Model createModel() {
    Model model;
    model.setName("moving_points_model");
    auto* pelvis = new Body("pelvis", 1.0, SimTK::Vec3(0), SimTK::Inertia(1));
    auto* spine = new Body("spine", 1.0, SimTK::Vec3(0), SimTK::Inertia(1));
    auto* thigh = new Body("thigh", 1.0, SimTK::Vec3(0), SimTK::Inertia(1));
    model.addBody(pelvis);
    model.addBody(spine);
    model.addBody(thigh);
    model.addJoint(new WeldJoint("ground_pelvis", model.getGround(), *pelvis));
    auto* lumbar = new CantileverFreeBeamJoint("lumbar",
            *pelvis, SimTK::Vec3(0.02, 0.1, -0.01),
                     SimTK::Vec3(-SimTK::Pi/2, 0.1, 0.2),
            *spine, SimTK::Vec3(0), SimTK::Vec3(0), 0.3);
    for (int i = 0; i < 3; ++i) {
        lumbar->upd_coordinates(i).setName("lumbar_coord_" + std::to_string(i));
    }
    model.addJoint(lumbar);
    model.addJoint(new PinJoint("hip",
            *pelvis, SimTK::Vec3(0.1, -0.1, 0.05), SimTK::Vec3(0),
            *thigh, SimTK::Vec3(0, 0.2, 0), SimTK::Vec3(0)));
    model.finalizeConnections();
    return model;
}

Scholz2015GeometryPath& addPathSpring(Model& model) {
    auto* spring = new PathSpring();
    spring->setName("spring");
    spring->setStiffness(1.0);
    spring->setRestingLength(0.0);
    spring->setDissipation(0.0);
    spring->set_path(Scholz2015GeometryPath());
    model.addComponent(spring);
    return spring->updPath<Scholz2015GeometryPath>();
}

ExpressionBasedFunction expression(const std::string& expr,
        const std::vector<std::string>& variables) {
    return ExpressionBasedFunction(expr, variables);
}

// A path with a moving origin on the lumbar joint's (offset) parent frame,
// driven by the lumbar coordinates; a moving via point on the spine, driven by
// the hip coordinate; and a fixed insertion on the thigh.
Scholz2015GeometryPath& addMovingPath(Model& model) {
    auto& path = addPathSpring(model);
    const std::vector<std::string> abc = {"a", "b", "c"};
    path.appendPathPoint(
            model.getComponent<CantileverFreeBeamJoint>("/jointset/lumbar")
                    .getParentFrame(),
            expression("0.03 + 0.02*sin(a) + 0.01*b*c", abc),
            expression("-0.02 + 0.03*cos(b) - 0.01*a", abc),
            expression("0.05 + 0.02*a*a + 0.01*sin(c)", abc),
            c_lumbarCoordinates);
    path.appendPathPoint(model.getComponent<Body>("/bodyset/spine"),
            expression("0.04 + 0.02*sin(h)", {"h"}),
            expression("0.1 + 0.01*h*h", {"h"}),
            expression("0.02*cos(h)", {"h"}),
            {"/jointset/hip/hip_coord_0"});
    path.appendPathPoint(model.getComponent<Body>("/bodyset/thigh"),
            SimTK::Vec3(0.05, -0.1, 0.02));
    return path;
}

void setCoordinates(const Model& model, SimTK::State& state,
        const SimTK::Vector& q, const SimTK::Vector& u) {
    const auto& coords = model.getCoordinateSet();
    for (int i = 0; i < coords.getSize(); ++i) {
        coords.get(i).setValue(state, q[i], false);
        coords.get(i).setSpeedValue(state, u[i]);
    }
    model.realizeVelocity(state);
}

const SimTK::Vector c_q = SimTK::Vector(SimTK::Vec4(0.3, -0.4, 0.2, 0.5));
const SimTK::Vector c_u = SimTK::Vector(SimTK::Vec4(0.7, -1.1, 0.5, 0.9));

} // namespace

TEST_CASE("Constant location functions match fixed path points") {
    // A moving point on the joint's offset parent frame F with constant
    // functions must match a fixed point on the pelvis at X_BF * location.
    const SimTK::Vec3 location(0.03, -0.02, 0.05);

    Model movingModel = createModel();
    auto& movingPath = addPathSpring(movingModel);
    const auto& frame = movingModel.getComponent<CantileverFreeBeamJoint>(
            "/jointset/lumbar").getParentFrame();
    auto constant = [](double value) {
        return expression(std::to_string(value) + " + 0*a", {"a"});
    };
    movingPath.appendPathPoint(frame, constant(location[0]),
            constant(location[1]), constant(location[2]),
            {c_lumbarCoordinates[0]});
    movingPath.appendPathPoint(movingModel.getComponent<Body>("/bodyset/thigh"),
            SimTK::Vec3(0.05, -0.1, 0.02));
    SimTK::State movingState = movingModel.initSystem();
    setCoordinates(movingModel, movingState, c_q, c_u);

    Model fixedModel = createModel();
    auto& fixedPath = addPathSpring(fixedModel);
    const SimTK::Vec3 locationInPelvis =
            fixedModel.getComponent<CantileverFreeBeamJoint>("/jointset/lumbar")
                    .getParentFrame().findTransformInBaseFrame() * location;
    fixedPath.appendPathPoint(fixedModel.getComponent<Body>("/bodyset/pelvis"),
            locationInPelvis);
    fixedPath.appendPathPoint(fixedModel.getComponent<Body>("/bodyset/thigh"),
            SimTK::Vec3(0.05, -0.1, 0.02));
    SimTK::State fixedState = fixedModel.initSystem();
    setCoordinates(fixedModel, fixedState, c_q, c_u);

    CHECK_THAT(movingPath.getLength(movingState),
            Catch::Matchers::WithinAbs(fixedPath.getLength(fixedState), 1e-12));
    CHECK_THAT(movingPath.getLengtheningSpeed(movingState),
            Catch::Matchers::WithinAbs(
                    fixedPath.getLengtheningSpeed(fixedState), 1e-12));
}

TEST_CASE("Lengthening speed and moment arms are consistent with length") {
    Model model = createModel();
    auto& path = addMovingPath(model);
    SimTK::State state = model.initSystem();
    setCoordinates(model, state, c_q, c_u);

    auto lengthAt = [&](const SimTK::Vector& q) {
        SimTK::State s = state;
        setCoordinates(model, s, q, c_u);
        return path.getLength(s);
    };

    // Lengthening speed (qdot = u for these coordinates).
    const double h = 1e-6;
    const double lengthDot =
            (lengthAt(c_q + h * c_u) - lengthAt(c_q - h * c_u)) / (2 * h);
    CHECK_THAT(path.getLengtheningSpeed(state),
            Catch::Matchers::WithinAbs(lengthDot, 1e-8));

    // Moment arms from the forces the path produces (MomentArmSolver) must
    // match -dL/dq, which requires the generalized forces of the moving
    // points.
    const auto& coords = model.getCoordinateSet();
    for (int i = 0; i < coords.getSize(); ++i) {
        SimTK::Vector dq(4, 0.0);
        dq[i] = h;
        const double expected = -(lengthAt(c_q + dq) - lengthAt(c_q - dq)) /
                (2 * h);
        CAPTURE(coords.get(i).getName());
        CHECK_THAT(path.computeMomentArm(state, coords.get(i)),
                Catch::Matchers::WithinAbs(expected, 1e-7));
    }

    // The independent coordinates include the lumbar coordinates, even though
    // the path's origin frame (on the pelvis) is not a child of that joint.
    CHECK(path.findIndependentCoordinates(state).size() == 4);
}

TEST_CASE("Serialization") {
    Model model = createModel();
    auto& path = addMovingPath(model);
    SimTK::State state = model.initSystem();
    setCoordinates(model, state, c_q, c_u);
    const double length = path.getLength(state);

    model.print("testScholz2015GeometryPathMovingPoints_model.osim");
    Model copy("testScholz2015GeometryPathMovingPoints_model.osim");
    SimTK::State copyState = copy.initSystem();
    setCoordinates(copy, copyState, c_q, c_u);
    CHECK_THAT(copy.getComponent<PathSpring>("/spring").getPath()
            .getLength(copyState),
            Catch::Matchers::WithinAbs(length, 1e-12));
}

TEST_CASE("Path point accessors") {
    // Fixed and moving points are both path points, with a parent frame and a
    // state-dependent location.
    Model model = createModel();
    auto& path = addMovingPath(model);
    SimTK::State state = model.initSystem();
    setCoordinates(model, state, c_q, c_u);
    CHECK(path.getNumPathPoints() == 3);

    const auto& origin = path.getOrigin();
    REQUIRE(dynamic_cast<const Scholz2015GeometryPathMovingPoint*>(&origin));
    CHECK(origin.getParentFrame().getAbsolutePathString() ==
            "/jointset/lumbar/pelvis_offset");
    const double a = c_q[0], b = c_q[1], c = c_q[2];
    const SimTK::Vec3 expected(0.03 + 0.02*std::sin(a) + 0.01*b*c,
            -0.02 + 0.03*std::cos(b) - 0.01*a,
            0.05 + 0.02*a*a + 0.01*std::sin(c));
    CHECK_THAT((origin.calcLocation(state) - expected).norm(),
            Catch::Matchers::WithinAbs(0., 1e-12));

    CHECK(&path.getPathPoint(0) == &origin);
    CHECK(dynamic_cast<const Scholz2015GeometryPathMovingPoint*>(
            &path.getPathPoint(1)));
    const auto* insertion = dynamic_cast<const Scholz2015GeometryPathPoint*>(
            &path.getInsertion());
    REQUIRE(insertion);
    CHECK(&path.getPathPoint(2) == insertion);
    CHECK(insertion->calcLocation(state) ==
            insertion->getPathPoint().get_location());
}

TEST_CASE("Invalid location functions") {
    Model model = createModel();
    auto& path = addPathSpring(model);
    CHECK_THROWS_WITH(
            path.appendPathPoint(model.getComponent<Body>("/bodyset/spine"),
                    expression("a", {"a"}), expression("a + b", {"a", "b"}),
                    expression("a", {"a"}), {"/jointset/hip/hip_coord_0"}),
            Catch::Matchers::ContainsSubstring(
                    "Expected the location functions to take one argument "
                    "per coordinate (1)"));
}
