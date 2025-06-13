/* -------------------------------------------------------------------------- *
 * OpenSim Moco: sandboxSandbox.cpp                                           *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2019 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Christopher Dembia                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0          *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

// This file provides a way to easily prototype or test temporary snippets of
// code during development.

#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/Actuators/DeGrooteFregly2016Muscle.h>
#include <OpenSim/Actuators/HyfydyMuscle.h>
#include <OpenSim/Actuators/Millard2012EquilibriumMuscle.h>

using namespace OpenSim;

enum MuscleType {
    Hyfydy,
    DeGrooteFregly2016,
    Millard2012Equilibrium,
};

Model createHangingMuscleModel(MuscleType type, 
            double optimalFiberLength, double tendonSlackLength) {
    Model model;
    model.setName("isometric_muscle");
    model.set_gravity(SimTK::Vec3(9.81, 0, 0));
    auto* body = new Body("body", 0.5, SimTK::Vec3(0), SimTK::Inertia(0));
    model.addComponent(body);

    // Allows translation along x.
    auto* joint = new SliderJoint("joint", model.getGround(), *body);
    auto& coord = joint->updCoordinate(SliderJoint::Coord::TranslationX);
    coord.setName("height");
    model.addComponent(joint);

    if (type == MuscleType::DeGrooteFregly2016) {
        auto* actu = new DeGrooteFregly2016Muscle();
        actu->setName("muscle");
        actu->set_max_isometric_force(1000.0);
        actu->set_optimal_fiber_length(optimalFiberLength);
        actu->set_tendon_slack_length(tendonSlackLength);
        actu->set_ignore_activation_dynamics(false);
        actu->set_ignore_tendon_compliance(true);
        actu->set_fiber_damping(0.01);
        actu->set_max_contraction_velocity(10);
        actu->set_pennation_angle_at_optimal(0);
        actu->addNewPathPoint("origin", model.updGround(), SimTK::Vec3(0));
        actu->addNewPathPoint("insertion", *body, SimTK::Vec3(0));
        model.addForce(actu);

    } else if (type == MuscleType::Hyfydy) {
        auto* actu = new HyfydyMuscle();
        actu->setName("muscle");
        actu->set_max_isometric_force(1000.0);
        actu->set_optimal_fiber_length(optimalFiberLength);
        actu->set_tendon_slack_length(tendonSlackLength);
        actu->set_ignore_activation_dynamics(false);
        actu->set_ignore_tendon_compliance(true);
        actu->set_fiber_damping(0.01);
        actu->set_max_contraction_velocity(10);
        actu->set_pennation_angle_at_optimal(0);
        actu->addNewPathPoint("origin", model.updGround(), SimTK::Vec3(0));
        actu->addNewPathPoint("insertion", *body, SimTK::Vec3(0));
        model.addForce(actu);
        
    } else if (type == MuscleType::Millard2012Equilibrium) {
        auto* actu = new Millard2012EquilibriumMuscle();
        actu->setName("muscle");
        actu->set_max_isometric_force(1000.0);
        actu->set_optimal_fiber_length(optimalFiberLength);
        actu->set_tendon_slack_length(tendonSlackLength);
        actu->set_ignore_activation_dynamics(false);
        actu->set_ignore_tendon_compliance(true);
        actu->set_fiber_damping(0.01);
        actu->set_max_contraction_velocity(10);
        actu->set_pennation_angle_at_optimal(0);
        actu->addNewPathPoint("origin", model.updGround(), SimTK::Vec3(0));
        actu->addNewPathPoint("insertion", *body, SimTK::Vec3(0));
        model.addForce(actu);
        
    } else {
        throw Exception("Unknown muscle type.");
    }
    body->attachGeometry(new Sphere(0.05));

    return model;
}

int main() {

    Logger::setLevel(Logger::Level::Error);

    const double optimalFiberLength = 0.1;
    const double tendonSlackLength = 0.05;
    SimTK::Real initHeight = 0.165;
    SimTK::Real finalHeight = 0.155;
    SimTK::Real initSpeed = 0;
    SimTK::Real finalSpeed = 0;
    MocoBounds heightBounds(0.14, 0.17);
    MocoBounds speedBounds(-10, 10);
    MocoBounds actuBounds(0.1, 1);
    int N = 1000000;

    SimTK::Random::Uniform randomHeight(
        heightBounds.getLower(), heightBounds.getUpper());
    SimTK::Random::Uniform randomSpeed(
        speedBounds.getLower(), speedBounds.getUpper());
    SimTK::Random::Uniform randomActu(
        actuBounds.getLower(), actuBounds.getUpper());

    // randomize state lambda function
    auto randomizeState = [&](SimTK::State& s) {
        s.updQ()[0] = randomHeight.getValue();
        s.updU()[0] = randomSpeed.getValue();
        s.updZ()[0] = randomActu.getValue();
    };

    std::cout << "DeGrooteFregly2016Muscle" << std::endl;
    std::cout << "------------------------" << std::endl;
    Model dgfModel = createHangingMuscleModel(
        MuscleType::DeGrooteFregly2016, optimalFiberLength, tendonSlackLength);
    SimTK::State s = dgfModel.initSystem();
    SimTK::Real realTime = SimTK::realTime();
    for (int i = 0; i < N; ++i) {
        randomizeState(s);
        dgfModel.realizeDynamics(s);
    }
    SimTK::Real dgfTime = SimTK::realTime() - realTime;
    std::cout << "Real time for " << N << " realizations: " << dgfTime
              << " seconds." << std::endl << std::endl;
    
    std::cout << "HyfydyMuscle" << std::endl;
    std::cout << "------------" << std::endl;
    Model hyfydyModel = createHangingMuscleModel(
        MuscleType::Hyfydy, optimalFiberLength, tendonSlackLength);
    s = hyfydyModel.initSystem();
    realTime = SimTK::realTime();
    for (int i = 0; i < N; ++i) {
        randomizeState(s);
        hyfydyModel.realizeDynamics(s);
    }
    SimTK::Real hyfydyTime = SimTK::realTime() - realTime;
    std::cout << "Real time for " << N << " realizations: " << hyfydyTime
              << " seconds." << std::endl << std::endl;

    std::cout << "Millard2012EquilibriumMuscle" << std::endl;
    std::cout << "----------------------------" << std::endl;
    Model millardModel = createHangingMuscleModel(
        MuscleType::Millard2012Equilibrium, optimalFiberLength, tendonSlackLength);
    s = millardModel.initSystem();
    realTime = SimTK::realTime();
    for (int i = 0; i < N; ++i) {
        randomizeState(s);
        millardModel.realizeDynamics(s);
    }
    SimTK::Real millardTime = SimTK::realTime() - realTime;
    std::cout << "Real time for " << N << " realizations: " << millardTime
              << " seconds." << std::endl << std::endl;

    std::cout << "HyfydyMuscles is " << (hyfydyTime < dgfTime ? "faster" : "slower")
              << " than DeGrooteFregly2016Muscles by a factor of "
              << (dgfTime / hyfydyTime) << "." << std::endl;
    std::cout << "Millard2012EquilibriumMuscles is "
              << (millardTime < dgfTime ? "faster" : "slower")
              << " than DeGrooteFregly2016Muscles by a factor of "
              << (dgfTime / millardTime) << "." << std::endl;
    std::cout << "HyfydyMuscles is "
              << (hyfydyTime < millardTime ? "faster" : "slower")
              << " than Millard2012EquilibriumMuscles by a factor of "
              << (millardTime / hyfydyTime) << "." << std::endl;
    std::cout << "Millard2012EquilibriumMuscles is "
              << (millardTime < hyfydyTime ? "faster" : "slower")
              << " than HyfydyMuscles by a factor of "
              << (hyfydyTime / millardTime) << "." << std::endl;
    std::cout << "DeGrooteFregly2016Muscles is "
              << (dgfTime < hyfydyTime ? "faster" : "slower")
              << " than HyfydyMuscles by a factor of "
              << (hyfydyTime / dgfTime) << "." << std::endl;
    std::cout << "DeGrooteFregly2016Muscles is "
              << (dgfTime < millardTime ? "faster" : "slower")
              << " than Millard2012EquilibriumMuscles by a factor of "
              << (millardTime / dgfTime) << "." << std::endl;

    return EXIT_SUCCESS;
}
