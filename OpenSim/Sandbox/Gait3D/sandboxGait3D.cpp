/* -------------------------------------------------------------------------- *
 *                        OpenSim:  sandboxGait3D.cpp                         *
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

#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/Simulation/VisualizerUtilities.h>
#include <OpenSim/Actuators/ModelOperators.h>
#include <OpenSim/Actuators/DeGrooteFregly2016Muscle.h>
#include <OpenSim/Actuators/Millard2012EquilibriumMuscle.h>
#include <OpenSim/Actuators/ActivationCoordinateActuator.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>
#include <OpenSim/Simulation/Model/ExponentialContactForce.h>
#include <OpenSim/Common/STOFileAdapter.h>
#include <OpenSim/Common/LinearFunction.h>

#include <OpenSim/Simulation/Model/CoordinateLinearStop.h>
#include <OpenSim/Simulation/Model/CoordinateLinearDamper.h>
#include <OpenSim/Simulation/Model/CoordinateLinearSpring.h>
#include <OpenSim/Simulation/Control/DiscreteController.h>

using namespace OpenSim;

// Enums
enum BodyType {LeftFoot=0, RightFoot, LeftShank, RightShank, LeftThigh,
               RightThigh, Pelvis, Torso};

enum MuscleType {GlutMed_R=0, AddMag_R, Hamstrings_R, Bifemsh_R, GlutMax_R,
                 Iliopsoas_R, RectFem_R, Vasti_R, Gastroc_R, Soleus_R, TibAnt_R,
                 GlutMed_L, AddMag_L, Hamstrings_L, Bifemsh_L, GlutMax_L,
                 Iliopsoas_L, RectFem_L, Vasti_L, Gastroc_L, Soleus_L, TibAnt_L};

enum Contact {Heel=0, LateralToe, MedialToe};

enum JointType {Custom=0, Ball};

// Data arrays
SimTK::Real massData[] = {1.25, 1.25, 3.7075, 3.7075, 9.3014, 9.3014,
                          11.777, 34.2366};

SimTK::Vec3 inertiaData[] = {SimTK::Vec3(0.0014, 0.0039, 0.0041),
                             SimTK::Vec3(0.0014, 0.0039, 0.0041),
                             SimTK::Vec3(0.0504, 0.0051, 0.0511),
                             SimTK::Vec3(0.0504, 0.0051, 0.0511),
                             SimTK::Vec3(0.1339, 0.0351, 0.1412),
                             SimTK::Vec3(0.1339, 0.0351, 0.1412),
                             SimTK::Vec3(0.1028, 0.0871, 0.0579),
                             SimTK::Vec3(1.4745, 0.7555, 1.4314)};

SimTK::Vec3 leftContactPoints[] = {SimTK::Vec3(-0.085, -0.015, 0.005),
                                   SimTK::Vec3(0.0425, -0.03, -0.041),
                                   SimTK::Vec3(0.085, -0.03, 0.0275)};

SimTK::Vec3 rightContactPoints[] = {SimTK::Vec3(-0.085, -0.015, -0.005),
                                    SimTK::Vec3(0.0425, -0.03, 0.041),
                                    SimTK::Vec3(0.085, -0.03, -0.0275)};

// Template-based muscle factory
template<typename MuscleType>
MuscleType* createMuscle(Model& model, const std::string& name,
                        double maxIsometricForce, double optimalFiberLength,
                        double tendonSlackLength, double pennationAngle) {
    auto* muscle = new MuscleType();
    muscle->setName(name);
    muscle->set_max_isometric_force(maxIsometricForce);
    muscle->set_optimal_fiber_length(optimalFiberLength);
    muscle->set_tendon_slack_length(tendonSlackLength);
    muscle->set_pennation_angle_at_optimal(pennationAngle);
    muscle->set_ignore_tendon_compliance(true);
    muscle->set_ignore_activation_dynamics(true);
    muscle->set_default_activation(0.01);
    model.addForce(muscle);
    return muscle;
}

// Specialization for DeGrooteFregly2016Muscle
template<>
DeGrooteFregly2016Muscle* createMuscle<DeGrooteFregly2016Muscle>(
        Model& model, const std::string& name,
        double maxIsometricForce, double optimalFiberLength,
        double tendonSlackLength, double pennationAngle) {
    auto* muscle = new DeGrooteFregly2016Muscle();
    muscle->setName(name);
    muscle->set_max_isometric_force(maxIsometricForce);
    muscle->set_optimal_fiber_length(optimalFiberLength);
    muscle->set_tendon_slack_length(tendonSlackLength);
    muscle->set_pennation_angle_at_optimal(pennationAngle);
    muscle->set_ignore_tendon_compliance(true);
    muscle->set_ignore_activation_dynamics(true);
    muscle->set_default_activation(0.01);
    model.addForce(muscle);
    return muscle;
}

// Utility functions
void addContactGeometry(Body* body, const SimTK::Vec3& offset,
                       const std::string& name, const double& radius) {
    PhysicalOffsetFrame* offsetFrame =
            new PhysicalOffsetFrame(name, *body, offset);
    offsetFrame->attachGeometry(new Sphere(radius));
    body->addComponent(offsetFrame);
}

void addContact(Model& model, const std::string& name, PhysicalFrame* frame,
               const SimTK::Vec3& location, const SimTK::Transform& transform,
               const SimTK::ExponentialSpringParameters& params) {
    auto* contact = new ExponentialContactForce(transform, *frame, location,
        params);
    contact->setName(name);
    model.addForce(contact);
}

template<typename MuscleType>
void createMuscles(Model& model, Body* pelvis, Body* torso, Body* leftThigh,
                   Body* leftShank, Body* leftFoot, Body* rightThigh,
                   Body* rightShank, Body* rightFoot, bool useObstacles) {

    // Wrap obstacles
    // --------------
    if (useObstacles) {
        auto* glut_max_r_obstacle = new ContactEllipsoid(SimTK::Vec3(0.04, 0.04, 0.1),
            SimTK::Vec3(-0.04, -0.085, 0.09),
            SimTK::Vec3(-0.2, 0.5, 0.),
            *pelvis);
        glut_max_r_obstacle->setName("glut_max_r_obstacle");
        pelvis->addComponent(glut_max_r_obstacle);

        auto* glut_max_l_obstacle = new ContactEllipsoid(SimTK::Vec3(0.04, 0.04, 0.1),
            SimTK::Vec3(-0.04, -0.085, -0.09),
            SimTK::Vec3(0.2, -0.5, 0.),
            *pelvis);
        glut_max_l_obstacle->setName("glut_max_l_obstacle");
        pelvis->addComponent(glut_max_l_obstacle);

        auto* gastroc_r_obstacle = new ContactEllipsoid(SimTK::Vec3(0.065, 0.065, 0.2),
            SimTK::Vec3(-0.00861, 0.05, 0.),
            SimTK::Vec3(2.9672, 0.028, -1.478),
            *rightShank);
        gastroc_r_obstacle->setName("gastroc_r_obstacle");
        rightShank->addComponent(gastroc_r_obstacle);

        auto* gastroc_l_obstacle = new ContactEllipsoid(SimTK::Vec3(0.065, 0.065, 0.2),
            SimTK::Vec3(-0.00861, 0.05, 0.),
            SimTK::Vec3(-2.9672, 0.028, 1.478),
            *leftShank);
        gastroc_l_obstacle->setName("gastroc_l_obstacle");
        leftShank->addComponent(gastroc_l_obstacle);
    }

    // Right leg muscles
    // -----------------
    auto* glut_med_r = createMuscle<MuscleType>(model, "glut_med_r", 2045,
                                               0.0733, 0.066, 0.3578);
    glut_med_r->set_path(Scholz2015GeometryPath());
    auto& glut_med_r_path = glut_med_r->template updPath<Scholz2015GeometryPath>();
    glut_med_r_path.setName("glut_med_r_path");
    glut_med_r_path.setOrigin(*pelvis, SimTK::Vec3(-0.0148, 0.0445, 0.0766));
    glut_med_r_path.setInsertion(*rightThigh, SimTK::Vec3(-0.0258, 0.1642, 0.0527));

    auto* add_mag_r = createMuscle<MuscleType>(model, "add_mag_r", 2268,
                                              0.087, 0.06, 0.0872665);
    add_mag_r->set_path(Scholz2015GeometryPath());
    auto& add_mag_r_path = add_mag_r->template updPath<Scholz2015GeometryPath>();
    add_mag_r_path.setName("add_mag_r_path");
    add_mag_r_path.setOrigin(*pelvis, SimTK::Vec3(-0.0025, -0.1174, 0.0255));
    add_mag_r_path.setInsertion(*rightThigh, SimTK::Vec3(-0.0045, 0.0489, 0.0339));

    auto* hamstrings_r = createMuscle<MuscleType>(model, "hamstrings_r", 2594,
                                                 0.0976, 0.319, 0.2025);
    hamstrings_r->set_path(Scholz2015GeometryPath());
    auto& hamstrings_r_path = hamstrings_r->template updPath<Scholz2015GeometryPath>();
    hamstrings_r_path.setName("hamstrings_r_path");
    hamstrings_r_path.setOrigin(*pelvis, SimTK::Vec3(-0.05526, -0.10257, 0.06944));
    hamstrings_r_path.setInsertion(*rightShank, SimTK::Vec3(-0.021, 0.1467, 0.0343));
    hamstrings_r_path.addViaPoint(*rightShank, SimTK::Vec3(-0.028, 0.1667, 0.02943));

    auto* bifemsh_r = createMuscle<MuscleType>(model, "bifemsh_r", 804,
                                              0.1103, 0.095, 0.2147);
    bifemsh_r->set_path(Scholz2015GeometryPath());
    auto& bifemsh_r_path = bifemsh_r->template updPath<Scholz2015GeometryPath>();
    bifemsh_r_path.setName("bifemsh_r_path");
    bifemsh_r_path.setOrigin(*rightThigh, SimTK::Vec3(0.005, -0.0411, 0.0234));
    bifemsh_r_path.setInsertion(*rightShank, SimTK::Vec3(-0.021, 0.1467, 0.0343));
    bifemsh_r_path.addViaPoint(*rightShank, SimTK::Vec3(-0.028, 0.1667, 0.02943));

    auto* glut_max_r = createMuscle<MuscleType>(model, "glut_max_r", 1944,
                                               0.1569, 0.111, 0.3822);
    glut_max_r->set_path(Scholz2015GeometryPath());
    auto& glut_max_r_path = glut_max_r->template updPath<Scholz2015GeometryPath>();
    glut_max_r_path.setName("glut_max_r_path");
    glut_max_r_path.setOrigin(*pelvis, SimTK::Vec3(-0.0642, 0.0176, 0.0563));
    glut_max_r_path.setInsertion(*rightThigh, SimTK::Vec3(-0.0156, 0.0684, 0.0419));
    if (useObstacles) {
        glut_max_r_path.addObstacle(
            pelvis->getComponent<ContactGeometry>("glut_max_r_obstacle"),
            SimTK::Vec3(-0.04, 0., 0.));
    } else {
        glut_max_r_path.addViaPoint(*pelvis, SimTK::Vec3(-0.0669, -0.052, 0.0914));
        glut_max_r_path.addViaPoint(*rightThigh, SimTK::Vec3(-0.0426, 0.117, 0.0293));
    }

    auto* iliopsoas_r = createMuscle<MuscleType>(model, "iliopsoas_r", 2186,
                                                0.1066, 0.152, 0.2496);
    iliopsoas_r->set_path(Scholz2015GeometryPath());
    auto& iliopsoas_r_path = iliopsoas_r->template updPath<Scholz2015GeometryPath>();
    iliopsoas_r_path.setName("iliopsoas_r_path");
    iliopsoas_r_path.setOrigin(*pelvis, SimTK::Vec3(0.006, 0.0887, 0.0289));
    iliopsoas_r_path.setInsertion(*rightThigh, SimTK::Vec3(-0.0188, 0.1103, 0.0104));
    iliopsoas_r_path.addViaPoint(*pelvis, SimTK::Vec3(0.0407, -0.01, 0.076));
    iliopsoas_r_path.addViaPoint(*rightThigh, SimTK::Vec3(0.033, 0.135, 0.0038));

    auto* rect_fem_r = createMuscle<MuscleType>(model, "rect_fem_r", 1169,
                                               0.0759, 0.3449, 0.2426);
    rect_fem_r->set_path(Scholz2015GeometryPath());
    auto& rect_fem_r_path = rect_fem_r->template updPath<Scholz2015GeometryPath>();
    rect_fem_r_path.setName("rect_fem_r_path");
    rect_fem_r_path.setOrigin(*pelvis, SimTK::Vec3(0.0412, -0.0311, 0.0968));
    rect_fem_r_path.setInsertion(*rightShank, SimTK::Vec3(0.038, 0.2117, 0.0018));
    rect_fem_r_path.addViaPoint(*rightThigh, SimTK::Vec3(0.038, -0.17, 0.004));

    auto* vasti_r = createMuscle<MuscleType>(model, "vasti_r", 4530,
                                            0.0993, 0.1231, 0.0785);
    vasti_r->set_path(Scholz2015GeometryPath());
    auto& vasti_r_path = vasti_r->template updPath<Scholz2015GeometryPath>();
    vasti_r_path.setName("vasti_r_path");
    vasti_r_path.setOrigin(*rightThigh, SimTK::Vec3(0.029, -0.0224, 0.031));
    vasti_r_path.setInsertion(*rightShank, SimTK::Vec3(0.038, 0.2117, 0.0018));
    vasti_r_path.addViaPoint(*rightThigh, SimTK::Vec3(0.038, -0.17, 0.007));

    auto* gastroc_r = createMuscle<MuscleType>(model, "gastroc_r", 2241,
                                              0.051, 0.384, 0.1728);
    gastroc_r->set_path(Scholz2015GeometryPath());
    auto& gastroc_r_path = gastroc_r->template updPath<Scholz2015GeometryPath>();
    gastroc_r_path.setName("gastroc_r_path");
    gastroc_r_path.setOrigin(*rightThigh, SimTK::Vec3(-0.02, -0.218, -0.024));
    gastroc_r_path.setInsertion(*rightFoot, SimTK::Vec3(-0.095, 0.001, -0.0053));
    if (useObstacles) {
        gastroc_r_path.addObstacle(
            rightShank->getComponent<ContactGeometry>("gastroc_r_obstacle"),
            SimTK::Vec3(0., -0.065, 0.));
    }

    auto* soleus_r = createMuscle<MuscleType>(model, "soleus_r", 3549,
                                             0.044, 0.248, 0.4939);
    soleus_r->set_path(Scholz2015GeometryPath());
    auto& soleus_r_path = soleus_r->template updPath<Scholz2015GeometryPath>();
    soleus_r_path.setName("soleus_r_path");
    soleus_r_path.setOrigin(*rightShank, SimTK::Vec3(-0.0024, 0.0334, 0.0071));
    soleus_r_path.setInsertion(*rightFoot, SimTK::Vec3(-0.095, 0.001, -0.0053));

    auto* tib_ant_r = createMuscle<MuscleType>(model, "tib_ant_r", 1579,
                                              0.0683, 0.243, 0.1676);
    tib_ant_r->set_path(Scholz2015GeometryPath());
    auto& tib_ant_r_path = tib_ant_r->template updPath<Scholz2015GeometryPath>();
    tib_ant_r_path.setName("tib_ant_r_path");
    tib_ant_r_path.setOrigin(*rightShank, SimTK::Vec3(0.0179, 0.0243, 0.0115));
    tib_ant_r_path.setInsertion(*rightFoot, SimTK::Vec3(0.0166, -0.0122, -0.0305));
    tib_ant_r_path.addViaPoint(*rightShank, SimTK::Vec3(0.0329, -0.2084, -0.0177));

    // Left leg muscles
    // -----------------
    auto* glut_med_l = createMuscle<MuscleType>(model, "glut_med_l", 2045,
                                               0.0733, 0.066, 0.3578);
    glut_med_l->set_path(Scholz2015GeometryPath());
    auto& glut_med_l_path = glut_med_l->template updPath<Scholz2015GeometryPath>();
    glut_med_l_path.setName("glut_med_l_path");
    glut_med_l_path.setOrigin(*pelvis, SimTK::Vec3(-0.0148, 0.0445, -0.0766));
    glut_med_l_path.setInsertion(*leftThigh, SimTK::Vec3(-0.0258, 0.1642, -0.0527));

    auto* add_mag_l = createMuscle<MuscleType>(model, "add_mag_l", 2268,
                                              0.087, 0.06, 0.0872665);
    add_mag_l->set_path(Scholz2015GeometryPath());
    auto& add_mag_l_path = add_mag_l->template updPath<Scholz2015GeometryPath>();
    add_mag_l_path.setName("add_mag_l_path");
    add_mag_l_path.setOrigin(*pelvis, SimTK::Vec3(-0.0025, -0.1174, -0.0255));
    add_mag_l_path.setInsertion(*leftThigh, SimTK::Vec3(-0.0045, 0.0489, -0.0339));

    auto* hamstrings_l = createMuscle<MuscleType>(model, "hamstrings_l", 2594,
                                                 0.0976, 0.319, 0.2025);
    hamstrings_l->set_path(Scholz2015GeometryPath());
    auto& hamstrings_l_path = hamstrings_l->template updPath<Scholz2015GeometryPath>();
    hamstrings_l_path.setName("hamstrings_l_path");
    hamstrings_l_path.setOrigin(*pelvis, SimTK::Vec3(-0.05526, -0.10257, -0.06944));
    hamstrings_l_path.setInsertion(*leftShank, SimTK::Vec3(-0.021, 0.1467, -0.0343));
    hamstrings_l_path.addViaPoint(*leftShank, SimTK::Vec3(-0.028, 0.1667, -0.02943));

    auto* bifemsh_l = createMuscle<MuscleType>(model, "bifemsh_l", 804,
                                              0.1103, 0.095, 0.2147);
    bifemsh_l->set_path(Scholz2015GeometryPath());
    auto& bifemsh_l_path = bifemsh_l->template updPath<Scholz2015GeometryPath>();
    bifemsh_l_path.setName("bifemsh_l_path");
    bifemsh_l_path.setOrigin(*leftThigh, SimTK::Vec3(0.005, -0.0411, -0.0234));
    bifemsh_l_path.setInsertion(*leftShank, SimTK::Vec3(-0.021, 0.1467, -0.0343));
    bifemsh_l_path.addViaPoint(*leftShank, SimTK::Vec3(-0.028, 0.1667, -0.02943));

    auto* glut_max_l = createMuscle<MuscleType>(model, "glut_max_l", 1944,
                                               0.1569, 0.111, 0.3822);
    glut_max_l->set_path(Scholz2015GeometryPath());
    auto& glut_max_l_path = glut_max_l->template updPath<Scholz2015GeometryPath>();
    glut_max_l_path.setName("glut_max_l_path");
    glut_max_l_path.setOrigin(*pelvis, SimTK::Vec3(-0.0642, 0.0176, -0.0563));
    glut_max_l_path.setInsertion(*leftThigh, SimTK::Vec3(-0.0156, 0.0684, -0.0419));
    if (useObstacles) {
        glut_max_l_path.addObstacle(
            pelvis->getComponent<ContactGeometry>("glut_max_l_obstacle"),
            SimTK::Vec3(-0.04, 0., 0.));
    } else {
        glut_max_l_path.addViaPoint(*pelvis, SimTK::Vec3(-0.0669, -0.052, -0.0914));
        glut_max_l_path.addViaPoint(*leftThigh, SimTK::Vec3(-0.0426, 0.117, -0.0293));
    }

    auto* iliopsoas_l = createMuscle<MuscleType>(model, "iliopsoas_l", 2186,
                                                0.1066, 0.152, 0.2496);
    iliopsoas_l->set_path(Scholz2015GeometryPath());
    auto& iliopsoas_l_path = iliopsoas_l->template updPath<Scholz2015GeometryPath>();
    iliopsoas_l_path.setName("iliopsoas_l_path");
    iliopsoas_l_path.setOrigin(*pelvis, SimTK::Vec3(0.006, 0.0887, -0.0289));
    iliopsoas_l_path.setInsertion(*leftThigh, SimTK::Vec3(-0.0188, 0.1103, -0.0104));
    iliopsoas_l_path.addViaPoint(*pelvis, SimTK::Vec3(0.0407, -0.01, -0.076));
    iliopsoas_l_path.addViaPoint(*leftThigh, SimTK::Vec3(0.033, 0.135, -0.0038));

    auto* rect_fem_l = createMuscle<MuscleType>(model, "rect_fem_l", 1169,
                                               0.0759, 0.3449, 0.2426);
    rect_fem_l->set_path(Scholz2015GeometryPath());
    auto& rect_fem_l_path = rect_fem_l->template updPath<Scholz2015GeometryPath>();
    rect_fem_l_path.setName("rect_fem_l_path");
    rect_fem_l_path.setOrigin(*pelvis, SimTK::Vec3(0.0412, -0.0311, -0.0968));
    rect_fem_l_path.setInsertion(*leftShank, SimTK::Vec3(0.038, 0.2117, -0.0018));
    rect_fem_l_path.addViaPoint(*leftThigh, SimTK::Vec3(0.038, -0.17, -0.004));

    auto* vasti_l = createMuscle<MuscleType>(model, "vasti_l", 4530,
                                            0.0993, 0.1231, 0.0785);
    vasti_l->set_path(Scholz2015GeometryPath());
    auto& vasti_l_path = vasti_l->template updPath<Scholz2015GeometryPath>();
    vasti_l_path.setName("vasti_l_path");
    vasti_l_path.setOrigin(*leftThigh, SimTK::Vec3(0.029, -0.0224, -0.031));
    vasti_l_path.setInsertion(*leftShank, SimTK::Vec3(0.038, 0.2117, -0.0018));
    vasti_l_path.addViaPoint(*leftThigh, SimTK::Vec3(0.038, -0.17, -0.007));

    auto* gastroc_l = createMuscle<MuscleType>(model, "gastroc_l", 2241,
                                              0.051, 0.384, 0.1728);
    gastroc_l->set_path(Scholz2015GeometryPath());
    auto& gastroc_l_path = gastroc_l->template updPath<Scholz2015GeometryPath>();
    gastroc_l_path.setName("gastroc_l_path");
    gastroc_l_path.setOrigin(*leftThigh, SimTK::Vec3(-0.02, -0.218, 0.024));
    gastroc_l_path.setInsertion(*leftFoot, SimTK::Vec3(-0.095, 0.001, 0.0053));
    if (useObstacles) {
        gastroc_l_path.addObstacle(
            leftShank->getComponent<ContactGeometry>("gastroc_l_obstacle"),
            SimTK::Vec3(-0.065, 0., 0.));
    }

    auto* soleus_l = createMuscle<MuscleType>(model, "soleus_l", 3549,
                                             0.044, 0.248, 0.4939);
    soleus_l->set_path(Scholz2015GeometryPath());
    auto& soleus_l_path = soleus_l->template updPath<Scholz2015GeometryPath>();
    soleus_l_path.setName("soleus_l_path");
    soleus_l_path.setOrigin(*leftShank, SimTK::Vec3(-0.0024, 0.0334, -0.0071));
    soleus_l_path.setInsertion(*leftFoot, SimTK::Vec3(-0.095, 0.001, 0.0053));

    auto* tib_ant_l = createMuscle<MuscleType>(model, "tib_ant_l", 1579,
                                              0.0683, 0.243, 0.1676);
    tib_ant_l->set_path(Scholz2015GeometryPath());
    auto& tib_ant_l_path = tib_ant_l->template updPath<Scholz2015GeometryPath>();
    tib_ant_l_path.setName("tib_ant_l_path");
    tib_ant_l_path.setOrigin(*leftShank, SimTK::Vec3(0.0179, 0.0243, -0.0115));
    tib_ant_l_path.setInsertion(*leftFoot, SimTK::Vec3(0.0166, -0.0122, -0.0305));
    tib_ant_l_path.addViaPoint(*leftShank, SimTK::Vec3(0.0329, -0.2084, -0.0177));
}

void addJointStiffnessAndDamping(Model& model, const std::string& coordName,
            double stiffness, double damping) {
    CoordinateLinearDamper* damper =
        new CoordinateLinearDamper(coordName, damping);
    damper->setName(coordName + "_damper");
    model.addForce(damper);

    CoordinateLinearSpring* spring =
        new CoordinateLinearSpring(coordName, stiffness);
    spring->setName(coordName + "_spring");
    model.addForce(spring);
}

void addJointStop(Model& model, const std::string& coordName, double stiffness,
            double damping, double qLow, double qHigh) {
    CoordinateLinearStop* stop =
        new CoordinateLinearStop(coordName, stiffness, damping,
                            SimTK::convertDegreesToRadians(qLow),
                            SimTK::convertDegreesToRadians(qHigh));
    stop->setName(coordName + "_stop");
    model.addForce(stop);
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool visualize = false;
    bool show = false;
    double controlValue = 0.1;  // Default control value
    bool randomizeSpeeds = true;  // Default to randomizing speeds
    JointType jointType = JointType::Custom;
    bool useObstacles = false;  // Default to using obstacles

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--visualize" || arg == "-v") {
            visualize = true;
        } else if (arg == "--show" || arg == "-s") {
            show = true;
        } else if (arg == "--control-value" || arg == "-c") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --control-value requires a value\n";
                return EXIT_FAILURE;
            }
            try {
                controlValue = std::stod(argv[++i]);
                if (controlValue < 0.0 || controlValue > 1.0) {
                    std::cerr << "Error: Control value must be between 0.0 and 1.0\n";
                    return EXIT_FAILURE;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid control value '" << argv[i] << "'\n";
                return EXIT_FAILURE;
            }
        } else if (arg == "--joint-type" || arg == "-j") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --joint-type requires a value\n";
                return EXIT_FAILURE;
            }
            std::string type = argv[++i];
            if (type == "custom") {
                jointType = JointType::Custom;
            } else if (type == "ball") {
                jointType = JointType::Ball;
            } else {
                std::cerr << "Error: Unknown hip joint type '" << type << "'\n";
                return EXIT_FAILURE;
            }
        } else if (arg == "--use-obstacles" || arg == "-o") {
            useObstacles = true;
        } else if (arg == "--randomize-speeds" || arg == "-r") {
            randomizeSpeeds = true;
        } else if (arg == "--no-randomize-speeds" || arg == "-n") {
            randomizeSpeeds = false;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
            std::cout << "Options:\n";
            std::cout << "  --visualize, -v              Enable visualization\n";
            std::cout << "  --show, -s                   Show the model in the visualizer\n";
            std::cout << "  --joint-type TYPE, -j TYPE    Set joint type (custom or ball)\n";
            std::cout << "  --use-obstacles, -o          Use obstacles in the simulation\n";
            std::cout << "  --control-value VALUE, -c VALUE  Set control value\n";
            std::cout << "                                   (0.0 to 1.0, default: 0.1)\n";
            std::cout << "  --randomize-speeds, -r        Randomize initial speeds (default)\n";
            std::cout << "  --no-randomize-speeds, -n     Use zero initial speeds\n";
            std::cout << "  --help, -h                   Show this help message\n";
            std::cout << "\nControl Value:\n";
            std::cout << "  Sets the constant excitation value\n";
            std::cout << "\nSpeed Initialization:\n";
            std::cout << "  --randomize-speeds: Start with random initial velocities (default)\n";
            std::cout << "  --no-randomize-speeds: Start with zero initial velocities\n";
            return EXIT_SUCCESS;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            std::cerr << "Use --help for usage information.\n";
            return EXIT_FAILURE;
        }
    }

    // Display control value
    std::cout << "Using control value: " << controlValue << "\n";

    // Display speed initialization setting
    if (randomizeSpeeds) {
        std::cout << "Using randomized initial speeds\n";
    } else {
        std::cout << "Using zero initial speeds\n";
    }


    // Create the model.
    // ----------------
    Model model;

    // Bodies
    // ------
    Body* pelvis = new Body("pelvis", massData[Pelvis], SimTK::Vec3(0),
                           SimTK::Inertia(inertiaData[Pelvis]));
    PhysicalOffsetFrame* pelvisOffset = new PhysicalOffsetFrame("pelvis_offset",
        *pelvis, SimTK::Vec3(-0.01, -0.05, 0));
    pelvisOffset->attachGeometry(new Ellipsoid(0.07, 0.07, 0.12));
    pelvis->addComponent(pelvisOffset);

    Body* torso = new Body("torso", massData[Torso], SimTK::Vec3(0),
                          SimTK::Inertia(inertiaData[Torso]));
    torso->attachGeometry(new Ellipsoid(0.1, 0.27, 0.1));
    addContactGeometry(torso, SimTK::Vec3(0, 0.38, 0), "torso_offset", 0.09);

    Body* leftThigh = new Body("leftThigh", massData[LeftThigh], SimTK::Vec3(0),
                              SimTK::Inertia(inertiaData[LeftThigh]));
    leftThigh->attachGeometry(new Ellipsoid(0.04, 0.2, 0.04));

    Body* leftShank = new Body("leftShank", massData[LeftShank], SimTK::Vec3(0),
                              SimTK::Inertia(inertiaData[LeftShank]));
    leftShank->attachGeometry(new Cylinder(0.02, 0.22));

    Body* leftFoot = new Body("leftFoot", massData[LeftFoot], SimTK::Vec3(0),
                             SimTK::Inertia(inertiaData[LeftFoot]));
    leftFoot->attachGeometry(new Ellipsoid(0.1, 0.03, 0.05));
    addContactGeometry(leftFoot, leftContactPoints[0], "heel_geometry", 0.02);
    addContactGeometry(leftFoot, leftContactPoints[1], "lateralToe_geometry", 0.02);
    addContactGeometry(leftFoot, leftContactPoints[2], "medialToe_geometry", 0.02);

    Body* rightThigh = new Body("rightThigh", massData[RightThigh], SimTK::Vec3(0),
                               SimTK::Inertia(inertiaData[RightThigh]));
    rightThigh->attachGeometry(new Ellipsoid(0.04, 0.2, 0.04));

    Body* rightShank = new Body("rightShank", massData[RightShank], SimTK::Vec3(0),
                               SimTK::Inertia(inertiaData[RightShank]));
    rightShank->attachGeometry(new Cylinder(0.02, 0.22));

    Body* rightFoot = new Body("rightFoot", massData[RightFoot], SimTK::Vec3(0),
                              SimTK::Inertia(inertiaData[RightFoot]));
    rightFoot->attachGeometry(new Ellipsoid(0.1, 0.03, 0.05));
    addContactGeometry(rightFoot, rightContactPoints[0], "heel_geometry", 0.02);
    addContactGeometry(rightFoot, rightContactPoints[1], "lateralToe_geometry", 0.02);
    addContactGeometry(rightFoot, rightContactPoints[2], "medialToe_geometry", 0.02);

    model.addBody(pelvis);
    model.addBody(torso);
    model.addBody(leftThigh);
    model.addBody(leftShank);
    model.addBody(leftFoot);
    model.addBody(rightThigh);
    model.addBody(rightShank);
    model.addBody(rightFoot);

    // Joints
    // ------
    FreeJoint* pelvisGround = new FreeJoint("pelvis_ground",
            model.getGround(), *pelvis);

    Joint* lumbar;
    Joint* leftHip;
    Joint* rightHip;
    if (jointType == JointType::Custom) {
        SpatialTransform lumbar_transform;
        lumbar_transform[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_coord_0", 1, 1));
        lumbar_transform[0].setFunction(new LinearFunction());
        lumbar_transform[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_coord_1", 1, 1));
        lumbar_transform[1].setFunction(new LinearFunction());
        lumbar_transform[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_coord_2", 1, 1));
        lumbar_transform[2].setFunction(new LinearFunction());
        lumbar = new CustomJoint("lumbar",
            *pelvis, SimTK::Vec3(0, 0.05, 0), SimTK::Vec3(0),
            *torso, SimTK::Vec3(0, -0.25, 0), SimTK::Vec3(0),
            lumbar_transform);

        SpatialTransform hip_l_transform;
        hip_l_transform[0].setCoordinateNames(OpenSim::Array<std::string>("hip_l_coord_0", 1, 1));
        hip_l_transform[0].setFunction(new LinearFunction());
        hip_l_transform[1].setCoordinateNames(OpenSim::Array<std::string>("hip_l_coord_1", 1, 1));
        hip_l_transform[1].setFunction(new LinearFunction());
        hip_l_transform[2].setCoordinateNames(OpenSim::Array<std::string>("hip_l_coord_2", 1, 1));
        hip_l_transform[2].setFunction(new LinearFunction());
        leftHip = new CustomJoint("hip_l",
            *pelvis, SimTK::Vec3(0, -0.0661, -0.0835), SimTK::Vec3(0),
            *leftThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0),
            hip_l_transform);

        SpatialTransform hip_r_transform;
        hip_r_transform[0].setCoordinateNames(OpenSim::Array<std::string>("hip_r_coord_0", 1, 1));
        hip_r_transform[0].setFunction(new LinearFunction());
        hip_r_transform[1].setCoordinateNames(OpenSim::Array<std::string>("hip_r_coord_1", 1, 1));
        hip_r_transform[1].setFunction(new LinearFunction());
        hip_r_transform[2].setCoordinateNames(OpenSim::Array<std::string>("hip_r_coord_2", 1, 1));
        hip_r_transform[2].setFunction(new LinearFunction());
        rightHip = new CustomJoint("hip_r",
            *pelvis, SimTK::Vec3(0, -0.0661, 0.0835), SimTK::Vec3(0),
            *rightThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0),
            hip_r_transform);

    } else {
        lumbar = new BallJoint("lumbar",
            *pelvis, SimTK::Vec3(0, 0.05, 0), SimTK::Vec3(0),
            *torso, SimTK::Vec3(0, -0.25, 0), SimTK::Vec3(0));

        leftHip = new BallJoint("hip_l",
            *pelvis, SimTK::Vec3(0, -0.0661, -0.0835), SimTK::Vec3(0),
            *leftThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0));

        rightHip = new BallJoint("hip_r",
            *pelvis, SimTK::Vec3(0, -0.0661, 0.0835), SimTK::Vec3(0),
            *rightThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0));
    }

    PinJoint* leftKnee = new PinJoint("knee_l",
            *leftThigh, SimTK::Vec3(0, -0.226, 0), SimTK::Vec3(0),
            *leftShank, SimTK::Vec3(0, 0.1867, 0), SimTK::Vec3(0));

    PinJoint* leftAnkle = new PinJoint("ankle_l",
            *leftShank, SimTK::Vec3(0, -0.2433, 0), SimTK::Vec3(0),
            *leftFoot, SimTK::Vec3(-0.05123, 0.01195, 0.00792), SimTK::Vec3(0));

    PinJoint* rightKnee = new PinJoint("knee_r",
            *rightThigh, SimTK::Vec3(0, -0.226, 0), SimTK::Vec3(0),
            *rightShank, SimTK::Vec3(0, 0.1867, 0), SimTK::Vec3(0));

    PinJoint* rightAnkle = new PinJoint("ankle_r",
            *rightShank, SimTK::Vec3(0, -0.2433, 0), SimTK::Vec3(0),
            *rightFoot, SimTK::Vec3(-0.05123, 0.01195, -0.00792), SimTK::Vec3(0));

    model.addJoint(pelvisGround);
    model.addJoint(lumbar);
    model.addJoint(leftHip);
    model.addJoint(leftKnee);
    model.addJoint(leftAnkle);
    model.addJoint(rightHip);
    model.addJoint(rightKnee);
    model.addJoint(rightAnkle);

    // Muscles
    // ------
    createMuscles<Millard2012EquilibriumMuscle>(model, pelvis, torso,
                                            leftThigh, leftShank, leftFoot,
                                            rightThigh, rightShank,
                                            rightFoot, useObstacles);

    // Lumbar torque actuators
    // -----------------------
    ActivationCoordinateActuator* lumbarTorqueX =
        new ActivationCoordinateActuator("lumbar_coord_0");
    lumbarTorqueX->setName("lumbar_coord_0_torque");
    lumbarTorqueX->setOptimalForce(100.0);
    lumbarTorqueX->setMinControl(-1.0);
    lumbarTorqueX->setMaxControl(1.0);
    lumbarTorqueX->set_activation_time_constant(0.05);
    model.addForce(lumbarTorqueX);

    ActivationCoordinateActuator* lumbarTorqueY =
        new ActivationCoordinateActuator("lumbar_coord_1");
    lumbarTorqueY->setName("lumbar_coord_1_torque");
    lumbarTorqueY->setOptimalForce(100.0);
    lumbarTorqueY->setMinControl(-1.0);
    lumbarTorqueY->setMaxControl(1.0);
    lumbarTorqueY->set_activation_time_constant(0.05);
    model.addForce(lumbarTorqueY);

    ActivationCoordinateActuator* lumbarTorqueZ =
        new ActivationCoordinateActuator("lumbar_coord_2");
    lumbarTorqueZ->setName("lumbar_coord_2_torque");
    lumbarTorqueZ->setOptimalForce(100.0);
    lumbarTorqueZ->setMinControl(-1.0);
    lumbarTorqueZ->setMaxControl(1.0);
    lumbarTorqueZ->set_activation_time_constant(0.05);
    model.addForce(lumbarTorqueZ);

    // Joint damping
    // -------------
    SimTK::Real stiffness = 5.0;
    SimTK::Real damping = 0.5;
    addJointStiffnessAndDamping(model, "lumbar_coord_0", stiffness, 10.0*damping);
    addJointStiffnessAndDamping(model, "lumbar_coord_1", stiffness, 10.0*damping);
    addJointStiffnessAndDamping(model, "lumbar_coord_2", stiffness, 10.0*damping);
    addJointStiffnessAndDamping(model, "hip_l_coord_0", stiffness, damping);
    addJointStiffnessAndDamping(model, "hip_l_coord_1", stiffness, damping);
    addJointStiffnessAndDamping(model, "hip_l_coord_2", stiffness, damping);
    addJointStiffnessAndDamping(model, "hip_r_coord_0", stiffness, damping);
    addJointStiffnessAndDamping(model, "hip_r_coord_1", stiffness, damping);
    addJointStiffnessAndDamping(model, "hip_r_coord_2", stiffness, damping);
    addJointStiffnessAndDamping(model, "knee_l_coord_0", stiffness, damping);
    addJointStiffnessAndDamping(model, "ankle_l_coord_0", stiffness, damping);
    addJointStiffnessAndDamping(model, "knee_r_coord_0", stiffness, damping);
    addJointStiffnessAndDamping(model, "ankle_r_coord_0", stiffness, damping);

    // Joint stops
    // ----------
    SimTK::Real stopStiffness = 100;
    SimTK::Real stopDamping = 5.0;
    if (jointType == JointType::Custom) {
        addJointStop(model, "lumbar_coord_0", stopStiffness, stopDamping,
                     -20.0, 20.0);
        addJointStop(model, "lumbar_coord_1", stopStiffness, stopDamping,
                     -20.0, 20.0);
        addJointStop(model, "lumbar_coord_2", stopStiffness, stopDamping,
                     -20.0, 20.0);
        addJointStop(model, "hip_l_coord_0", stopStiffness, stopDamping,
                     -30.0, 30.0);
        addJointStop(model, "hip_l_coord_1", stopStiffness, stopDamping,
                     -30.0, 30.0);
        addJointStop(model, "hip_l_coord_2", stopStiffness, stopDamping,
                     -30.0, 30.0);
        addJointStop(model, "hip_r_coord_0", stopStiffness, stopDamping,
                     -30.0, 30.0);
        addJointStop(model, "hip_r_coord_1", stopStiffness, stopDamping,
                     -30.0, 30.0);
        addJointStop(model, "hip_r_coord_2", stopStiffness, stopDamping,
                     -30.0, 30.0);
    }

    addJointStop(model, "knee_l_coord_0", stopStiffness, stopDamping,
                 -120.0, -3.0);
    addJointStop(model, "knee_r_coord_0", stopStiffness, stopDamping,
                 -120.0, -3.0);
    addJointStop(model, "ankle_l_coord_0", stopStiffness, stopDamping,
                 -60.0, 25.0);
    addJointStop(model, "ankle_r_coord_0", stopStiffness, stopDamping,
                 -60.0, 25.0);

    // Contact
    // -------
    SimTK::Transform transform(SimTK::Rotation(-0.5*SimTK::Pi, SimTK::XAxis),
                              SimTK::Vec3(0));

    SimTK::ExponentialSpringParameters params;
    params.setNormalViscosity(1.0);
    params.setInitialMuStatic(0.9);
    params.setInitialMuKinetic(0.6);
    params.setSettleVelocity(0.1);

    addContact(model, "left_heel_contact", leftFoot,
               leftContactPoints[0], transform, params);
    addContact(model, "left_lateralToe_contact", leftFoot,
               leftContactPoints[1], transform, params);
    addContact(model, "left_medialToe_contact", leftFoot,
               leftContactPoints[2], transform, params);

    addContact(model, "right_heel_contact", rightFoot,
               rightContactPoints[0], transform, params);
    addContact(model, "right_lateralToe_contact", rightFoot,
               rightContactPoints[1], transform, params);
    addContact(model, "right_medialToe_contact", rightFoot,
               rightContactPoints[2], transform, params);

    // Torso contact points (around the body)
    addContact(model, "torso_front_contact", torso,
               SimTK::Vec3(0, 0.1, 0), transform, params);
    addContact(model, "torso_back_contacpt", torso,
               SimTK::Vec3(0, -0.1, 0), transform, params);
    addContact(model, "torso_left_contact", torso,
               SimTK::Vec3(0, 0, 0.05), transform, params);
    addContact(model, "torso_right_contact", torso,
               SimTK::Vec3(0, 0, -0.05), transform, params);

    // Pelvis contact points (around the pelvis)
    addContact(model, "pelvis_front_contact", pelvis,
               SimTK::Vec3(0, 0.05, 0), transform, params);
    addContact(model, "pelvis_back_contact", pelvis,
               SimTK::Vec3(0, -0.05, 0), transform, params);
    addContact(model, "pelvis_left_contact", pelvis,
               SimTK::Vec3(0, 0, 0.06), transform, params);
    addContact(model, "pelvis_right_contact", pelvis,
               SimTK::Vec3(0, 0, -0.06), transform, params);

    // Left hip contact points
    addContact(model, "left_hip_front_contact", leftThigh,
               SimTK::Vec3(0, 0.15, 0), transform, params);
    addContact(model, "left_hip_back_contact", leftThigh,
               SimTK::Vec3(0, 0.19, 0), transform, params);
    addContact(model, "left_hip_medial_contact", leftThigh,
               SimTK::Vec3(0, 0.17, 0.02), transform, params);
    addContact(model, "left_hip_lateral_contact", leftThigh,
               SimTK::Vec3(0, 0.17, -0.02), transform, params);

    // Right hip contact points
    addContact(model, "right_hip_front_contact", rightThigh,
               SimTK::Vec3(0, 0.15, 0), transform, params);
    addContact(model, "right_hip_back_contact", rightThigh,
               SimTK::Vec3(0, 0.19, 0), transform, params);
    addContact(model, "right_hip_medial_contact", rightThigh,
               SimTK::Vec3(0, 0.17, -0.02), transform, params);
    addContact(model, "right_hip_lateral_contact", rightThigh,
               SimTK::Vec3(0, 0.17, 0.02), transform, params);

    // Left knee contact points
    addContact(model, "left_knee_front_contact", leftShank,
               SimTK::Vec3(0, 0.15, 0), transform, params);
    addContact(model, "left_knee_back_contact", leftShank,
               SimTK::Vec3(0, 0.22, 0), transform, params);
    addContact(model, "left_knee_medial_contact", leftShank,
               SimTK::Vec3(0, 0.1867, 0.015), transform, params);
    addContact(model, "left_knee_lateral_contact", leftShank,
               SimTK::Vec3(0, 0.1867, -0.015), transform, params);

    // Right knee contact points
    addContact(model, "right_knee_front_contact", rightShank,
               SimTK::Vec3(0, 0.15, 0), transform, params);
    addContact(model, "right_knee_back_contact", rightShank,
               SimTK::Vec3(0, 0.22, 0), transform, params);
    addContact(model, "right_knee_medial_contact", rightShank,
               SimTK::Vec3(0, 0.1867, -0.015), transform, params);
    addContact(model, "right_knee_lateral_contact", rightShank,
               SimTK::Vec3(0, 0.1867, 0.015), transform, params);

    // Left ankle contact points
    addContact(model, "left_ankle_front_contact", leftFoot,
               SimTK::Vec3(-0.03, 0.012, 0.008), transform, params);
    addContact(model, "left_ankle_back_contact", leftFoot,
               SimTK::Vec3(-0.07, 0.012, 0.008), transform, params);
    addContact(model, "left_ankle_medial_contact", leftFoot,
               SimTK::Vec3(-0.05123, 0.012, 0.013), transform, params);
    addContact(model, "left_ankle_lateral_contact", leftFoot,
               SimTK::Vec3(-0.05123, 0.012, 0.003), transform, params);

    // Right ankle contact points
    addContact(model, "right_ankle_front_contact", rightFoot,
               SimTK::Vec3(-0.03, 0.012, -0.008), transform, params);
    addContact(model, "right_ankle_back_contact", rightFoot,
               SimTK::Vec3(-0.07, 0.012, -0.008), transform, params);
    addContact(model, "right_ankle_medial_contact", rightFoot,
               SimTK::Vec3(-0.05123, 0.012, -0.013), transform, params);
    addContact(model, "right_ankle_lateral_contact", rightFoot,
               SimTK::Vec3(-0.05123, 0.012, -0.003), transform, params);

    // Force aggregators
    // -----------------
    // model.initSystem();
    // ForceAggregator* limitTorqueAggregator = new ForceAggregator();
    // limitTorqueAggregator->setName("limit_torque_aggregator");
    // for (const auto& limit_torque : model.getComponentList<CoordinateLinearStop>()) {
    //     limitTorqueAggregator->addForce(limit_torque);
    // }
    // model.addComponent(limitTorqueAggregator);

    // ForceAggregator* contactForceAggregator = new ForceAggregator();
    // contactForceAggregator->setName("contact_force_aggregator");
    // contactForceAggregator->addForce(model.getComponent<ExponentialContactForce>("/forceset/left_heel_contact"));
    // contactForceAggregator->addForce(model.getComponent<ExponentialContactForce>("/forceset/left_lateralToe_contact"));
    // contactForceAggregator->addForce(model.getComponent<ExponentialContactForce>("/forceset/left_medialToe_contact"));
    // contactForceAggregator->addForce(model.getComponent<ExponentialContactForce>("/forceset/right_heel_contact"));
    // contactForceAggregator->addForce(model.getComponent<ExponentialContactForce>("/forceset/right_lateralToe_contact"));
    // contactForceAggregator->addForce(model.getComponent<ExponentialContactForce>("/forceset/right_medialToe_contact"));
    // model.addComponent(contactForceAggregator);

    // Controller
    // ---------
    DiscreteController* controller = new DiscreteController();
    controller->setName("controller");
    model.addController(controller);

    // Construct system
    // ---------------
    model.finalizeConnections();
    SimTK::State state = model.initSystem();

    // Default state
    // -------------
    model.getComponent<Coordinate>(
        "/jointset/pelvis_ground/pelvis_ground_coord_4").setValue(state, 1.05);
    model.getComponent<Coordinate>(
        "/jointset/lumbar/lumbar_coord_2").setValue(state, -SimTK::Pi/8);
    model.getComponent<Coordinate>(
        "/jointset/hip_l/hip_l_coord_0").setValue(state,
            SimTK::convertDegreesToRadians(10.0));
    model.getComponent<Coordinate>(
        "/jointset/hip_l/hip_l_coord_2").setValue(state,
            SimTK::convertDegreesToRadians(30.0));
    model.getComponent<Coordinate>(
        "/jointset/hip_r/hip_r_coord_0").setValue(state,
            SimTK::convertDegreesToRadians(-10.0));
    model.getComponent<Coordinate>(
        "/jointset/hip_r/hip_r_coord_2").setValue(state,
            SimTK::convertDegreesToRadians(30.0));
    model.getComponent<Coordinate>(
        "/jointset/knee_l/knee_l_coord_0").setValue(state,
            SimTK::convertDegreesToRadians(-60.0));
    model.getComponent<Coordinate>(
        "/jointset/ankle_l/ankle_l_coord_0").setValue(state,
            SimTK::convertDegreesToRadians(20.0));
    model.getComponent<Coordinate>(
        "/jointset/knee_r/knee_r_coord_0").setValue(state,
            SimTK::convertDegreesToRadians(-60.0));
    model.getComponent<Coordinate>(
        "/jointset/ankle_r/ankle_r_coord_0").setValue(state,
            SimTK::convertDegreesToRadians(20.0));

    // Initialize speeds based on command line option
    if (randomizeSpeeds) {
        // Randomize the initial speeds.
        state.updU() = SimTK::Test::randVector(state.getNU());
    } else {
        // Use zero initial speeds.
        state.updU() = SimTK::Vector(state.getNU(), 0.0);
    }

    // Default controls
    // ----------------
    SimTK::Vector controls(model.getNumControls(), controlValue);
    model.getComponent<DiscreteController>("/controllerset/controller").
        setDiscreteControls(state, controls);

    // Show Model
    // ----------
    if (show) {
        VisualizerUtilities::showModel(model);
    }

    // Simulate
    // --------
    SimTK::Real finalTime = 20.0;
    Manager manager(model);
    manager.setWriteToStorage(visualize);
    manager.setPerformAnalyses(visualize);
    // manager.setRecordStatesTrajectory(true);
    manager.setIntegratorMethod(Manager::IntegratorMethod::CPodes);
    manager.setIntegratorAccuracy(1e-2);
    // manager.setIntegratorMaximumStepSize(1e-2);
    manager.initialize(state);
    manager.integrate(finalTime);

    // StatesTrajectory statesTraj = manager.getStatesTrajectory();
    // const ForceAggregator& limitTorqueAgg =
    //         model.getComponent<ForceAggregator>("/limit_torque_aggregator");
    // const ForceAggregator& contactForceAgg =
    //         model.getComponent<ForceAggregator>("/contact_force_aggregator");

    // for (const auto& s : statesTraj) {
    //     model.realizeAcceleration(s);
    //     const SimTK::Vector& generalizedForces = limitTorqueAgg.getGeneralizedForces(s);
    //     const SimTK::Vector_<SimTK::SpatialVec>& bodyForces = limitTorqueAgg.getBodyForces(s);
    //     std::cout << "limit Torque Aggregator" << std::endl;
    //     std::cout << "Generalized Forces: " << generalizedForces << std::endl;
    //     std::cout << "Body Forces: " << bodyForces << std::endl;

    //     const SimTK::Vector& contactGeneralizedForces = contactForceAgg.getGeneralizedForces(s);
    //     const SimTK::Vector_<SimTK::SpatialVec>& contactBodyForces = contactForceAgg.getBodyForces(s);
    //     std::cout << "Contact Force Aggregator" << std::endl;
    //     std::cout << "Generalized Forces: " << contactGeneralizedForces << std::endl;
    //     std::cout << "Body Forces: " << contactBodyForces << std::endl;
    // }




    // Visualize
    // ---------
    if (visualize) {
        TimeSeriesTable table = manager.getStatesTable();

        std::string statesFileName = fmt::format("states_{}_control{}.sto",
                                                 randomizeSpeeds ? "random" : "zero",
                                                 controlValue);
        STOFileAdapter::write(table, statesFileName);

        VisualizerUtilities::showMotion(model, table);
    }


    return EXIT_SUCCESS;
}
