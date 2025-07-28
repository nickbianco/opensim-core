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
#include <OpenSim/Actuators/HyfydyMuscle.h>
#include <OpenSim/Actuators/DeGrooteFregly2016Muscle.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>
#include <OpenSim/Simulation/Model/ExponentialContactForce.h>
#include <OpenSim/Common/STOFileAdapter.h>

// Include custom components
#include "PointPathMuscle.h"
#include "MobilityLinearDamper.h"
#include "MobilityLinearStop.h"
#include "PointBasedPath.h"
#include "DiscreteController.h"

using namespace OpenSim;

// Enums
enum BodyType {LeftFoot=0, RightFoot, LeftShank, RightShank, LeftThigh, 
               RightThigh, Pelvis, Torso};

enum MuscleType {GlutMed_R=0, AddMag_R, Hamstrings_R, Bifemsh_R, GlutMax_R, 
                 Iliopsoas_R, RectFem_R, Vasti_R, Gastroc_R, Soleus_R, TibAnt_R, 
                 GlutMed_L, AddMag_L, Hamstrings_L, Bifemsh_L, GlutMax_L, 
                 Iliopsoas_L, RectFem_L, Vasti_L, Gastroc_L, Soleus_L, TibAnt_L};

enum Contact {Heel=0, LateralToe, MedialToe};

enum MuscleModelType {PointPath=0, Hyfydy, DeGrooteFregly};

enum PathType {Scholz2015=0, PointBased, Geometry};

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

// Specialization for HyfydyMuscle
template<>
HyfydyMuscle* createMuscle<HyfydyMuscle>(Model& model, const std::string& name, 
                                        double maxIsometricForce, 
                                        double optimalFiberLength, 
                                        double tendonSlackLength, 
                                        double pennationAngle) {
    auto* muscle = new HyfydyMuscle();
    muscle->setName(name);
    muscle->set_max_isometric_force(maxIsometricForce);
    muscle->set_optimal_fiber_length(optimalFiberLength);
    muscle->set_tendon_slack_length(tendonSlackLength);
    muscle->set_pennation_angle_at_optimal(pennationAngle);
    muscle->set_ignore_tendon_compliance(true);
    muscle->set_ignore_activation_dynamics(true);
    muscle->set_default_activation(0.0);
    model.addForce(muscle);
    return muscle;
}

// Template-based path addition function
template<typename MuscleType>
void addPathToMuscle(MuscleType* muscle, const std::string& name, 
                    const std::vector<std::pair<PhysicalFrame*, SimTK::Vec3>>& 
                    points, PathType pathType = Scholz2015) {
    if (pathType == PointBased) {
        // Use PointBasedPath
        auto* path = new PointBasedPath();
        path->setName(name);
        
        for (size_t i = 0; i < points.size(); ++i) {
            path->addStation(name + "_point_" + std::to_string(i), 
                          *points[i].first, points[i].second);
        }
        muscle->setPath(path);
    } else if (pathType == PathType::Geometry) {
        auto* path = new GeometryPath();
        path->setName(name);
        for (size_t i = 0; i < points.size(); ++i) {
            path->appendNewPathPoint(name + "_point_" + std::to_string(i), 
                          *points[i].first, points[i].second);
        }
        muscle->setPath(path);
    } else {
        // Use Scholz2015GeometryPath (default)
        auto* path = new Scholz2015GeometryPath();
        path->setName(name);

        std::vector<Station*> stations;
        for (size_t i = 0; i < points.size(); ++i) {
            auto* station = new Station(*points[i].first, points[i].second);
            station->setName(name + "_point_" + std::to_string(i));
            path->addComponent(station);
            stations.push_back(station);
        }

        int numSegments = static_cast<int>(stations.size()) - 1;
        for (int i = 0; i < numSegments; ++i) {
            if (i == 0) {
                path->createInitialPathSegment("segment_" + std::to_string(i), 
                                             *stations[i], *stations[i + 1]);
            } else {
                path->appendPathSegment("segment_" + std::to_string(i), 
                                      *stations[i + 1]);
            }
        }
        muscle->setPath(path);
    }
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
    auto* station = new Station(*frame, location);
    station->setName(name);
    frame->addComponent(station);
    auto* contact = new ExponentialContactForce(transform, *station, params);
    contact->setName(name);
    model.addForce(contact);
}

void createPointPathMuscles(Model& model, Body* pelvis, Body* torso, 
        Body* leftThigh, Body* leftShank, Body* leftFoot,
        Body* rightThigh, Body* rightShank, Body* rightFoot, 
        double controlValue, bool ignoreActivationDynamics) {
    PointPathMuscle* muscle;

    // Right leg muscles
    // glut_med_r
    muscle = new PointPathMuscle("glut_med_r", 2045, 0.0733, 0.066, 0.3578);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.0148, 0.0445, 0.0766));
    muscle->addPoint("point2", *rightThigh, 
                     SimTK::Vec3(-0.0258, 0.1642, 0.0527));
    model.addForce(muscle);

    // add_mag_r
    muscle = new PointPathMuscle("add_mag_r", 2268, 0.087, 0.06, 0.0872665);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.0025, -0.1174, 0.0255));
    muscle->addPoint("point2", *rightThigh, 
                     SimTK::Vec3(-0.0045, 0.0489, 0.0339));
    model.addForce(muscle);

    // hamstrings_r
    muscle = new PointPathMuscle("hamstrings_r", 2594, 0.0976, 0.319, 0.2025);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.05526, -0.10257, 0.06944));
    muscle->addPoint("point2", *rightShank, 
                     SimTK::Vec3(-0.028, 0.1667, 0.02943));
    muscle->addPoint("point3", *rightShank, 
                     SimTK::Vec3(-0.021, 0.1467, 0.0343));
    model.addForce(muscle);

    // bifemsh_r
    muscle = new PointPathMuscle("bifemsh_r", 804, 0.1103, 0.095, 0.2147);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *rightThigh, 
                     SimTK::Vec3(0.005, -0.0411, 0.0234));
    muscle->addPoint("point2", *rightShank, 
                     SimTK::Vec3(-0.028, 0.1667, 0.02943));
    muscle->addPoint("point3", *rightShank, 
                     SimTK::Vec3(-0.021, 0.1467, 0.0343));
    model.addForce(muscle);

    // glut_max_r
    muscle = new PointPathMuscle("glut_max_r", 1944, 0.1569, 0.111, 0.3822);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.0642, 0.0176, 0.0563));
    muscle->addPoint("point2", *pelvis, 
                     SimTK::Vec3(-0.0669, -0.052, 0.0914));
    muscle->addPoint("point3", *rightThigh, 
                     SimTK::Vec3(-0.0426, 0.117, 0.0293));
    muscle->addPoint("point4", *rightThigh, 
                     SimTK::Vec3(-0.0156, 0.0684, 0.0419));
    model.addForce(muscle);

    // iliopsoas_r
    muscle = new PointPathMuscle("iliopsoas_r", 2186, 0.1066, 0.152, 0.2496);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(0.006, 0.0887, 0.0289));
    muscle->addPoint("point2", *pelvis, 
                     SimTK::Vec3(0.0407, -0.01, 0.076));
    muscle->addPoint("point3", *rightThigh, 
                     SimTK::Vec3(0.033, 0.135, 0.0038));
    muscle->addPoint("point4", *rightThigh, 
                     SimTK::Vec3(-0.0188, 0.1103, 0.0104));
    model.addForce(muscle);

    // rect_fem_r
    muscle = new PointPathMuscle("rect_fem_r", 1169, 0.0759, 0.3449, 0.2426);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(0.0412, -0.0311, 0.0968));
    muscle->addPoint("point2", *rightThigh, 
                     SimTK::Vec3(0.038, -0.17, 0.004));
    muscle->addPoint("point3", *rightShank, 
                     SimTK::Vec3(0.038, 0.2117, 0.0018));
    model.addForce(muscle);

    // vasti_r
    muscle = new PointPathMuscle("vasti_r", 4530, 0.0993, 0.1231, 0.0785);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *rightThigh, 
                     SimTK::Vec3(0.029, -0.0224, 0.031));
    muscle->addPoint("point2", *rightThigh, 
                     SimTK::Vec3(0.038, -0.17, 0.007));
    muscle->addPoint("point3", *rightShank, 
                     SimTK::Vec3(0.038, 0.2117, 0.0018));
    model.addForce(muscle);

    // gastroc_r
    muscle = new PointPathMuscle("gastroc_r", 2241, 0.051, 0.384, 0.1728);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *rightThigh, 
                     SimTK::Vec3(-0.02, -0.218, -0.024));
    muscle->addPoint("point2", *rightFoot, 
                     SimTK::Vec3(-0.095, 0.001, -0.0053));
    model.addForce(muscle);

    // soleus_r
    muscle = new PointPathMuscle("soleus_r", 3549, 0.044, 0.248, 0.4939);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *rightShank, 
                     SimTK::Vec3(-0.0024, 0.0334, 0.0071));
    muscle->addPoint("point2", *rightFoot, 
                     SimTK::Vec3(-0.095, 0.001, -0.0053));
    model.addForce(muscle);

    // tib_ant_r
    muscle = new PointPathMuscle("tib_ant_r", 1579, 0.0683, 0.243, 0.1676);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *rightShank, 
                     SimTK::Vec3(0.0179, 0.0243, 0.0115));
    muscle->addPoint("point2", *rightShank, 
                     SimTK::Vec3(0.0329, -0.2084, -0.0177));
    muscle->addPoint("point3", *rightFoot, 
                     SimTK::Vec3(0.0166, -0.0122, -0.0305));
    model.addForce(muscle);

    // Left leg muscles
    // glut_med_l
    muscle = new PointPathMuscle("glut_med_l", 2045, 0.0733, 0.066, 0.3578);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.0148, 0.0445, -0.0766));
    muscle->addPoint("point2", *leftThigh, 
                     SimTK::Vec3(-0.0258, 0.1642, -0.0527));
    model.addForce(muscle);

    // add_mag_l
    muscle = new PointPathMuscle("add_mag_l", 2268, 0.087, 0.06, 0.0872665);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.0025, -0.1174, -0.0255));
    muscle->addPoint("point2", *leftThigh, 
                     SimTK::Vec3(-0.0045, 0.0489, -0.0339));
    model.addForce(muscle);

    // hamstrings_l
    muscle = new PointPathMuscle("hamstrings_l", 2594, 0.0976, 0.319, 0.2025);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.05526, -0.10257, -0.06944));
    muscle->addPoint("point2", *leftShank, 
                     SimTK::Vec3(-0.028, 0.1667, -0.02943));
    muscle->addPoint("point3", *leftShank, 
                     SimTK::Vec3(-0.021, 0.1467, -0.0343));
    model.addForce(muscle);

    // bifemsh_l
    muscle = new PointPathMuscle("bifemsh_l", 804, 0.1103, 0.095, 0.2147);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *leftThigh, 
                     SimTK::Vec3(0.005, -0.0411, -0.0234));
    muscle->addPoint("point2", *leftShank, 
                     SimTK::Vec3(-0.028, 0.1667, -0.02943));
    muscle->addPoint("point3", *leftShank, 
                     SimTK::Vec3(-0.021, 0.1467, -0.0343));
    model.addForce(muscle);

    // glut_max_l
    muscle = new PointPathMuscle("glut_max_l", 1944, 0.1569, 0.111, 0.3822);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(-0.0642, 0.0176, -0.0563));
    muscle->addPoint("point2", *pelvis, 
                     SimTK::Vec3(-0.0669, -0.052, -0.0914));
    muscle->addPoint("point3", *leftThigh, 
                     SimTK::Vec3(-0.0426, 0.117, -0.0293));
    muscle->addPoint("point4", *leftThigh, 
                     SimTK::Vec3(-0.0156, 0.0684, -0.0419));
    model.addForce(muscle);

    // iliopsoas_l
    muscle = new PointPathMuscle("iliopsoas_l", 2186, 0.1066, 0.152, 0.2496);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(0.006, 0.0887, -0.0289));
    muscle->addPoint("point2", *pelvis, 
                     SimTK::Vec3(0.0407, -0.01, -0.076));
    muscle->addPoint("point3", *leftThigh, 
                     SimTK::Vec3(0.033, 0.135, -0.0038));
    muscle->addPoint("point4", *leftThigh, 
                     SimTK::Vec3(-0.0188, 0.1103, -0.0104));
    model.addForce(muscle);

    // rect_fem_l
    muscle = new PointPathMuscle("rect_fem_l", 1169, 0.0759, 0.3449, 0.2426);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *pelvis, 
                     SimTK::Vec3(0.0412, -0.0311, -0.0968));
    muscle->addPoint("point2", *leftThigh, 
                     SimTK::Vec3(0.038, -0.17, -0.004));
    muscle->addPoint("point3", *leftShank, 
                     SimTK::Vec3(0.038, 0.2117, -0.0018));
    model.addForce(muscle);

    // vasti_l
    muscle = new PointPathMuscle("vasti_l", 4530, 0.0993, 0.1231, 0.0785);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *leftThigh, 
                     SimTK::Vec3(0.029, -0.0224, -0.031));
    muscle->addPoint("point2", *leftThigh, 
                     SimTK::Vec3(0.038, -0.17, -0.007));
    muscle->addPoint("point3", *leftShank, 
                     SimTK::Vec3(0.038, 0.2117, -0.0018));
    model.addForce(muscle);

    // gastroc_l
    muscle = new PointPathMuscle("gastroc_l", 2241, 0.051, 0.384, 0.1728);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *leftThigh, 
                     SimTK::Vec3(-0.02, -0.218, 0.024));
    muscle->addPoint("point2", *leftFoot, 
                     SimTK::Vec3(-0.095, 0.001, 0.0053));
    model.addForce(muscle);

    // soleus_l
    muscle = new PointPathMuscle("soleus_l", 3549, 0.044, 0.248, 0.4939);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *leftShank, 
                     SimTK::Vec3(-0.0024, 0.0334, -0.0071));
    muscle->addPoint("point2", *leftFoot, 
                     SimTK::Vec3(-0.095, 0.001, 0.0053));
    model.addForce(muscle);

    // tib_ant_l
    muscle = new PointPathMuscle("tib_ant_l", 1579, 0.0683, 0.243, 0.1676);
    muscle->set_ignore_activation_dynamics(ignoreActivationDynamics);
    muscle->setControlValue(controlValue);
    muscle->addPoint("point1", *leftShank, 
                     SimTK::Vec3(0.0179, 0.0243, -0.0115));
    muscle->addPoint("point2", *leftShank, 
                     SimTK::Vec3(0.0329, -0.2084, 0.0177));
    muscle->addPoint("point3", *leftFoot, 
                     SimTK::Vec3(0.0166, -0.0122, 0.0305));
    model.addForce(muscle);
}

template<typename MuscleType>
void createMuscles(Model& model, Body* pelvis, Body* torso, Body* leftThigh, 
                   Body* leftShank, Body* leftFoot, Body* rightThigh, 
                   Body* rightShank, Body* rightFoot, PathType pathType) {
    std::vector<std::pair<PhysicalFrame*, SimTK::Vec3>> points;
    
    // Right leg muscles
    auto* glut_med_r = createMuscle<MuscleType>(model, "glut_med_r", 2045, 
                                               0.0733, 0.066, 0.3578);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0148, 0.0445, 0.0766));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0258, 0.1642, 0.0527));
    addPathToMuscle(glut_med_r, "glut_med_r_path", points, pathType);

    auto* add_mag_r = createMuscle<MuscleType>(model, "add_mag_r", 2268, 
                                              0.087, 0.06, 0.0872665);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0025, -0.1174, 0.0255));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0045, 0.0489, 0.0339));
    addPathToMuscle(add_mag_r, "add_mag_r_path", points, pathType);

    auto* hamstrings_r = createMuscle<MuscleType>(model, "hamstrings_r", 2594, 
                                                 0.0976, 0.319, 0.2025);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.05526, -0.10257, 0.06944));
    points.emplace_back(rightShank, SimTK::Vec3(-0.028, 0.1667, 0.02943));
    points.emplace_back(rightShank, SimTK::Vec3(-0.021, 0.1467, 0.0343));
    addPathToMuscle(hamstrings_r, "hamstrings_r_path", points, pathType);

    auto* bifemsh_r = createMuscle<MuscleType>(model, "bifemsh_r", 804, 
                                              0.1103, 0.095, 0.2147);
    points.clear();
    points.emplace_back(rightThigh, SimTK::Vec3(0.005, -0.0411, 0.0234));
    points.emplace_back(rightShank, SimTK::Vec3(-0.028, 0.1667, 0.02943));
    points.emplace_back(rightShank, SimTK::Vec3(-0.021, 0.1467, 0.0343));
    addPathToMuscle(bifemsh_r, "bifemsh_r_path", points, pathType);

    auto* glut_max_r = createMuscle<MuscleType>(model, "glut_max_r", 1944, 
                                               0.1569, 0.111, 0.3822);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0642, 0.0176, 0.0563));
    points.emplace_back(pelvis, SimTK::Vec3(-0.0669, -0.052, 0.0914));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0426, 0.117, 0.0293));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0156, 0.0684, 0.0419));
    addPathToMuscle(glut_max_r, "glut_max_r_path", points, pathType);

    auto* iliopsoas_r = createMuscle<MuscleType>(model, "iliopsoas_r", 2186, 
                                                0.1066, 0.152, 0.2496);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.006, 0.0887, 0.0289));
    points.emplace_back(pelvis, SimTK::Vec3(0.0407, -0.01, 0.076));
    points.emplace_back(rightThigh, SimTK::Vec3(0.033, 0.135, 0.0038));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0188, 0.1103, 0.0104));
    addPathToMuscle(iliopsoas_r, "iliopsoas_r_path", points, pathType);

    auto* rect_fem_r = createMuscle<MuscleType>(model, "rect_fem_r", 1169, 
                                               0.0759, 0.3449, 0.2426);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.0412, -0.0311, 0.0968));
    points.emplace_back(rightThigh, SimTK::Vec3(0.038, -0.17, 0.004));
    points.emplace_back(rightShank, SimTK::Vec3(0.038, 0.2117, 0.0018));
    addPathToMuscle(rect_fem_r, "rect_fem_r_path", points, pathType);

    auto* vasti_r = createMuscle<MuscleType>(model, "vasti_r", 4530, 
                                            0.0993, 0.1231, 0.0785);
    points.clear();
    points.emplace_back(rightThigh, SimTK::Vec3(0.029, -0.0224, 0.031));
    points.emplace_back(rightThigh, SimTK::Vec3(0.038, -0.17, 0.007));
    points.emplace_back(rightShank, SimTK::Vec3(0.038, 0.2117, 0.0018));
    addPathToMuscle(vasti_r, "vasti_r_path", points, pathType);

    auto* gastroc_r = createMuscle<MuscleType>(model, "gastroc_r", 2241, 
                                              0.051, 0.384, 0.1728);
    points.clear();
    points.emplace_back(rightThigh, SimTK::Vec3(-0.02, -0.218, -0.024));
    points.emplace_back(rightFoot, SimTK::Vec3(-0.095, 0.001, -0.0053));
    addPathToMuscle(gastroc_r, "gastroc_r_path", points, pathType);

    auto* soleus_r = createMuscle<MuscleType>(model, "soleus_r", 3549, 
                                             0.044, 0.248, 0.4939);
    points.clear();
    points.emplace_back(rightShank, SimTK::Vec3(-0.0024, 0.0334, 0.0071));
    points.emplace_back(rightFoot, SimTK::Vec3(-0.095, 0.001, -0.0053));
    addPathToMuscle(soleus_r, "soleus_r_path", points, pathType);

    auto* tib_ant_r = createMuscle<MuscleType>(model, "tib_ant_r", 1579, 
                                              0.0683, 0.243, 0.1676);
    points.clear();
    points.emplace_back(rightShank, SimTK::Vec3(0.0179, 0.0243, 0.0115));
    points.emplace_back(rightShank, SimTK::Vec3(0.0329, -0.2084, -0.0177));
    points.emplace_back(rightFoot, SimTK::Vec3(0.0166, -0.0122, -0.0305));
    addPathToMuscle(tib_ant_r, "tib_ant_r_path", points, pathType);

    // Left leg muscles
    auto* glut_med_l = createMuscle<MuscleType>(model, "glut_med_l", 2045, 
                                               0.0733, 0.066, 0.3578);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0148, 0.0445, -0.0766));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0258, 0.1642, -0.0527));
    addPathToMuscle(glut_med_l, "glut_med_l_path", points, pathType);

    auto* add_mag_l = createMuscle<MuscleType>(model, "add_mag_l", 2268, 
                                              0.087, 0.06, 0.0872665);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0025, -0.1174, -0.0255));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0045, 0.0489, -0.0339));
    addPathToMuscle(add_mag_l, "add_mag_l_path", points, pathType);

    auto* hamstrings_l = createMuscle<MuscleType>(model, "hamstrings_l", 2594, 
                                                 0.0976, 0.319, 0.2025);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.05526, -0.10257, -0.06944));
    points.emplace_back(leftShank, SimTK::Vec3(-0.028, 0.1667, -0.02943));
    points.emplace_back(leftShank, SimTK::Vec3(-0.021, 0.1467, -0.0343));
    addPathToMuscle(hamstrings_l, "hamstrings_l_path", points, pathType);

    auto* bifemsh_l = createMuscle<MuscleType>(model, "bifemsh_l", 804, 
                                              0.1103, 0.095, 0.2147);
    points.clear();
    points.emplace_back(leftThigh, SimTK::Vec3(0.005, -0.0411, -0.0234));
    points.emplace_back(leftShank, SimTK::Vec3(-0.028, 0.1667, -0.02943));
    points.emplace_back(leftShank, SimTK::Vec3(-0.021, 0.1467, -0.0343));
    addPathToMuscle(bifemsh_l, "bifemsh_l_path", points, pathType);

    auto* glut_max_l = createMuscle<MuscleType>(model, "glut_max_l", 1944, 
                                               0.1569, 0.111, 0.3822);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0642, 0.0176, -0.0563));
    points.emplace_back(pelvis, SimTK::Vec3(-0.0669, -0.052, -0.0914));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0426, 0.117, -0.0293));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0156, 0.0684, -0.0419));
    addPathToMuscle(glut_max_l, "glut_max_l_path", points, pathType);

    auto* iliopsoas_l = createMuscle<MuscleType>(model, "iliopsoas_l", 2186, 
                                                0.1066, 0.152, 0.2496);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.006, 0.0887, -0.0289));
    points.emplace_back(pelvis, SimTK::Vec3(0.0407, -0.01, -0.076));
    points.emplace_back(leftThigh, SimTK::Vec3(0.033, 0.135, -0.0038));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0188, 0.1103, -0.0104));
    addPathToMuscle(iliopsoas_l, "iliopsoas_l_path", points, pathType);

    auto* rect_fem_l = createMuscle<MuscleType>(model, "rect_fem_l", 1169, 
                                               0.0759, 0.3449, 0.2426);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.0412, -0.0311, -0.0968));
    points.emplace_back(leftThigh, SimTK::Vec3(0.038, -0.17, -0.004));
    points.emplace_back(leftShank, SimTK::Vec3(0.038, 0.2117, -0.0018));
    addPathToMuscle(rect_fem_l, "rect_fem_l_path", points, pathType);

    auto* vasti_l = createMuscle<MuscleType>(model, "vasti_l", 4530, 
                                            0.0993, 0.1231, 0.0785);
    points.clear();
    points.emplace_back(leftThigh, SimTK::Vec3(0.029, -0.0224, -0.031));
    points.emplace_back(leftThigh, SimTK::Vec3(0.038, -0.17, -0.007));
    points.emplace_back(leftShank, SimTK::Vec3(0.038, 0.2117, -0.0018));
    addPathToMuscle(vasti_l, "vasti_l_path", points, pathType);

    auto* gastroc_l = createMuscle<MuscleType>(model, "gastroc_l", 2241, 
                                              0.051, 0.384, 0.1728);
    points.clear();
    points.emplace_back(leftThigh, SimTK::Vec3(-0.02, -0.218, 0.024));
    points.emplace_back(leftFoot, SimTK::Vec3(-0.095, 0.001, 0.0053));
    addPathToMuscle(gastroc_l, "gastroc_l_path", points, pathType);

    auto* soleus_l = createMuscle<MuscleType>(model, "soleus_l", 3549, 
                                             0.044, 0.248, 0.4939);
    points.clear();
    points.emplace_back(leftShank, SimTK::Vec3(-0.0024, 0.0334, -0.0071));
    points.emplace_back(leftFoot, SimTK::Vec3(-0.095, 0.001, 0.0053));
    addPathToMuscle(soleus_l, "soleus_l_path", points, pathType);

    auto* tib_ant_l = createMuscle<MuscleType>(model, "tib_ant_l", 1579, 
                                              0.0683, 0.243, 0.1676);
    points.clear();
    points.emplace_back(leftShank, SimTK::Vec3(0.0179, 0.0243, -0.0115));
    points.emplace_back(leftShank, SimTK::Vec3(0.0329, -0.2084, 0.0177));
    points.emplace_back(leftFoot, SimTK::Vec3(0.0166, -0.0122, 0.0305));
    addPathToMuscle(tib_ant_l, "tib_ant_l_path", points, pathType);
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool visualize = false;
    MuscleModelType muscleModel = PointPath;
    PathType pathType = Scholz2015;
    double controlValue = 0.1;  // Default control value
    bool randomizeSpeeds = true;  // Default to randomizing speeds
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--visualize" || arg == "-v") {
            visualize = true;
        } else if (arg == "--muscle-model" || arg == "-m") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --muscle-model requires a value\n";
                return EXIT_FAILURE;
            }
            std::string modelType = argv[++i];
            if (modelType == "pointpath") {
                muscleModel = PointPath;
            } else if (modelType == "hyfydy") {
                muscleModel = Hyfydy;
            } else if (modelType == "degroote") {
                muscleModel = DeGrooteFregly;
            } else {
                std::cerr << "Error: Unknown muscle model '" << modelType << "'\n";
                std::cerr << "Valid options: pointpath, hyfydy, degroote\n";
                return EXIT_FAILURE;
            }
        } else if (arg == "--path-type" || arg == "-p") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --path-type requires a value\n";
                return EXIT_FAILURE;
            }
            std::string pathTypeStr = argv[++i];
            if (pathTypeStr == "scholz2015") {
                pathType = Scholz2015;
            } else if (pathTypeStr == "pointbased") {
                pathType = PointBased;
            } else if (pathTypeStr == "geometry") {
                pathType = PathType::Geometry;
            } else {
                std::cerr << "Error: Unknown path type '" << pathTypeStr << "'\n";
                std::cerr << "Valid options: scholz2015, pointbased\n";
                return EXIT_FAILURE;
            }
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
        } else if (arg == "--randomize-speeds" || arg == "-r") {
            randomizeSpeeds = true;
        } else if (arg == "--no-randomize-speeds" || arg == "-n") {
            randomizeSpeeds = false;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
            std::cout << "Options:\n";
            std::cout << "  --visualize, -v              Enable visualization\n";
            std::cout << "  --muscle-model MODEL, -m MODEL  Select muscle model\n";
            std::cout << "                                (pointpath, hyfydy, degroote)\n";
            std::cout << "  --path-type TYPE, -p TYPE     Select path type\n";
            std::cout << "                                (scholz2015, pointbased)\n";
            std::cout << "  --control-value VALUE, -c VALUE  Set control value for PointPathMuscle\n";
            std::cout << "                                (0.0 to 1.0, default: 0.1)\n";
            std::cout << "  --randomize-speeds, -r        Randomize initial speeds (default)\n";
            std::cout << "  --no-randomize-speeds, -n     Use zero initial speeds\n";
            std::cout << "  --help, -h                   Show this help message\n";
            std::cout << "\nMuscle Models:\n";
            std::cout << "  pointpath     Use PointPathMuscle (default)\n";
            std::cout << "  hyfydy        Use HyfydyMuscle\n";
            std::cout << "  degroote      Use DeGrooteFregly2016Muscle\n";
            std::cout << "\nPath Types:\n";
            std::cout << "  scholz2015    Use Scholz2015GeometryPath (default)\n";
            std::cout << "  pointbased    Use PointBasedPath\n";
            std::cout << "\nControl Value:\n";
            std::cout << "  Sets the constant excitation value for PointPathMuscle (0.0 = no activation, 1.0 = full activation)\n";
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
    
    // Select muscle type based on command line option
    bool usePointPathMuscle = (muscleModel == PointPath);
    
    // Set up muscle type based on selection
    if (muscleModel == DeGrooteFregly) {
        std::cout << "Using DeGrooteFregly2016Muscle model\n";
    } else if (muscleModel == Hyfydy) {
        std::cout << "Using HyfydyMuscle model\n";
    } else {
        std::cout << "Using PointPathMuscle model\n";
    }
    
    // Set up path type based on selection
    if (pathType == PointBased) {
        std::cout << "Using PointBasedPath\n";
    } else {
        std::cout << "Using Scholz2015GeometryPath\n";
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

    BallJoint* lumbar = new BallJoint("lumbar", 
            *pelvis, SimTK::Vec3(0, 0.05, 0), SimTK::Vec3(0), 
            *torso, SimTK::Vec3(0, -0.25, 0), SimTK::Vec3(0));

    BallJoint* leftHip = new BallJoint("hip_l", 
            *pelvis, SimTK::Vec3(0, -0.0661, -0.0835), SimTK::Vec3(0), 
            *leftThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0));

    PinJoint* leftKnee = new PinJoint("knee_l", 
            *leftThigh, SimTK::Vec3(0, -0.226, 0), SimTK::Vec3(0),
            *leftShank, SimTK::Vec3(0, 0.1867, 0), SimTK::Vec3(0));

    PinJoint* leftAnkle = new PinJoint("ankle_l",
            *leftShank, SimTK::Vec3(0, -0.2433, 0), SimTK::Vec3(0),
            *leftFoot, SimTK::Vec3(-0.05123, 0.01195, 0.00792), SimTK::Vec3(0));

    BallJoint* rightHip = new BallJoint("hip_r",
            *pelvis, SimTK::Vec3(0, -0.0661, 0.0835), SimTK::Vec3(0),
            *rightThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0));

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
    if (usePointPathMuscle) {
        createPointPathMuscles(model, pelvis, torso, leftThigh, leftShank, 
                              leftFoot, rightThigh, rightShank, rightFoot, 
                              controlValue, true);
    } else {
        // Use template-based muscle creation based on selected model
        if (muscleModel == Hyfydy) {
            createMuscles<HyfydyMuscle>(model, pelvis, torso, leftThigh, 
                                       leftShank, leftFoot, rightThigh, 
                                       rightShank, rightFoot, pathType);
        } else if (muscleModel == DeGrooteFregly) {
            createMuscles<DeGrooteFregly2016Muscle>(model, pelvis, torso, 
                                                   leftThigh, leftShank, leftFoot, 
                                                   rightThigh, rightShank, 
                                                   rightFoot, pathType);
        }
    }

    // Joint damping
    // -------------
    SimTK::Real damping = 1.0;
    MobilityLinearDamper* lumbarDamperX = 
        new MobilityLinearDamper("lumbar_coord_0", 10.0*damping);
    lumbarDamperX->setName("lumbar_coord_0_damper");
    model.addForce(lumbarDamperX);

    MobilityLinearDamper* lumbarDamperY = 
        new MobilityLinearDamper("lumbar_coord_1", 10.0*damping);
    lumbarDamperY->setName("lumbar_coord_1_damper");
    model.addForce(lumbarDamperY);

    MobilityLinearDamper* lumbarDamperZ = 
        new MobilityLinearDamper("lumbar_coord_2", 10.0*damping);
    lumbarDamperZ->setName("lumbar_coord_2_damper");
    model.addForce(lumbarDamperZ);

    MobilityLinearDamper* leftHipDamperX = 
        new MobilityLinearDamper("hip_l_coord_0", damping);
    leftHipDamperX->setName("hip_l_coord_0_damper");
    model.addForce(leftHipDamperX);

    MobilityLinearDamper* leftHipDamperY = 
        new MobilityLinearDamper("hip_l_coord_1", damping);
    leftHipDamperY->setName("hip_l_coord_1_damper");
    model.addForce(leftHipDamperY);

    MobilityLinearDamper* leftHipDamperZ = 
        new MobilityLinearDamper("hip_l_coord_2", damping);
    leftHipDamperZ->setName("hip_l_coord_2_damper");
    model.addForce(leftHipDamperZ);

    MobilityLinearDamper* rightHipDamperX = 
        new MobilityLinearDamper("hip_r_coord_0", damping);
    rightHipDamperX->setName("hip_r_coord_0_damper");
    model.addForce(rightHipDamperX);

    MobilityLinearDamper* rightHipDamperY = 
        new MobilityLinearDamper("hip_r_coord_1", damping);
    rightHipDamperY->setName("hip_r_coord_1_damper");
    model.addForce(rightHipDamperY);

    MobilityLinearDamper* rightHipDamperZ = 
        new MobilityLinearDamper("hip_r_coord_2", damping);
    rightHipDamperZ->setName("hip_r_coord_2_damper");
    model.addForce(rightHipDamperZ);

    MobilityLinearDamper* leftKneeDamper = 
        new MobilityLinearDamper("knee_l_coord_0", damping);
    leftKneeDamper->setName("knee_l_coord_0_damper");
    model.addForce(leftKneeDamper);

    MobilityLinearDamper* leftAnkleDamper = 
        new MobilityLinearDamper("ankle_l_coord_0", damping);
    leftAnkleDamper->setName("ankle_l_coord_0_damper");
    model.addForce(leftAnkleDamper);

    MobilityLinearDamper* rightKneeDamper = 
        new MobilityLinearDamper("knee_r_coord_0", damping);
    rightKneeDamper->setName("knee_r_coord_0_damper");
    model.addForce(rightKneeDamper);

    MobilityLinearDamper* rightAnkleDamper = 
        new MobilityLinearDamper("ankle_r_coord_0", damping);
    rightAnkleDamper->setName("ankle_r_coord_0_damper");
    model.addForce(rightAnkleDamper);

    // Joint stops
    // ----------
    MobilityLinearStop* leftKneeStop = 
        new MobilityLinearStop("knee_l_coord_0", 500, 2.95953, 
                              SimTK::convertDegreesToRadians(-120.0), 
                              SimTK::convertDegreesToRadians(-3.0));
    leftKneeStop->setName("knee_l_coord_0_stop");
    model.addForce(leftKneeStop);

    MobilityLinearStop* rightKneeStop = 
        new MobilityLinearStop("knee_r_coord_0", 500, 2.95953, 
                              SimTK::convertDegreesToRadians(-120.0), 
                              SimTK::convertDegreesToRadians(-3.0));
    rightKneeStop->setName("knee_r_coord_0_stop");
    model.addForce(rightKneeStop);

    MobilityLinearStop* leftAnkleStopX = 
        new MobilityLinearStop("ankle_l_coord_0", 500, 1.41762, 
                              SimTK::convertDegreesToRadians(-60.0), 
                              SimTK::convertDegreesToRadians(25.0));
    leftAnkleStopX->setName("ankle_l_coord_0_stop");
    model.addForce(leftAnkleStopX);

    MobilityLinearStop* rightAnkleStopX = 
        new MobilityLinearStop("ankle_r_coord_0", 500, 1.41762, 
                              SimTK::convertDegreesToRadians(-60.0), 
                              SimTK::convertDegreesToRadians(25.0));
    rightAnkleStopX->setName("ankle_r_coord_0_stop");
    model.addForce(rightAnkleStopX);


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
    
    // // Torso contact points (around the body)
    // addContact(model, "torso_front_contact", torso, 
    //            SimTK::Vec3(0, 0.1, 0), transform, params);
    // addContact(model, "torso_back_contacpt", torso, 
    //            SimTK::Vec3(0, -0.1, 0), transform, params);
    // addContact(model, "torso_left_contact", torso, 
    //            SimTK::Vec3(0, 0, 0.05), transform, params);
    // addContact(model, "torso_right_contact", torso, 
    //            SimTK::Vec3(0, 0, -0.05), transform, params);
    
    // // Pelvis contact points (around the pelvis)
    // addContact(model, "pelvis_front_contact", pelvis, 
    //            SimTK::Vec3(0, 0.05, 0), transform, params);
    // addContact(model, "pelvis_back_contact", pelvis, 
    //            SimTK::Vec3(0, -0.05, 0), transform, params);
    // addContact(model, "pelvis_left_contact", pelvis, 
    //            SimTK::Vec3(0, 0, 0.06), transform, params);
    // addContact(model, "pelvis_right_contact", pelvis, 
    //            SimTK::Vec3(0, 0, -0.06), transform, params);
    
    // // Left hip contact points
    // addContact(model, "left_hip_front_contact", leftThigh, 
    //            SimTK::Vec3(0, 0.15, 0), transform, params);
    // addContact(model, "left_hip_back_contact", leftThigh, 
    //            SimTK::Vec3(0, 0.19, 0), transform, params);
    // addContact(model, "left_hip_medial_contact", leftThigh, 
    //            SimTK::Vec3(0, 0.17, 0.02), transform, params);
    // addContact(model, "left_hip_lateral_contact", leftThigh, 
    //            SimTK::Vec3(0, 0.17, -0.02), transform, params);
    
    // // Right hip contact points
    // addContact(model, "right_hip_front_contact", rightThigh, 
    //            SimTK::Vec3(0, 0.15, 0), transform, params);
    // addContact(model, "right_hip_back_contact", rightThigh, 
    //            SimTK::Vec3(0, 0.19, 0), transform, params);
    // addContact(model, "right_hip_medial_contact", rightThigh, 
    //            SimTK::Vec3(0, 0.17, -0.02), transform, params);
    // addContact(model, "right_hip_lateral_contact", rightThigh, 
    //            SimTK::Vec3(0, 0.17, 0.02), transform, params);
    
    // // Left knee contact points
    // addContact(model, "left_knee_front_contact", leftShank, 
    //            SimTK::Vec3(0, 0.15, 0), transform, params);
    // addContact(model, "left_knee_back_contact", leftShank, 
    //            SimTK::Vec3(0, 0.22, 0), transform, params);
    // addContact(model, "left_knee_medial_contact", leftShank, 
    //            SimTK::Vec3(0, 0.1867, 0.015), transform, params);
    // addContact(model, "left_knee_lateral_contact", leftShank, 
    //            SimTK::Vec3(0, 0.1867, -0.015), transform, params);
    
    // // Right knee contact points
    // addContact(model, "right_knee_front_contact", rightShank, 
    //            SimTK::Vec3(0, 0.15, 0), transform, params);
    // addContact(model, "right_knee_back_contact", rightShank, 
    //            SimTK::Vec3(0, 0.22, 0), transform, params);
    // addContact(model, "right_knee_medial_contact", rightShank, 
    //            SimTK::Vec3(0, 0.1867, -0.015), transform, params);
    // addContact(model, "right_knee_lateral_contact", rightShank, 
    //            SimTK::Vec3(0, 0.1867, 0.015), transform, params);
    
    // // Left ankle contact points
    // addContact(model, "left_ankle_front_contact", leftFoot, 
    //            SimTK::Vec3(-0.03, 0.012, 0.008), transform, params);
    // addContact(model, "left_ankle_back_contact", leftFoot, 
    //            SimTK::Vec3(-0.07, 0.012, 0.008), transform, params);
    // addContact(model, "left_ankle_medial_contact", leftFoot, 
    //            SimTK::Vec3(-0.05123, 0.012, 0.013), transform, params);
    // addContact(model, "left_ankle_lateral_contact", leftFoot, 
    //            SimTK::Vec3(-0.05123, 0.012, 0.003), transform, params);
    
    // // Right ankle contact points
    // addContact(model, "right_ankle_front_contact", rightFoot, 
    //            SimTK::Vec3(-0.03, 0.012, -0.008), transform, params);
    // addContact(model, "right_ankle_back_contact", rightFoot, 
    //            SimTK::Vec3(-0.07, 0.012, -0.008), transform, params);
    // addContact(model, "right_ankle_medial_contact", rightFoot, 
    //            SimTK::Vec3(-0.05123, 0.012, -0.013), transform, params);
    // addContact(model, "right_ankle_lateral_contact", rightFoot, 
    //            SimTK::Vec3(-0.05123, 0.012, -0.003), transform, params);

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
    controller->setDiscreteControls(state, controls);

    // Simulate
    // --------
    SimTK::Real finalTime = 20.0;
    Manager manager(model);
    manager.setWriteToStorage(visualize);
    manager.setPerformAnalyses(visualize);
    manager.setIntegratorMethod(Manager::IntegratorMethod::CPodes);
    manager.setIntegratorAccuracy(1e-2);
    manager.initialize(state);
    manager.integrate(finalTime);
    
    // Visualize    
    // --------
    if (visualize) {
        TimeSeriesTable table = manager.getStatesTable();
        std::string pathTypeStr;
        if (pathType == Scholz2015) {
            pathTypeStr = "_scholz2015path";
        } else if (pathType == PointBased) {
            pathTypeStr = "_pointbasedpath";
        } else if (pathType == PathType::Geometry) {
            pathTypeStr = "_geometrypath";
        }

        std::string muscleModelStr;
        if (usePointPathMuscle) {
            muscleModelStr = "pointpath";
            pathTypeStr = "";
        } else {
            if (muscleModel == Hyfydy) {
                muscleModelStr = "hyfydy";
            } else {
                muscleModelStr = "degroote";
            }
        }

        std::string statesFileName = fmt::format("states_{}_{}_control{}{}.sto", 
                                                 muscleModelStr, 
                                                 randomizeSpeeds ? "random" : "zero",
                                                 controlValue,
                                                 pathTypeStr);
        STOFileAdapter::write(table, statesFileName);

        VisualizerUtilities::showMotion(model, table);
    }

    return EXIT_SUCCESS;
}
