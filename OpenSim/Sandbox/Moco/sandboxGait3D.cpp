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
#include <OpenSim/Simulation/VisualizerUtilities.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>
#include <OpenSim/Simulation/Model/ExponentialContactForce.h>

using namespace OpenSim;

class MobilityLinearDamper : public Force {
    OpenSim_DECLARE_CONCRETE_OBJECT(MobilityLinearDamper, Force);
public:
    OpenSim_DECLARE_PROPERTY(coordinate, std::string,
            "Coordinate (name) to apply force to.");
    OpenSim_DECLARE_PROPERTY(damping, double, 
            "Damping coefficient (N-m/(rad/s) for rotational, "
            "N/(m/s) for translational)");

    MobilityLinearDamper() {
        constructProperties();
    }

    MobilityLinearDamper(const std::string& coordinateName, double damping) {
        constructProperties();
        set_coordinate(coordinateName);
        set_damping(damping);
    }

private:
    void constructProperties() {
        constructProperty_coordinate("");
        constructProperty_damping(0.0);
    }

    void extendAddToSystem(SimTK::MultibodySystem& system) const override {
        Super::extendAddToSystem(system);

        const std::string& coordName = get_coordinate();
        const Coordinate& coord = getModel().getCoordinateSet().get(coordName);
        SimTK::MobilizerQIndex qIndex = coord.getMobilizerQIndex();
        SimTK::MobilizedBodyIndex mobodIndex = coord.getBodyIndex();
        const SimTK::MobilizedBody& mobod = 
                system.getMatterSubsystem().getMobilizedBody(mobodIndex);

        SimTK::GeneralForceSubsystem& forces = _model->updForceSubsystem();
        SimTK::Force::MobilityLinearDamper damper(forces, mobod,
            qIndex, get_damping());
    }
};

class MobilityLinearStop : public Force {
    OpenSim_DECLARE_CONCRETE_OBJECT(MobilityLinearStop, Force);
public:
    OpenSim_DECLARE_PROPERTY(coordinate, std::string,
        "Coordinate (name) to apply force to.");
    OpenSim_DECLARE_PROPERTY(stiffness, double, 
            "Stiffness coefficient (N-m/rad for rotational, "
            "N/m for translational)");
    OpenSim_DECLARE_PROPERTY(damping, double, 
            "Damping coefficient (N-m/(rad/s) for rotational, "
            "N/(m/s) for translational)");
    OpenSim_DECLARE_PROPERTY(q_high, double, 
            "High position (rad for rotational, m for translational)");
    OpenSim_DECLARE_PROPERTY(q_low, double, 
            "Low position (rad for rotational, m for translational)");

    MobilityLinearStop() {
        constructProperties();
    }

    MobilityLinearStop(const std::string& coordinateName, double stiffness, 
            double damping, double q_low, double q_high) {
        constructProperties();
        set_coordinate(coordinateName);
        set_stiffness(stiffness);
        set_damping(damping); 
        set_q_low(q_low);
        set_q_high(q_high);
    }

private:
    void constructProperties() {
        constructProperty_coordinate("");
        constructProperty_stiffness(0.0);
        constructProperty_damping(0.0);
        constructProperty_q_high(0.0);
        constructProperty_q_low(0.0);
    }

    void extendAddToSystem(SimTK::MultibodySystem& system) const override {
        Super::extendAddToSystem(system);

        const std::string& coordName = get_coordinate();
        const Coordinate& coord = getModel().getCoordinateSet().get(coordName);
        SimTK::MobilizerQIndex qIndex = coord.getMobilizerQIndex();
        SimTK::MobilizedBodyIndex mobodIndex = coord.getBodyIndex();
        const SimTK::MobilizedBody& mobod = 
                system.getMatterSubsystem().getMobilizedBody(mobodIndex);

        SimTK::GeneralForceSubsystem& forces = _model->updForceSubsystem();
        SimTK::Force::MobilityLinearStop stop(forces, mobod,
            qIndex, get_stiffness(), get_damping(), get_q_low(), get_q_high());
    }
};




enum BodyType {LeftFoot=0, RightFoot, LeftShank, RightShank, LeftThigh, 
               RightThigh, Pelvis, Torso};
enum Muscle {GlutMed_R=0, AddMag_R, Hamstrings_R, Bifemsh_R, GlutMax_R, 
             Iliopsoas_R, RectFem_R, Vasti_R, Gastroc_R, Soleus_R, TibAnt_R, 
             GlutMed_L, AddMag_L, Hamstrings_L, Bifemsh_L, GlutMax_L, 
             Iliopsoas_L, RectFem_L, Vasti_L, Gastroc_L, Soleus_L, TibAnt_L};

enum Contact {Heel=0, LateralToe, MedialToe};

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

void addContactGeometry(Body* body, const SimTK::Vec3& offset, 
                        const std::string& name, const double& radius) {
    PhysicalOffsetFrame* offsetFrame = 
            new PhysicalOffsetFrame(name, *body, offset);
    offsetFrame->attachGeometry(new Sphere(radius));
    body->addComponent(offsetFrame); 
}

DeGrooteFregly2016Muscle* addMuscle(Model& model, const std::string& name, 
        double maxIsometricForce, double optimalFiberLength, 
        double tendonSlackLength, double pennationAngle) {
    auto* muscle = new DeGrooteFregly2016Muscle();
    muscle->setName(name);
    muscle->set_max_isometric_force(maxIsometricForce);
    muscle->set_optimal_fiber_length(optimalFiberLength);
    muscle->set_tendon_slack_length(tendonSlackLength);
    muscle->set_pennation_angle_at_optimal(pennationAngle);
    muscle->set_ignore_tendon_compliance(true);
    muscle->set_default_activation(0.1);
    model.addForce(muscle);
    return muscle;
}

void addPath(DeGrooteFregly2016Muscle* muscle, const std::string& name, 
             const std::vector<std::pair<PhysicalFrame*, SimTK::Vec3>>& points) {
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
int main() {
    // Logger::setLevel(Logger::Level::Warn);
    
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

    BallJoint* leftHip = new BallJoint("leftHip", 
            *pelvis, SimTK::Vec3(0, -0.0661, -0.0835), SimTK::Vec3(0), 
            *leftThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0));

    PinJoint* leftKnee = new PinJoint("leftKnee", 
            *leftThigh, SimTK::Vec3(0, -0.226, 0), SimTK::Vec3(0),
            *leftShank, SimTK::Vec3(0, 0.1867, 0), SimTK::Vec3(0));

    PinJoint* leftAnkle = new PinJoint("leftAnkle",
            *leftShank, SimTK::Vec3(0, -0.2433, 0), SimTK::Vec3(0),
            *leftFoot, SimTK::Vec3(-0.05123, 0.01195, 0.00792), SimTK::Vec3(0));

    BallJoint* rightHip = new BallJoint("rightHip",
            *pelvis, SimTK::Vec3(0, -0.0661, 0.0835), SimTK::Vec3(0),
            *rightThigh, SimTK::Vec3(0, 0.17, 0), SimTK::Vec3(0));

    PinJoint* rightKnee = new PinJoint("rightKnee",
            *rightThigh, SimTK::Vec3(0, -0.226, 0), SimTK::Vec3(0),
            *rightShank, SimTK::Vec3(0, 0.1867, 0), SimTK::Vec3(0));

    PinJoint* rightAnkle = new PinJoint("rightAnkle",
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
    auto* glut_med_r = addMuscle(model, "glut_med_r", 2045, 0.0733, 0.066, 0.3578);
    std::vector<std::pair<PhysicalFrame*, SimTK::Vec3>> points;
    points.emplace_back(pelvis, SimTK::Vec3(-0.0148, 0.0445, 0.0766));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0258, 0.1642, 0.0527));
    addPath(glut_med_r, "glut_med_r_path", points);

    // add_mag_r
    auto* add_mag_r = addMuscle(model, "add_mag_r", 2268, 0.087, 0.06, 0.0872665);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0025, -0.1174, 0.0255));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0045, 0.0489, 0.0339));
    addPath(add_mag_r, "add_mag_r_path", points);

    // hamstrings_r
    auto* hamstrings_r = addMuscle(model, "hamstrings_r", 2594, 0.0976, 0.319, 0.2025);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.05526, -0.10257, 0.06944));
    points.emplace_back(rightShank, SimTK::Vec3(-0.028, 0.1667, 0.02943));
    points.emplace_back(rightShank, SimTK::Vec3(-0.021, 0.1467, 0.0343));
    addPath(hamstrings_r, "hamstrings_r_path", points);

    // bifemsh_r
    auto* bifemsh_r = addMuscle(model, "bifemsh_r", 804, 0.1103, 0.095, 0.2147);
    points.clear();
    points.emplace_back(rightThigh, SimTK::Vec3(0.005, -0.0411, 0.0234));
    points.emplace_back(rightShank, SimTK::Vec3(-0.028, 0.1667, 0.02943));
    points.emplace_back(rightShank, SimTK::Vec3(-0.021, 0.1467, 0.0343));
    addPath(bifemsh_r, "bifemsh_r_path", points);

    // glut_max_r
    auto* glut_max_r = addMuscle(model, "glut_max_r", 1944, 0.1569, 0.111, 0.3822);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0642, 0.0176, 0.0563));
    points.emplace_back(pelvis, SimTK::Vec3(-0.0669, -0.052, 0.0914));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0426, 0.117, 0.0293));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0156, 0.0684, 0.0419));
    addPath(glut_max_r, "glut_max_r_path", points);

    // iliopsoas_r
    auto* iliopsoas_r = addMuscle(model, "iliopsoas_r", 2186, 0.1066, 0.152, 0.2496);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.006, 0.0887, 0.0289));
    points.emplace_back(pelvis, SimTK::Vec3(0.0407, -0.01, 0.076));
    points.emplace_back(rightThigh, SimTK::Vec3(0.033, 0.135, 0.0038));
    points.emplace_back(rightThigh, SimTK::Vec3(-0.0188, 0.1103, 0.0104));
    addPath(iliopsoas_r, "iliopsoas_r_path", points);

    // rect_fem_r
    auto* rect_fem_r = addMuscle(model, "rect_fem_r", 1169, 0.0759, 0.3449, 0.2426);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.0412, -0.0311, 0.0968));
    points.emplace_back(rightThigh, SimTK::Vec3(0.038, -0.17, 0.004));
    points.emplace_back(rightShank, SimTK::Vec3(0.038, 0.2117, 0.0018));
    addPath(rect_fem_r, "rect_fem_r_path", points);

    // vasti_r
    auto* vasti_r = addMuscle(model, "vasti_r", 4530, 0.0993, 0.1231, 0.0785);
    points.clear();
    points.emplace_back(rightThigh, SimTK::Vec3(0.029, -0.0224, 0.031));
    points.emplace_back(rightThigh, SimTK::Vec3(0.038, -0.17, 0.007));
    points.emplace_back(rightShank, SimTK::Vec3(0.038, 0.2117, 0.0018));
    addPath(vasti_r, "vasti_r_path", points);

    // gastroc_r
    auto* gastroc_r = addMuscle(model, "gastroc_r", 2241, 0.051, 0.384, 0.1728);
    points.clear();
    points.emplace_back(rightThigh, SimTK::Vec3(-0.02, -0.218, -0.024));
    points.emplace_back(rightFoot, SimTK::Vec3(-0.095, 0.001, -0.0053));
    addPath(gastroc_r, "gastroc_r_path", points);

    // soleus_r
    auto* soleus_r = addMuscle(model, "soleus_r", 3549, 0.044, 0.248, 0.4939);
    points.clear();
    points.emplace_back(rightShank, SimTK::Vec3(-0.0024, 0.0334, 0.0071));
    points.emplace_back(rightFoot, SimTK::Vec3(-0.095, 0.001, -0.0053));
    addPath(soleus_r, "soleus_r_path", points);

    // tib_ant_r
    auto* tib_ant_r = addMuscle(model, "tib_ant_r", 1579, 0.0683, 0.243, 0.1676);
    points.clear();
    points.emplace_back(rightShank, SimTK::Vec3(0.0179, 0.0243, 0.0115));
    points.emplace_back(rightShank, SimTK::Vec3(0.0329, -0.2084, -0.0177));
    points.emplace_back(rightFoot, SimTK::Vec3(0.0166, -0.0122, -0.0305));
    addPath(tib_ant_r, "tib_ant_r_path", points);

    // glut_med_l
    auto* glut_med_l = addMuscle(model, "glut_med_l", 2045, 0.0733, 0.066, 0.3578);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0148, 0.0445, -0.0766));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0258, 0.1642, -0.0527));
    addPath(glut_med_l, "glut_med_l_path", points);

    // add_mag_l
    auto* add_mag_l = addMuscle(model, "add_mag_l", 2268, 0.087, 0.06, 0.0872665);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0025, -0.1174, -0.0255));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0045, 0.0489, -0.0339));
    addPath(add_mag_l, "add_mag_l_path", points);

    // hamstrings_l
    auto* hamstrings_l = addMuscle(model, "hamstrings_l", 2594, 0.0976, 0.319, 0.2025);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.05526, -0.10257, -0.06944));
    points.emplace_back(leftShank, SimTK::Vec3(-0.028, 0.1667, -0.02943));
    points.emplace_back(leftShank, SimTK::Vec3(-0.021, 0.1467, -0.0343));
    addPath(hamstrings_l, "hamstrings_l_path", points);

    // bifemsh_l
    auto* bifemsh_l = addMuscle(model, "bifemsh_l", 804, 0.1103, 0.095, 0.2147);
    points.clear();
    points.emplace_back(leftThigh, SimTK::Vec3(0.005, -0.0411, -0.0234));
    points.emplace_back(leftShank, SimTK::Vec3(-0.028, 0.1667, -0.02943));
    points.emplace_back(leftShank, SimTK::Vec3(-0.021, 0.1467, -0.0343));
    addPath(bifemsh_l, "bifemsh_l_path", points);

    // glut_max_l
    auto* glut_max_l = addMuscle(model, "glut_max_l", 1944, 0.1569, 0.111, 0.3822);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(-0.0642, 0.0176, -0.0563));
    points.emplace_back(pelvis, SimTK::Vec3(-0.0669, -0.052, -0.0914));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0426, 0.117, -0.0293));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0156, 0.0684, -0.0419));
    addPath(glut_max_l, "glut_max_l_path", points);

    // iliopsoas_l
    auto* iliopsoas_l = addMuscle(model, "iliopsoas_l", 2186, 0.1066, 0.152, 0.2496);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.006, 0.0887, -0.0289));
    points.emplace_back(pelvis, SimTK::Vec3(0.0407, -0.01, -0.076));
    points.emplace_back(leftThigh, SimTK::Vec3(0.033, 0.135, -0.0038));
    points.emplace_back(leftThigh, SimTK::Vec3(-0.0188, 0.1103, -0.0104));
    addPath(iliopsoas_l, "iliopsoas_l_path", points);

    // rect_fem_l
    auto* rect_fem_l = addMuscle(model, "rect_fem_l", 1169, 0.0759, 0.3449, 0.2426);
    points.clear();
    points.emplace_back(pelvis, SimTK::Vec3(0.0412, -0.0311, -0.0968));
    points.emplace_back(leftThigh, SimTK::Vec3(0.038, -0.17, -0.004));
    points.emplace_back(leftShank, SimTK::Vec3(0.038, 0.2117, -0.0018));
    addPath(rect_fem_l, "rect_fem_l_path", points);

    // vasti_l
    auto* vasti_l = addMuscle(model, "vasti_l", 4530, 0.0993, 0.1231, 0.0785);
    points.clear();
    points.emplace_back(leftThigh, SimTK::Vec3(0.029, -0.0224, -0.031));
    points.emplace_back(leftThigh, SimTK::Vec3(0.038, -0.17, -0.007));
    points.emplace_back(leftShank, SimTK::Vec3(0.038, 0.2117, -0.0018));
    addPath(vasti_l, "vasti_l_path", points);

    // gastroc_l
    auto* gastroc_l = addMuscle(model, "gastroc_l", 2241, 0.051, 0.384, 0.1728);
    points.clear();
    points.emplace_back(leftThigh, SimTK::Vec3(-0.02, -0.218, 0.024));
    points.emplace_back(leftFoot, SimTK::Vec3(-0.095, 0.001, 0.0053));
    addPath(gastroc_l, "gastroc_l_path", points);

    // soleus_l
    auto* soleus_l = addMuscle(model, "soleus_l", 3549, 0.044, 0.248, 0.4939);
    points.clear();
    points.emplace_back(leftShank, SimTK::Vec3(-0.0024, 0.0334, -0.0071));
    points.emplace_back(leftFoot, SimTK::Vec3(-0.095, 0.001, 0.0053));
    addPath(soleus_l, "soleus_l_path", points);

    // tib_ant_l
    auto* tib_ant_l = addMuscle(model, "tib_ant_l", 1579, 0.0683, 0.243, 0.1676);
    points.clear();
    points.emplace_back(leftShank, SimTK::Vec3(0.0179, 0.0243, -0.0115));
    points.emplace_back(leftShank, SimTK::Vec3(0.0329, -0.2084, 0.0177));
    points.emplace_back(leftFoot, SimTK::Vec3(0.0166, -0.0122, 0.0305));
    addPath(tib_ant_l, "tib_ant_l_path", points);

    // Joint damping
    // -------------
    SimTK::Real damping = 0.1;
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
        new MobilityLinearDamper("leftHip_coord_0", damping);
    leftHipDamperX->setName("leftHip_coord_0_damper");
    model.addForce(leftHipDamperX);

    MobilityLinearDamper* leftHipDamperY = 
        new MobilityLinearDamper("leftHip_coord_1", damping);
    leftHipDamperY->setName("leftHip_coord_1_damper");
    model.addForce(leftHipDamperY);

    MobilityLinearDamper* leftHipDamperZ = 
        new MobilityLinearDamper("leftHip_coord_2", damping);
    leftHipDamperZ->setName("leftHip_coord_2_damper");
    model.addForce(leftHipDamperZ);

    MobilityLinearDamper* rightHipDamperX = 
        new MobilityLinearDamper("rightHip_coord_0", damping);
    rightHipDamperX->setName("rightHip_coord_0_damper");
    model.addForce(rightHipDamperX);

    MobilityLinearDamper* rightHipDamperY = 
        new MobilityLinearDamper("rightHip_coord_1", damping);
    rightHipDamperY->setName("rightHip_coord_1_damper");
    model.addForce(rightHipDamperY);

    MobilityLinearDamper* rightHipDamperZ = 
        new MobilityLinearDamper("rightHip_coord_2", damping);
    rightHipDamperZ->setName("rightHip_coord_2_damper");
    model.addForce(rightHipDamperZ);

    MobilityLinearDamper* leftKneeDamper = 
        new MobilityLinearDamper("leftKnee_coord_0", damping);
    leftKneeDamper->setName("leftKnee_coord_0_damper");
    model.addForce(leftKneeDamper);

    MobilityLinearDamper* leftAnkleDamper = 
        new MobilityLinearDamper("leftAnkle_coord_0", damping);
    leftAnkleDamper->setName("leftAnkle_coord_0_damper");
    model.addForce(leftAnkleDamper);

    MobilityLinearDamper* rightKneeDamper = 
        new MobilityLinearDamper("rightKnee_coord_0", damping);
    rightKneeDamper->setName("rightKnee_coord_0_damper");
    model.addForce(rightKneeDamper);

    MobilityLinearDamper* rightAnkleDamper = 
        new MobilityLinearDamper("rightAnkle_coord_0", damping);
    rightAnkleDamper->setName("rightAnkle_coord_0_damper");
    model.addForce(rightAnkleDamper);

    // Joint stops
    // ----------
    // MobilityLinearStop* leftHipStopX = 
    //     new MobilityLinearStop("leftHip_coord_0", 20, 1.22629, 
    //     SimTK::convertDegreesToRadians(-20.0), 
    //     SimTK::convertDegreesToRadians(45.0));
    // leftHipStopX->setName("leftHip_coord_0_stop");
    // model.addForce(leftHipStopX);

    // MobilityLinearStop* leftHipStopZ = 
    //     new MobilityLinearStop("leftHip_coord_2", 20, 1.22629, 
    //     SimTK::convertDegreesToRadians(-20.0), 
    //     SimTK::convertDegreesToRadians(45.0));
    // leftHipStopZ->setName("leftHip_coord_2_stop");
    // model.addForce(leftHipStopZ);

    // MobilityLinearStop* rightHipStopX = 
    //     new MobilityLinearStop("rightHip_coord_0", 20, 1.22629, 
    //     SimTK::convertDegreesToRadians(-20.0), 
    //     SimTK::convertDegreesToRadians(45.0));
    // rightHipStopX->setName("rightHip_coord_0_stop");
    // model.addForce(rightHipStopX);
    
    // MobilityLinearStop* rightHipStopZ = 
    //     new MobilityLinearStop("rightHip_coord_2", 20, 1.22629, 
    //     SimTK::convertDegreesToRadians(-20.0), 
    //     SimTK::convertDegreesToRadians(45.0));
    // rightHipStopZ->setName("rightHip_coord_2_stop");
    // model.addForce(rightHipStopZ);

    MobilityLinearStop* leftKneeStop = 
        new MobilityLinearStop("leftKnee_coord_0", 500, 2.95953, 
        SimTK::convertDegreesToRadians(-120.0), 
        SimTK::convertDegreesToRadians(-3.0));
    leftKneeStop->setName("leftKnee_coord_0_stop");
    model.addForce(leftKneeStop);

    MobilityLinearStop* rightKneeStop = 
        new MobilityLinearStop("rightKnee_coord_0", 500, 2.95953, 
        SimTK::convertDegreesToRadians(-120.0), 
        SimTK::convertDegreesToRadians(-3.0));
    rightKneeStop->setName("rightKnee_coord_0_stop");
    model.addForce(rightKneeStop);

    MobilityLinearStop* leftAnkleStopX = 
        new MobilityLinearStop("leftAnkle_coord_0", 500, 1.41762, 
        SimTK::convertDegreesToRadians(-60.0), 
        SimTK::convertDegreesToRadians(25.0));
    leftAnkleStopX->setName("leftAnkle_coord_0_stop");
    model.addForce(leftAnkleStopX);

    MobilityLinearStop* rightAnkleStopX = 
        new MobilityLinearStop("rightAnkle_coord_0", 500, 1.41762, 
        SimTK::convertDegreesToRadians(-60.0), 
        SimTK::convertDegreesToRadians(25.0));
    rightAnkleStopX->setName("rightAnkle_coord_0_stop");
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

    // Construct system
    // ---------------
    model.finalizeConnections();
    SimTK::State state = model.initSystem();

    // model.printSubcomponentInfo();

    // Default state
    // -------------
    state.updQ()[4] = 1.05;
    // model.getVisualizer().show(state);

    // Simulate
    // --------
    SimTK::Real finalTime = 20.0;
    Manager manager(model);
    manager.setWriteToStorage(false);
    manager.setPerformAnalyses(false);
    manager.setIntegratorMethod(Manager::IntegratorMethod::CPodes);
    manager.setIntegratorAccuracy(1e-2);
    manager.initialize(state);
    SimTK::Real cpuTime = SimTK::cpuTime();
    SimTK::Real realTime = SimTK::realTime();   
    manager.integrate(finalTime);
    cpuTime = SimTK::cpuTime() - cpuTime;
    realTime = SimTK::realTime() - realTime;
    std::cout << "CPU time: " << cpuTime << " seconds" << std::endl;
    std::cout << "Real time: " << realTime << " seconds" << std::endl;
    std::cout << "Real time factor: " << finalTime / realTime << std::endl;
    
    // Visualize    
    // --------
    // TimeSeriesTable table = manager.getStatesTable();
    // VisualizerUtilities::showMotion(model, table);

    return EXIT_SUCCESS;
}
