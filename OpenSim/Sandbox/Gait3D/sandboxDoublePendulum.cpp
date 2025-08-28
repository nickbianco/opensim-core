/* -------------------------------------------------------------------------- *
 * OpenSim Moco: sandboxDoublePendulum.cpp                                    *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2019 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Nicholas Bianco                                                 *
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

#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>
#include <OpenSim/Simulation/VisualizerUtilities.h>
#include <OpenSim/Simulation/Control/Controller.h>
#include <OpenSim/Simulation/Model/ContactGeometry.h>
#include <OpenSim/Simulation/Model/ExponentialContactForce.h>
#include <OpenSim/Actuators/ModelFactory.h>
#include <OpenSim/Actuators/CoordinateActuator.h>
#include <OpenSim/Common/Sine.h>

using namespace OpenSim;

class DiscreteController : public Controller {
    OpenSim_DECLARE_CONCRETE_OBJECT(DiscreteController, Controller);
public:
    DiscreteController() = default;
    void setDiscreteControls(SimTK::State& s,
            const SimTK::Vector& controls) const;
    SimTK::Vector& updDiscreteControls(SimTK::State& s) const;
    const SimTK::Vector& getDiscreteControls(const SimTK::State& s) const;
    void computeControls(
            const SimTK::State& s, SimTK::Vector& controls) const override;
protected:
    void extendRealizeTopology(SimTK::State&) const override;
    mutable SimTK::DiscreteVariableIndex m_discreteVarIndex;

};

void DiscreteController::setDiscreteControls(SimTK::State& s,
        const SimTK::Vector& controls) const {
    updDiscreteControls(s) = controls;
}

SimTK::Vector& DiscreteController::updDiscreteControls(SimTK::State& s) const {
    const SimTK::Subsystem& subSys = getSystem().getDefaultSubsystem();
    auto& dv = subSys.updDiscreteVariable(s, m_discreteVarIndex);
    auto& discreteControls = SimTK::Value<SimTK::Vector>::updDowncast(dv).upd();
    return discreteControls;
}

const SimTK::Vector& DiscreteController::getDiscreteControls(
        const SimTK::State& s) const {
    const SimTK::Subsystem& subSys = getSystem().getDefaultSubsystem();
    auto& dv = subSys.getDiscreteVariable(s, m_discreteVarIndex);
    auto& discreteControls = SimTK::Value<SimTK::Vector>::downcast(dv).get();
    return discreteControls;
}

void DiscreteController::computeControls(
        const SimTK::State& s, SimTK::Vector& controls) const {
    const SimTK::Subsystem& subSys = getSystem().getDefaultSubsystem();
    const auto& dv = subSys.getDiscreteVariable(s, m_discreteVarIndex) ;
    const auto& discreteControls =
            SimTK::Value<SimTK::Vector>::downcast(dv).get();
    controls += discreteControls;
}

void DiscreteController::extendRealizeTopology(SimTK::State& state) const {
    Super::extendRealizeTopology(state);
    const SimTK::Subsystem& subSys = getSystem().getDefaultSubsystem();
    m_discreteVarIndex =
            subSys.allocateDiscreteVariable(state, SimTK::Stage::Dynamics,
                    new SimTK::Value<SimTK::Vector>(
                            SimTK::Vector(getModel().getNumControls(), 0.0)));
}

void addContact(Model& model, const std::string& name, PhysicalFrame* frame, 
               const SimTK::Vec3& location, const SimTK::Transform& transform, 
               const SimTK::ExponentialSpringParameters& params) {
    auto* contact = new ExponentialContactForce(transform, *frame, location, 
        params);
    contact->setName(name);
    model.addForce(contact);
}


Model createNLinkPendulum(int numLinks) {
    Model model;
    OPENSIM_THROW_IF(numLinks < 0, Exception, "numLinks must be nonnegative.");
    std::string name;
    if (numLinks == 0) {
        name = "empty_model";
    } else if (numLinks == 1) {
        name = "pendulum";
    } else if (numLinks == 2) {
        name = "double_pendulum";
    } else {
        name = std::to_string(numLinks) + "_link_pendulum";
    }
    model.setName(name);
    const auto& ground = model.getGround();

    using SimTK::Inertia;
    using SimTK::Vec3;

    auto* root = new OpenSim::Body("root", 1, Vec3(0), Inertia(1));
    model.addBody(root);

    // Assume each body is 1 m long.
    auto* free = new FreeJoint("free", ground, Vec3(0), Vec3(0), *root,
            Vec3(0, 0, 0), Vec3(0));
    model.addJoint(free);
    const PhysicalFrame* prevBody = root;

    Ellipsoid bodyGeometry(0.5, 0.1, 0.1);
    bodyGeometry.setColor(SimTK::Gray);

    for (int i = 0; i < numLinks; ++i) {
        const std::string istr = std::to_string(i);
        auto* bi = new OpenSim::Body("b" + istr, 1, Vec3(0), Inertia(1));
        model.addBody(bi);

        // Assume each body is 1 m long.
        auto* ji = new PinJoint("j" + istr, *prevBody, Vec3(0), Vec3(0), *bi,
                Vec3(-1, 0, 0), Vec3(0));
        auto& qi = ji->updCoordinate();
        qi.setName("q" + istr);
        model.addJoint(ji);

        // auto* taui = new CoordinateActuator();
        // taui->setCoordinate(&ji->updCoordinate());
        // taui->setName("tau" + istr);
        // taui->setOptimalForce(1);
        // model.addComponent(taui);

        auto* marker = new Marker("marker" + istr, *bi, Vec3(0));
        model.addMarker(marker);

        // Attach an ellipsoid to a frame located at the center of each body.
        PhysicalOffsetFrame* bicenter = new PhysicalOffsetFrame(
                "b" + istr + "center", *bi, SimTK::Transform(Vec3(-0.5, 0, 0)));
        bi->addComponent(bicenter);
        bicenter->attachGeometry(bodyGeometry.clone());

        prevBody = bi;
    }

    model.finalizeConnections();

    return model;
}

void simulateScholz2015GeometryPath() {

    Model model = createNLinkPendulum(2);
    model.setUseVisualizer(true);

    // Create a PathActuator with a Scholz2015GeometryPath.
    auto* actu = new PathActuator();
    actu->set_path(Scholz2015GeometryPath());
    model.addComponent(actu);   

    // Add contact forces.
    SimTK::Transform transform(SimTK::Rotation(-0.5*SimTK::Pi, SimTK::XAxis), 
                              SimTK::Vec3(0));
    SimTK::ExponentialSpringParameters params;
    params.setNormalViscosity(1.0);
    params.setInitialMuStatic(0.9);
    params.setInitialMuKinetic(0.6);
    params.setSettleVelocity(0.1);
    addContact(model, "contact0", &model.updComponent<Body>("/bodyset/root"), 
               SimTK::Vec3(0), transform, params);
    // addContact(model, "contact1", &model.updComponent<Body>("/bodyset/b0"), 
    //            SimTK::Vec3(0), transform, params);
    // addContact(model, "contact2", &model.updComponent<Body>("/bodyset/b1"), 
    //            SimTK::Vec3(0), transform, params);

    // Set the path's origin and insertion.
    Scholz2015GeometryPath& path = actu->updPath<Scholz2015GeometryPath>();
    path.setOrigin(model.getComponent<Body>("/bodyset/root"), SimTK::Vec3(-0.5, 0, 0));
    path.setInsertion(model.getComponent<Body>("/bodyset/b0"), 
            SimTK::Vec3(-0.5, 0, 0));

    auto* obstacle = new ContactEllipsoid(SimTK::Vec3(0.1, 0.1, 0.3),
        SimTK::Vec3(0), SimTK::Vec3(0), model.getComponent<Body>("/bodyset/root"));
    // auto* obstacle = new ContactCylinder(0.1,
    //     SimTK::Vec3(0), SimTK::Vec3(0), model.getComponent<Body>("/bodyset/root"));
    model.addComponent(obstacle);
    path.addObstacle(*obstacle, SimTK::Vec3(0.0, 0.1, 0.0));

    // Add a discrete controller to the model.
    DiscreteController* controller = new DiscreteController();
    controller->setName("controller");
    model.addController(controller);

    // Initialize the system.
    SimTK::State state = model.initSystem();
    SimTK::Vector controls(model.getNumControls(), 10);
    controller->setDiscreteControls(state, controls);
    // state.updQ()[4] = 0.;

    // Simulate.
    Manager manager(model);
    // manager.setIntegratorMethod(Manager::IntegratorMethod::SemiExplicitEuler2);
    manager.setIntegratorMaximumStepSize(1e-3);
    manager.setIntegratorAccuracy(1e-3);
    manager.initialize(state);

    auto start = std::chrono::high_resolution_clock::now();
    manager.integrate(20.0);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Integration took " << elapsed.count() << " seconds." << std::endl;
 
    // TimeSeriesTable table = manager.getStatesTable();
    // VisualizerUtilities::showMotion(model, table);
}

int main() {
    simulateScholz2015GeometryPath();
    return 0;
}

