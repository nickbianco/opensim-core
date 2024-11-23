/* -------------------------------------------------------------------------- *
 * OpenSim Moco: StationPlaneContactForce.cpp                                 *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2017 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Nicholas Bianco, Christopher Dembia                             *
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

#include "StationPlaneContactForce.h"

using namespace OpenSim;

OpenSim::Array<std::string> StationPlaneContactForce::getRecordLabels() const {
    OpenSim::Array<std::string> labels;
    const auto stationName = getConnectee("station").getName();
    labels.append(getName() + "." + stationName + ".force.X");
    labels.append(getName() + "." + stationName + ".force.Y");
    labels.append(getName() + "." + stationName + ".force.Z");
    return labels;
}

OpenSim::Array<double> StationPlaneContactForce::getRecordValues(
        const SimTK::State& s) const {
    OpenSim::Array<double> values;
    const SimTK::Vec3& force = getContactForceOnStation(s);
    values.append(force[0]);
    values.append(force[1]);
    values.append(force[2]);
    return values;
}

void StationPlaneContactForce::generateDecorations(bool fixed, 
        const ModelDisplayHints& hints, const SimTK::State& s,
        SimTK::Array_<SimTK::DecorativeGeometry>& geoms) const {
    Super::generateDecorations(fixed, hints, s, geoms);

    if (!fixed) {
        getModel().realizeVelocity(s);
        // Normalize contact force vector by body weight so that the line
        // is 1 meter long if the contact force magnitude is equal to
        // body weight.
        const double arrowLengthPerForce = 1.0 / 1000.0; // meters / Newton
        //const double mg =
        //        getModel().getTotalMass(s) * getModel().getGravity().norm();
        // TODO avoid recalculating.
        const auto& pt = getConnectee<Station>("station");
        const auto pt1 = pt.getLocationInGround(s);
        const SimTK::Vec3& force = getContactForceOnStation(s);
        const SimTK::Vec3 pt2 = pt1 + force * arrowLengthPerForce; // mg;
        SimTK::DecorativeLine line(pt1, pt2);
        line.setColor(SimTK::Green);
        line.setLineThickness(0.10);
        geoms.push_back(line);

        // TODO move to fixed.
        SimTK::DecorativeSphere sphere;
        sphere.setColor(SimTK::Green);
        sphere.setRadius(0.01);
        sphere.setBodyId(pt.getParentFrame().getMobilizedBodyIndex());
        sphere.setRepresentation(SimTK::DecorativeGeometry::DrawWireframe);
        sphere.setTransform(SimTK::Transform(pt.get_location()));
        geoms.push_back(sphere);
    }
}

void StationPlaneContactForce::extendAddToSystem(
        SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);

    this->_forceOnStationCV = addCacheVariable<SimTK::Vec3>(
            "force_on_station", SimTK::Vec3(0),
            SimTK::Stage::Velocity);
}

MeyerFregly2016Force::MeyerFregly2016Force() {
    constructProperties();
}

SimTK::Vec3 MeyerFregly2016Force::calcContactForceOnStation(
        const SimTK::State& state) const {
    SimTK::Vec3 force(0);
    const auto& station = getConnectee<Station>("station");
    const auto& pos = station.getLocationInGround(state);
    const auto& vel = station.getVelocityInGround(state);

    // Limit maximum height used in force calculation.
    SimTK::Real y = pos[1] - get_spring_resting_length();
    // TODO smoothly transition to this value.
    // https://github.com/rcnl-org/nmsm-core/blob/e88992fcb22bad30b589ff46a1af5d069a2c3831/src/GroundContactPersonalization/Optimizations/ModelCalculation/calcModeledVerticalGroundReactionForce.m#L52
    if (y > 0.354237930036971) {
        y = 0.354237930036971;
    }

    // Stiffness constants.
    const SimTK::Real k = get_stiffness();
    const SimTK::Real v = (k + klow) / (k - klow);
    const SimTK::Real s = (k - klow) / 2.0;
    const SimTK::Real constant =
            -s * (v * ymax - c * log(cosh((ymax + h) / c)));

    // Vertical spring-damper force.
    const SimTK::Real Fspring =
            -s * (v * y    - c * log(cosh((y    + h) / c))) - constant;
    force[1] = Fspring * (1 - get_dissipation() * vel[1]);

    // Friction force.
    const SimTK::Real slidingVelocity = 
            sqrt(vel[0] * vel[0] + vel[2] * vel[2]);
    const SimTK::Real horizontalForce = force[1] * (
            get_dynamic_friction() * 
            tanh(slidingVelocity /  get_latch_velocity()) + 
            get_viscous_friction() * slidingVelocity
    );
    
    force[0] = -vel[0] / (slidingVelocity + slipOffset) * horizontalForce;
    force[2] = -vel[2] / (slidingVelocity + slipOffset) * horizontalForce;

    return force;
}

void MeyerFregly2016Force::constructProperties() {
    constructProperty_stiffness(1.0e4);
    constructProperty_dissipation(1.0e-2);
    constructProperty_spring_resting_length(0);
    constructProperty_dynamic_friction(0);
    constructProperty_viscous_friction(5.0);
    constructProperty_latch_velocity(0.05);
}
