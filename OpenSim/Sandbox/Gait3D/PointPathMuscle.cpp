/* -------------------------------------------------------------------------- *
 *                        OpenSim:  PointPathMuscle.cpp                       *
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

#include "PointPathMuscle.h"
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PhysicalFrame.h>
#include <OpenSim/Simulation/Model/Station.h>
#include <OpenSim/Simulation/Model/ForceConsumer.h>
#include <SimTKsimbody.h>

using namespace OpenSim;

PointPathMuscle::PointPathMuscle(const std::string& name,
                                 const SimTK::Real& maxIsometricForce,
                                 const SimTK::Real& optimalFiberLength,
                                 const SimTK::Real& tendonSlackLength,
                                 const SimTK::Real& pennationAngleAtOptimal) :
    m_maxIsometricForce(maxIsometricForce),
    m_optimalFiberLength(optimalFiberLength),
    m_tendonSlackLength(tendonSlackLength),
    m_pennationAngleAtOptimal(pennationAngleAtOptimal) {
    setName(name);
    constructProperties();
}

void PointPathMuscle::addPoint(const std::string& name, 
                               const PhysicalFrame& frame, const SimTK::Vec3& location) {
    Station* station = new Station(frame, location);
    station->setName(name);
    addComponent(station);
    appendSocketConnectee_points(*station);
}

double PointPathMuscle::getActivation(const SimTK::State& s) const {
    if (get_ignore_activation_dynamics()) {
        return get_control_value();
    } else {
        return getStateVariableValue(s, "activation");
    }
}

SimTK::Real PointPathMuscle::getFiberWidth() const {
    const SimTK::Real normFiberWidth = std::sin(m_pennationAngleAtOptimal);
    return m_optimalFiberLength * normFiberWidth;
}

SimTK::Real PointPathMuscle::getSquareFiberWidth() const {
    return SimTK::square(getFiberWidth());
}

SimTK::Real PointPathMuscle::calcTendonForce(const SimTK::Real& tendonLength) const {
    const SimTK::Real tendonStrain = tendonLength > m_tendonSlackLength ? 
            (tendonLength - m_tendonSlackLength) / m_tendonSlackLength : 0;    
    return m_maxIsometricForce * tendonStrain * (cT1*tendonStrain + cT2);
}

SimTK::Real PointPathMuscle::calcActiveForceLengthMultiplier(
        const SimTK::Real& normalizedFiberLength) const {
    if (r1 < normalizedFiberLength && normalizedFiberLength < r2) {
        SimTK::Real fiberStrain = normalizedFiberLength - 1.0;
        SimTK::Real fiberStrainSquared = fiberStrain*fiberStrain;
        SimTK::Real fiberStrainCubed = fiberStrainSquared*fiberStrain;
        return cL1*fiberStrainCubed + cL2*fiberStrainSquared + 1.0;
    } else {
        return 0;       
    }
}

SimTK::Real PointPathMuscle::calcActiveForceVelocityMultiplier(
        const SimTK::Real& normalizedFiberVelocity) const {
    if (normalizedFiberVelocity >= 0) {
        return (Fvmax*normalizedFiberVelocity + cV2) / 
                (cV2 + normalizedFiberVelocity);
    } else if (-1.0 < normalizedFiberVelocity && normalizedFiberVelocity < 0.0) {
        return cV1*(normalizedFiberVelocity + 1.0) / 
                (cV1 - normalizedFiberVelocity);
    } else {
        return 0;
    }
}

SimTK::Real PointPathMuscle::calcPassiveForceLengthMultiplier(
        const SimTK::Real& normalizedFiberLength) const {
    if (normalizedFiberLength > 1.0) {
        SimTK::Real fiberStrain = normalizedFiberLength - 1.0;
        SimTK::Real fiberStrainSquared = fiberStrain*fiberStrain;
        SimTK::Real fiberStrainCubed = fiberStrainSquared*fiberStrain;
        return cP1*fiberStrainCubed + cP2*fiberStrainSquared;
    } else {
        return 0;       
    }
}

SimTK::Real PointPathMuscle::calcFiberForce(const SimTK::Real& activation,
        const SimTK::Real& normalizedFiberLength,
        const SimTK::Real& normalizedFiberVelocity) const {
    const SimTK::Real fl = calcActiveForceLengthMultiplier(normalizedFiberLength);
    const SimTK::Real fv = calcActiveForceVelocityMultiplier(normalizedFiberVelocity);
    const SimTK::Real fp = calcPassiveForceLengthMultiplier(normalizedFiberLength);
    const SimTK::Real fd = m_fiberDamping * normalizedFiberVelocity;

    return m_maxIsometricForce * (activation*fl*fv + fp + fd);
}

void PointPathMuscle::implProduceForces(const SimTK::State& state,
        ForceConsumer& forceConsumer) const {
    SimTK::Real length = 0;
    SimTK::Real lengthDot = 0;

    const auto& points = getSocket<Station>("points");
    int numPoints = points.getNumConnectees();

    SimTK::Array_<SimTK::UnitVec3> dirs(numPoints-1);
    SimTK::Array_<SimTK::Vec3> s_G(numPoints);

    const auto& point0 = points.getConnectee(0);
    const SimTK::MobilizedBody& body0 = point0.getParentFrame().getMobilizedBody();
    const SimTK::Transform& X_GB0 = body0.getBodyTransform(state);
    const SimTK::Vec3& p0 = point0.get_location();
    SimTK::Vec3 s0_G = X_GB0.R() * p0;
    s_G[0] = s0_G;
    SimTK::Vec3 p0_G = X_GB0.p() + s0_G;
    SimTK::Vec3 v0_G = body0.findStationVelocityInGround(state, p0);
    for (int i = 1; i < numPoints; ++i) {
        const auto& pointi = points.getConnectee(i);
        const SimTK::MobilizedBody& bodyi = pointi.getParentFrame().getMobilizedBody();
        const SimTK::Transform& X_GBi = bodyi.getBodyTransform(state);
        const SimTK::Vec3& pi = pointi.get_location();
        SimTK::Vec3 si_G = X_GBi.R() * pi;
        s_G[i] = si_G;
        SimTK::Vec3 pi_G = X_GBi.p() + si_G;
        SimTK::Vec3 vi_G = bodyi.findStationVelocityInGround(state, pi);
        const SimTK::Vec3& ri_G = pi_G - p0_G; // vector from point0 to pointi
        const SimTK::UnitVec3 dir(ri_G);
        dirs[i-1] = dir; 

        const SimTK::Real dist = ri_G.norm();  // distance between the points
        if( dist < SimTK::SignificantReal ) return;
        length += dist;
        lengthDot += dot(vi_G - v0_G, dir); // relative velocity    

        s0_G = si_G;
        p0_G = pi_G;
        v0_G = vi_G;
    }

    // Rigid tendon. 
    const SimTK::Real fiberLengthAlongTendon = length - m_tendonSlackLength;
    const SimTK::Real fiberLength = sqrt(
        SimTK::square(fiberLengthAlongTendon) + getSquareFiberWidth());
    const SimTK::Real cosPennationAngle = 
        fiberLengthAlongTendon / fiberLength;
    const SimTK::Real normalizedFiberLength = fiberLength / m_optimalFiberLength;
    const SimTK::Real normalizedFiberVelocity = 
        lengthDot / (m_vmax * m_optimalFiberLength);

    const SimTK::Real activation = getActivation(state);
    const SimTK::Real fiberForce = calcFiberForce(activation, 
        normalizedFiberLength, normalizedFiberVelocity);
    // rigid tendon
    const SimTK::Real tendonForce = fiberForce * cosPennationAngle;

    for (int i = 1; i < numPoints; ++i) {
        const SimTK::Vec3 f_G = tendonForce * dirs[i-1];
        SimTK::SpatialVec spatialVecLeft = SimTK::SpatialVec(s_G[i-1] % f_G, f_G);
        SimTK::SpatialVec spatialVecRight = -SimTK::SpatialVec(s_G[i] % f_G, f_G);

        forceConsumer.consumeBodySpatialVec(state, 
                points.getConnectee(i-1).getParentFrame(), 
                spatialVecLeft);
        forceConsumer.consumeBodySpatialVec(state, 
                points.getConnectee(i).getParentFrame(), 
                spatialVecRight);
    }
}

void PointPathMuscle::extendAddToSystem(SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);
    if (!get_ignore_activation_dynamics()) {
        addStateVariable("activation", SimTK::Stage::Dynamics);
    }
}

void PointPathMuscle::computeStateVariableDerivatives(const SimTK::State& s) const {
    if (!get_ignore_activation_dynamics()) {
        const auto& activation = getActivation(s);  
        const auto& excitation = get_control_value();
        static const double c2 = m_deactivationRate;
        static const double c1 = m_activationRate - c2;
        const SimTK::Real derivative = 
                (excitation - activation) * (c1*excitation + c2);
        setStateVariableDerivativeValue(s, "activation", derivative); 
    }
}

void PointPathMuscle::generateDecorations(bool fixed, 
                                         const ModelDisplayHints& hints,
                                         const SimTK::State& state,
                                         SimTK::Array_<SimTK::DecorativeGeometry>& 
                                         appendToThis) const {
    // Don't generate fixed decorations
    if (fixed) return;

    const auto& points = getSocket<Station>("points");
    int numPoints = points.getNumConnectees();

    if (numPoints < 2) return;

    // Get the first point's location in ground
    const auto& point0 = points.getConnectee(0);
    const SimTK::MobilizedBody& body0 = 
        point0.getParentFrame().getMobilizedBody();
    const SimTK::Transform& X_GB0 = body0.getBodyTransform(state);
    const SimTK::Vec3& p0 = point0.get_location();
    SimTK::Vec3 p0_G = X_GB0.p() + X_GB0.R() * p0;

    // Draw line segments between consecutive points
    for (int i = 1; i < numPoints; ++i) {
        const auto& pointi = points.getConnectee(i);
        const SimTK::MobilizedBody& bodyi = 
            pointi.getParentFrame().getMobilizedBody();
        const SimTK::Transform& X_GBi = bodyi.getBodyTransform(state);
        const SimTK::Vec3& pi = pointi.get_location();
        SimTK::Vec3 pi_G = X_GBi.p() + X_GBi.R() * pi;

        // Add line segment
        appendToThis.push_back(SimTK::DecorativeLine(p0_G, pi_G)
            .setLineThickness(4)
            .setColor(SimTK::Red)
            .setBodyId(0)
            .setIndexOnBody(i));

        // Update for next iteration
        p0_G = pi_G;
    }

    // Add spheres at each point for better visibility
    for (int i = 0; i < numPoints; ++i) {
        const auto& pointi = points.getConnectee(i);
        const SimTK::MobilizedBody& bodyi = 
            pointi.getParentFrame().getMobilizedBody();
        const SimTK::Transform& X_GBi = bodyi.getBodyTransform(state);
        const SimTK::Vec3& pi = pointi.get_location();
        SimTK::Vec3 pi_G = X_GBi.p() + X_GBi.R() * pi;

        // Add sphere at each point
        appendToThis.push_back(SimTK::DecorativeSphere(0.005)
            .setTransform(SimTK::Transform(pi_G))
            .setColor(SimTK::Blue)
            .setBodyId(0)
            .setIndexOnBody(numPoints + i));
    }
}

// Control value methods
void PointPathMuscle::setControlValue(double controlValue) {
    set_control_value(controlValue);
}

double PointPathMuscle::getControlValue() const {
    return get_control_value();
}

void PointPathMuscle::constructProperties() {
    constructProperty_control_value(0.1);
    constructProperty_ignore_activation_dynamics(false);
} 