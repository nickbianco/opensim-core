/* -------------------------------------------------------------------------- *
 *                        OpenSim: ForceAggregator.cpp                          *
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

#include "ForceAggregator.h"
#include "Model.h"

using namespace OpenSim;

//=============================================================================
// METHODS
//=============================================================================
void ForceAggregator::addForce(const Force& force) {
    OPENSIM_THROW_IF_FRMOBJ(!force.getForceIndex().isValid(), Exception,
        "Force '" + force.getName() + "' has an invalid index.");
    _forceIndexes.push_back(force.getForceIndex());
}

const SimTK::Vector& ForceAggregator::getGeneralizedForces(
        const SimTK::State& state) const {
    computeForces(state);
    return getCacheVariableValue(state, _generalizedForcesCV);
}

SimTK::Real ForceAggregator::getGeneralizedForcesSum(
        const SimTK::State& state) const {
    return getGeneralizedForces(state).sum();
}

const SimTK::Vector_<SimTK::SpatialVec>& ForceAggregator::getBodyForces(
        const SimTK::State& state) const {
    computeForces(state);
    return getCacheVariableValue(state, _bodyForcesCV);
}

SimTK::SpatialVec ForceAggregator::getBodyForcesSum(
        const SimTK::State& state) const {
    return getBodyForces(state).sum();
}

void ForceAggregator::computeForces(const SimTK::State& state) const {
    if (isCacheVariableValid(state, _generalizedForcesCV) &&
            isCacheVariableValid(state, _bodyForcesCV)) {
        return;
    }

    const Model& model = getModel();
    SimTK::Vector generalizedForces;
    SimTK::Vector_<SimTK::SpatialVec> bodyForces;
    model.calcForceContributionsSum(state, _forceIndexes, bodyForces,
            generalizedForces);

    setCacheVariableValue(state, _generalizedForcesCV, generalizedForces);
    setCacheVariableValue(state, _bodyForcesCV, bodyForces);
}

//=============================================================================
// MODEL COMPONENT INTERFACE
//=============================================================================
void ForceAggregator::extendAddToSystem(
        SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);

    SimTK::Vector initGeneralizedForces(
            system.getMatterSubsystem().getNumMobilities(), 0.0);
    SimTK::Vector_<SimTK::SpatialVec> initBodyForces(
            system.getMatterSubsystem().getNumBodies(), SimTK::SpatialVec(0));
    _generalizedForcesCV = addCacheVariable<SimTK::Vector>("generalized_forces",
            initGeneralizedForces, SimTK::Stage::Dynamics);
    _bodyForcesCV = addCacheVariable<SimTK::Vector_<SimTK::SpatialVec>>(
            "body_forces", initBodyForces, SimTK::Stage::Dynamics);
}
