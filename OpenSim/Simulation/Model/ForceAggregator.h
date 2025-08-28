#ifndef OPENSIM_FORCE_AGGREGATOR_H_
#define OPENSIM_FORCE_AGGREGATOR_H_
/* -------------------------------------------------------------------------- *
 *                        OpenSim: ForceAggregator.h                          *
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

#include <OpenSim/Simulation/osimSimulationDLL.h>
#include "Force.h"

namespace OpenSim {

class OSIMSIMULATION_API ForceAggregator : public ModelComponent {
OpenSim_DECLARE_CONCRETE_OBJECT(ForceAggregator, ModelComponent);

public:
//=============================================================================
// OUTPUTS
//=============================================================================
    // OpenSim_DECLARE_OUTPUT(generalized_forces, SimTK::Vector,
    //         getGeneralizedForces, SimTK::Stage::Dynamics);
    // OpenSim_DECLARE_OUTPUT(generalized_forces_sum, SimTK::Real,
    //         getGeneralizedForcesSum, SimTK::Stage::Dynamics);
    // OpenSim_DECLARE_OUTPUT(body_forces, SimTK::Vector_<SimTK::SpatialVec>,
    //         getBodyForces, SimTK::Stage::Dynamics);
    // OpenSim_DECLARE_OUTPUT(body_forces_sum, SimTK::SpatialVec,
    //         getBodyForcesSum, SimTK::Stage::Dynamics);

//=============================================================================
// METHODS
//=============================================================================

    // CONSTRUCTION
    ForceAggregator() = default;

    void addForce(const Force& force);
    const SimTK::Vector& getGeneralizedForces(
            const SimTK::State& state) const;
    SimTK::Real getGeneralizedForcesSum(const SimTK::State& state) const;

    const SimTK::Vector_<SimTK::SpatialVec>& getBodyForces(
            const SimTK::State& state) const;
    SimTK::SpatialVec getBodyForcesSum(const SimTK::State& state) const;

private:
    // MODEL COMPONENT INTERFACE
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;

    // CONVENIENCE METHODS
    void computeForces(const SimTK::State& state) const;

    // MEMBER VARIABLES
    SimTK::Array_<SimTK::ForceIndex> _forceIndexes;

    // CACHE VARIABLES
    mutable CacheVariable<SimTK::Vector> _generalizedForcesCV;
    mutable CacheVariable<SimTK::Vector_<SimTK::SpatialVec>> _bodyForcesCV;
};

} // namespace OpenSim

 #endif // OPENSIM_FORCE_AGGREGATOR_H_