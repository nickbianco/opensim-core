/* -------------------------------------------------------------------------- *
 *                     OpenSim:  CoordinateLinearStop.cpp                     *
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

#include "CoordinateLinearStop.h"

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/CoordinateSet.h>

using namespace OpenSim;

CoordinateLinearStop::CoordinateLinearStop() {
    constructProperties();
}

CoordinateLinearStop::CoordinateLinearStop(const std::string& coordinateName,
        double stiffness, double damping, double q_low, double q_high) {
    constructProperties();
    set_coordinate(coordinateName);
    set_stiffness(stiffness);
    set_damping(damping);
    set_q_low(q_low);
    set_q_high(q_high);
}

void CoordinateLinearStop::constructProperties() {
    constructProperty_coordinate("");
    constructProperty_stiffness(0.0);
    constructProperty_damping(0.0);
    constructProperty_q_high(0.0);
    constructProperty_q_low(0.0);
}

void CoordinateLinearStop::extendAddToSystem(
        SimTK::MultibodySystem& system) const {
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

    CoordinateLinearStop* mutableThis =
        const_cast<CoordinateLinearStop*>(this);
    mutableThis->_index = stop.getForceIndex();
}
