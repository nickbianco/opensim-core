/* -------------------------------------------------------------------------- *
 *                     OpenSim:  CoordinateLinearSpring.cpp                   *
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

#include "CoordinateLinearSpring.h"

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/CoordinateSet.h>

using namespace OpenSim;

CoordinateLinearSpring::CoordinateLinearSpring() {
    constructProperties();
}

CoordinateLinearSpring::CoordinateLinearSpring(
        const std::string& coordinateName, double default_stiffness) {
    constructProperties();
    set_coordinate(coordinateName);
    set_default_stiffness(default_stiffness);
}

CoordinateLinearSpring::CoordinateLinearSpring(
        const std::string& coordinateName, double default_stiffness,
        double rest_length) {
    constructProperties();
    set_coordinate(coordinateName);
    set_default_stiffness(default_stiffness);
    set_rest_length(rest_length);
}

void CoordinateLinearSpring::constructProperties() {
    constructProperty_coordinate("");
    constructProperty_default_stiffness(0.0);
    constructProperty_rest_length(0.0);
}

void CoordinateLinearSpring::extendAddToSystem(
        SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);

    const std::string& coordName = get_coordinate();
    const Coordinate& coord = getModel().getCoordinateSet().get(coordName);
    SimTK::MobilizerQIndex qIndex = coord.getMobilizerQIndex();
    SimTK::MobilizedBodyIndex mobodIndex = coord.getBodyIndex();
    const SimTK::MobilizedBody& mobod =
            system.getMatterSubsystem().getMobilizedBody(mobodIndex);

    SimTK::GeneralForceSubsystem& forces = _model->updForceSubsystem();
    SimTK::Force::MobilityLinearSpring spring(forces, mobod,
        qIndex, get_default_stiffness(), get_rest_length());

    CoordinateLinearSpring* mutableThis =
        const_cast<CoordinateLinearSpring*>(this);
    mutableThis->_index = spring.getForceIndex();
}
