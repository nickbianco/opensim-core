/* -------------------------------------------------------------------------- *
 *                     OpenSim:  CoordinateLimitForce.cpp                     *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Ajay Seth                                                       *
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

#include "ExponentialCoordinateForce.h"

#include <OpenSim/Simulation/Model/ForceConsumer.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>

using namespace OpenSim;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
ExponentialCoordinateForce::ExponentialCoordinateForce() {
    constructProperties();
}

void ExponentialCoordinateForce::constructProperties() {
    constructProperty_coordinate("");
    constructProperty_k0(0.0);
    constructProperty_k1(0.0);
    constructProperty_k2(0.0);
    constructProperty_k3(0.0);
    constructProperty_k4(0.0);
    constructProperty_theta(0.0);
    constructProperty_phi(0.0);
}

//=============================================================================
// METHODS
//=============================================================================
void ExponentialCoordinateForce::extendConnectToModel(Model& aModel) {
    Super::extendConnectToModel(aModel);

    const auto& coordinateNameOrPath = get_coordinate();
    if (_model->getCoordinateSet().contains(coordinateNameOrPath)) {
        _coord = &_model->getCoordinateSet().get(coordinateNameOrPath);
    } else if (_model->hasComponent<Coordinate>(coordinateNameOrPath)) {
        _coord = &_model->getComponent<Coordinate>(coordinateNameOrPath);
    } else {
        OPENSIM_THROW_FRMOBJ(Exception,
            "Received '{}' from property 'coordinate', but no coordinate "
            "with this name or path was found in the model.",
            coordinateNameOrPath);
    }
}

//=============================================================================
// COMPUTATIONS
//=============================================================================
double ExponentialCoordinateForce::calcForce(const SimTK::State& s) const {
    double q = _coord->getValue(s);
    const auto& k0 = get_k0();
    const auto& k1 = get_k1();
    const auto& k2 = get_k2();
    const auto& k3 = get_k3();
    const auto& k4 = get_k4();
    const auto& theta = get_theta();
    const auto& phi = get_phi();

    return k0 + k1*std::exp(k2*(q - theta)) + k3*std::exp(k4*(q - phi));
}

double ExponentialCoordinateForce::computePotentialEnergy(
        const SimTK::State& s) const {

    double q = _coord->getValue(s);
    const auto& k1 = get_k1();
    const auto& k2 = get_k2();
    const auto& k3 = get_k3();
    const auto& k4 = get_k4();
    const auto& theta = get_theta();
    const auto& phi = get_phi();

    // df/dq
    return k1*k2*std::exp(k2*(q - theta)) + k3*k4*std::exp(k4*(q - phi));
}

//=============================================================================
// FORCE PRODUCER INTERFACE
//=============================================================================
void ExponentialCoordinateForce::implProduceForces(const SimTK::State& s,
    ForceConsumer& forceConsumer) const
{
    forceConsumer.consumeGeneralizedForce(s, *_coord, calcForce(s));
}

//=============================================================================
// REPORTING
//=============================================================================
Array<std::string> ExponentialCoordinateForce::getRecordLabels() const {
    OpenSim::Array<std::string> labels("");
    labels.append(getName());
    labels.append("PotentialEnergy");
    return labels;
}

Array<double> ExponentialCoordinateForce::getRecordValues(
        const SimTK::State& state) const {
    OpenSim::Array<double> values(0.0, 0, 2);
    values.append(calcForce(state));
    values.append(computePotentialEnergy(state));
    return values;
}
