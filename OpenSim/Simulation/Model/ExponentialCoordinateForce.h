#ifndef OPENSIM_EXPONENTIAL_COORDINATE_FORCE_H_
#define OPENSIM_EXPONENTIAL_COORDINATE_FORCE_H_
/* -------------------------------------------------------------------------- *
 *                    OpenSim:  ExponentialCoordinateForce.h                  *
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

#include <OpenSim/Simulation/Model/ForceProducer.h>

namespace OpenSim {

class OSIMSIMULATION_API ExponentialCoordinateForce : public ForceProducer {
OpenSim_DECLARE_CONCRETE_OBJECT(ExponentialCoordinateForce, ForceProducer);
public:
//==============================================================================
// PROPERTIES
//==============================================================================
    OpenSim_DECLARE_PROPERTY(coordinate, std::string,
        "The name or full path of the coordinate this force is applied to.");
    OpenSim_DECLARE_PROPERTY(k0, double, "TODO.");
    OpenSim_DECLARE_PROPERTY(k1, double, "TODO.");
    OpenSim_DECLARE_PROPERTY(k2, double, "TODO.");
    OpenSim_DECLARE_PROPERTY(k3, double, "TODO.");
    OpenSim_DECLARE_PROPERTY(k4, double, "TODO.");
    OpenSim_DECLARE_PROPERTY(theta, double, "TODO.");
    OpenSim_DECLARE_PROPERTY(phi, double, "TODO.");

//=============================================================================
// METHODS
//=============================================================================
    /** Default constructor */
    ExponentialCoordinateForce();

    // COMPUTATIONS
    double calcForce(const SimTK::State& s) const;
    double computePotentialEnergy(const SimTK::State& s) const override;

    // REPORTING
    Array<std::string> getRecordLabels() const override;
    Array<double> getRecordValues(const SimTK::State& state) const override;

protected:
    // MODEL COMPONENT INTERFACE
    void extendConnectToModel(Model& aModel) override;

private:
    // FORCE PRODUCER INTERFACE
    void implProduceForces(const SimTK::State&, ForceConsumer&) const override;

    // HELPERS
    void constructProperties();

    SimTK::ReferencePtr<const Coordinate> _coord;

};  // class ExponentialCoordinateForce

} // namespace OpenSim

#endif // #ifndef OPENSIM_EXPONENTIAL_COORDINATE_FORCE_H_
