/* -------------------------------------------------------------------------- *
 *                         OpenSim:  BeamJoint.cpp                            *
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

#include "BeamJoint.h"
#include <OpenSim/Simulation/Model/Model.h>
#include "simbody/internal/MobilizedBody_Beam.h"

using namespace OpenSim;

BeamJoint::BeamJoint() : Super() {
    constructProperties();
}

BeamJoint::BeamJoint(const std::string&    name,
                     const PhysicalFrame&  parent,
                     const PhysicalFrame&  child,
                     const SimTK::Real&    beamLength) :
                     Super(name, parent, child) {
    constructProperties();
    set_beam_length(beamLength);
}

BeamJoint::BeamJoint(const std::string&    name,
                     const PhysicalFrame&  parent,
                     const SimTK::Vec3&    locationInParent,
                     const SimTK::Vec3&    orientationInParent,
                     const PhysicalFrame&  child,
                     const SimTK::Vec3&    locationInChild,
                     const SimTK::Vec3&    orientationInChild,
                     const SimTK::Real&    beamLength) :
    Super(name, parent, locationInParent, orientationInParent,
          child, locationInChild, orientationInChild) {
    constructProperties();
    set_beam_length(beamLength);
}

void BeamJoint::constructProperties() {
    constructProperty_beam_length(1.0);
}

void BeamJoint::extendAddToSystem(SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);
    SimTK::MobilizedBody::Beam mobod =
        createMobilizedBody<SimTK::MobilizedBody::Beam>(system);
    mobod.setDefaultLength(get_beam_length());
}
