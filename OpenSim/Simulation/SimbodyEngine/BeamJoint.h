#ifndef OPENSIM_BEAM_JOINT_H_
#define OPENSIM_BEAM_JOINT_H_
/* -------------------------------------------------------------------------- *
 *                          OpenSim:  BeamJoint.h                             *
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

#include "Joint.h"

namespace OpenSim {

class Model;

class OSIMSIMULATION_API BeamJoint : public Joint {
OpenSim_DECLARE_CONCRETE_OBJECT(BeamJoint, Joint);

public:
    enum class Coord: unsigned {
        Rotation1X = 0u, ///< 0
        Rotation2Y = 1u, ///< 1
        Rotation3Z = 2u  ///< 2
    };

private:
    /** Specify the Coordinates of the BeamJoint. */
    CoordinateIndex rx{ constructCoordinate(Coordinate::MotionType::Rotational,
                                   static_cast<unsigned>(Coord::Rotation1X)) };
    CoordinateIndex ry{ constructCoordinate(Coordinate::MotionType::Rotational,
                                   static_cast<unsigned>(Coord::Rotation2Y)) };
    CoordinateIndex rz{ constructCoordinate(Coordinate::MotionType::Rotational,
                                   static_cast<unsigned>(Coord::Rotation3Z)) };

public:
//==============================================================================
// PROPERTIES
//==============================================================================
    OpenSim_DECLARE_PROPERTY(beam_length, SimTK::Real,
        "The length of the beam.");

//=============================================================================
// METHODS
//=============================================================================
    BeamJoint();

    BeamJoint(const std::string&    name,
                const PhysicalFrame&  parent,
                const PhysicalFrame&  child,
                const SimTK::Real&    beamLength);

    BeamJoint(const std::string&    name,
              const PhysicalFrame&  parent,
              const SimTK::Vec3&    locationInParent,
              const SimTK::Vec3&    orientationInParent,
              const PhysicalFrame&  child,
              const SimTK::Vec3&    locationInChild,
              const SimTK::Vec3&    orientationInChild,
              const SimTK::Real&    beamLength);

    /** Use Joint's constructors. @see Joint */
    using Joint::Joint;

    /** Exposes getCoordinate() method defined in base class (overloaded below).
        @see Joint */
    using Joint::getCoordinate;

    /** Exposes updCoordinate() method defined in base class (overloaded below).
        @see Joint */
    using Joint::updCoordinate;

    /** Get a const reference to a Coordinate associated with this Joint.
        @see Coord */
    const Coordinate& getCoordinate(Coord idx) const {
        return get_coordinates( static_cast<unsigned>(idx) );
    }

    /** Get a writable reference to a Coordinate associated with this Joint.
        @see Coord */
    Coordinate& updCoordinate(Coord idx) {
        return upd_coordinates( static_cast<unsigned>(idx) );
    }

protected:
    // MODEL COMPONENT INTERFACE
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;

    void generateDecorations(bool fixed, const ModelDisplayHints& hints,
        const SimTK::State& state,
        SimTK::Array_<SimTK::DecorativeGeometry>& geometry)
        const override;

private:
    void constructProperties();
};

} // namespace OpenSim

#endif // OPENSIM_BEAM_JOINT_H_
