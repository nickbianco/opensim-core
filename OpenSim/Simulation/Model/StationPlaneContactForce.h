#ifndef OPENSIM_STATIONPLANECONTACTFORCE_H
#define OPENSIM_STATIONPLANECONTACTFORCE_H
/* -------------------------------------------------------------------------- *
 * OpenSim: StationPlaneContactForce.h                                        *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2024 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Nicholas Bianco, Christopher Dembia, Spencer Williams           *
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

#include <OpenSim/Simulation/Model/ForceConsumer.h>
#include <OpenSim/Simulation/Model/ForceProducer.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Station.h>

namespace OpenSim {

/** This class models compliant point contact with a ground plane y=0.

Vertical force is calculated using stiffness and dissipation coefficients, 
while horizontal force is calculated using a friction model. */
class OSIMSIMULATION_API StationPlaneContactForce : public ForceProducer {
OpenSim_DECLARE_ABSTRACT_OBJECT(StationPlaneContactForce, ForceProducer);

public:
//=============================================================================
// OUTPUTS
//=============================================================================
    OpenSim_DECLARE_OUTPUT(force_on_station, SimTK::Vec3,
            getContactForceOnStation, SimTK::Stage::Velocity);

//=============================================================================
// SOCKETS
//=============================================================================
    OpenSim_DECLARE_SOCKET(station, Station,
            "The body-fixed point that can contact the plane.");

//=============================================================================
// METHODS
//=============================================================================

    // CONSTRUCTION
    StationPlaneContactForce() = default;

    /** TODO */
    SimTK::Vec3 getContactForceOnStation(const SimTK::State& s) const {
        if (!isCacheVariableValid(s, _forceOnStationCV)) {
            setContactForceOnStation(s, calcContactForceOnStation(s));
        }
        return getCacheVariableValue<SimTK::Vec3>(s, _forceOnStationCV);
    }

    /** TODO */
    void setContactForceOnStation(const SimTK::State& s,
            const SimTK::Vec3& force) const {
        setCacheVariableValue(s, _forceOnStationCV, force);
    }

    // FORCE INTERFACE
    OpenSim::Array<std::string> getRecordLabels() const override;
    OpenSim::Array<double> getRecordValues(const SimTK::State& s) const override;

    // COMPONENT INTERFACE
    void generateDecorations(bool fixed, const ModelDisplayHints& hints,
            const SimTK::State& s,
            SimTK::Array_<SimTK::DecorativeGeometry>& geoms) const override;

protected:
    /** TODO */
    virtual SimTK::Vec3 calcContactForceOnStation(
            const SimTK::State& s) const = 0;

    // COMPONENT INTERFACE
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    
private:
    mutable CacheVariable<SimTK::Vec3> _forceOnStationCV;

    // FORCE PRODUCER INTERFACE
    void implProduceForces(const SimTK::State& s, ForceConsumer& forceConsumer) 
            const override {

        const SimTK::Vec3& force = getContactForceOnStation(s);
        const auto& pt = getConnectee<Station>("station");
        const auto& pos = pt.getLocationInGround(s);
        const auto& frame = pt.getParentFrame();

        forceConsumer.consumePointForce(s, frame, pt.get_location(), force);
        forceConsumer.consumePointForce(s, getModel().getGround(), pos, -force);
    }
};

/** This contact model is from the following paper:
Meyer A. J., Eskinazi, I., Jackson, J. N., Rao, A. V., Patten, C., & Fregly,
B. J. (2016). Muscle Synergies Facilitate Computational Prediction of
Subject-Specific Walking Motions. Frontiers in Bioengineering and
Biotechnology, 4, 1055–27. http://doi.org/10.3389/fbioe.2016.00077 

Following OpenSim convention, this contact element assumes that the y-direction
is vertical. Vertical contact force is calculated based on vertical position 
and velocity relative to a floor at y=0 using these equations:

\f[
v = \frac{k_{val} + k_{low}}{k_{val} - k_{low}}
\f]
\f[
s = \frac{k_{val} - k_{low}}{2}
\f]
\f[
R = -s * (v * y_{max} - c * log(\frac{cosh(y_{max} + h)}{c})
\f]
\f[
F = -s * (v * y - c * log(\frac{cosh(y + h)}{c}) - R
\f]

With the following values:
- \f$ k_{val} \f$: stiffness coefficient of contact element
- \f$ k_{low} \f$: small out-of-contact stiffness to assist optimizations
- \f$ y_{max} \f$: y value were out-of-contact force becomes zero
- \f$ c \f$: transition curvature of transition between linear regions
- \f$ h \f$: horizontal offset defining slope transition location
- \f$ y \f$: current y (vertical) position of contact element

Velocity is then used to incorporate non-linear damping:

\f[
F_{damped} = F * (1 - C_{val} * y_{vel})
\f]

With the following values:
- \f$ C_{val} \f$: damping coefficient of contact element
- \f$ y_{vel} \f$: current y (vertical) velocity of contact element

The force equation produces a force-penetration curve similar to a leaky 
rectified linear function with a smooth transition region near zero. This 
produces a smooth curve with a slight out-of-contact slope to assist 
gradient-based optimizations when an element is inappropriately out of contact.

Horizontal forces are then calculated based on this vertical force and the 
horizontal velocity components of the contact element. Both dynamic (modeled 
with a tanh function) and viscous (modeled with a linear function) friction 
models may be used. */
class OSIMSIMULATION_API MeyerFregly2016Force : public StationPlaneContactForce {
OpenSim_DECLARE_CONCRETE_OBJECT(MeyerFregly2016Force, StationPlaneContactForce);

public:
//=============================================================================
// PROPERTIES
//=============================================================================
    OpenSim_DECLARE_PROPERTY(stiffness, double,
            "Spring stiffness in N/m (default: 1e4).");
    OpenSim_DECLARE_PROPERTY(dissipation, double,
            "Dissipation coefficient in s/m (default: 0.01).");
    OpenSim_DECLARE_PROPERTY(spring_resting_length, double,
            "Spring resting length in m (default: 0).");
    OpenSim_DECLARE_PROPERTY(dynamic_friction, double,
            "Dynamic friction coefficient (default: 0).");
    OpenSim_DECLARE_PROPERTY(viscous_friction, double,
            "Viscous friction coefficient (default: 5.0).");
    OpenSim_DECLARE_PROPERTY(latch_velocity, double,
            "Latching velocity in m/s (default: 0.05).");

//=============================================================================
// METHODS
//=============================================================================

    // CONSTRUCTION
    MeyerFregly2016Force();

private:
    // STATION PLANE CONTACT FORCE INTERFACE
    SimTK::Vec3 calcContactForceOnStation(
            const SimTK::State& state) const override; 

    // PROPERTIES
    void constructProperties();

    // CONSTANTS
    constexpr static SimTK::Real klow = 1e-1;
    constexpr static SimTK::Real h = 1e-3;
    constexpr static SimTK::Real c = 5e-4;
    constexpr static SimTK::Real ymax = 1e-2;
    // Prevents dividing by zero when sliding velocity is zero.
    constexpr static SimTK::Real slipOffset = 1e-4;
};

} // namespace OpenSim

#endif // OPENSIM_STATIONPLANECONTACTFORCE_H
