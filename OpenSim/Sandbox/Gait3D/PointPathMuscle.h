/* -------------------------------------------------------------------------- *
 *                        OpenSim:  PointPathMuscle.h                         *
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

#ifndef OPENSIM_POINTPATHMUSCLE_H
#define OPENSIM_POINTPATHMUSCLE_H

#include <OpenSim/Simulation/Model/ForceProducer.h>
#include <OpenSim/Simulation/Model/Station.h>

namespace OpenSim {

class PointPathMuscle : public ForceProducer {
    OpenSim_DECLARE_CONCRETE_OBJECT(PointPathMuscle, ForceProducer);
public:
    OpenSim_DECLARE_LIST_SOCKET(points, Station, "TODO");
    OpenSim_DECLARE_PROPERTY(control_value, double, 
            "Constant control value (excitation) for the muscle (0.0 to 1.0)");
    OpenSim_DECLARE_PROPERTY(ignore_activation_dynamics, bool, 
            "If true, the activation dynamics are ignored and the control value is used as the activation");    

    PointPathMuscle(const std::string& name,
                    const SimTK::Real& maxIsometricForce,
                    const SimTK::Real& optimalFiberLength,
                    const SimTK::Real& tendonSlackLength,
                    const SimTK::Real& pennationAngleAtOptimal);

    void addPoint(const std::string& name, 
                  const PhysicalFrame& frame, const SimTK::Vec3& location);

    double getActivation(const SimTK::State& s) const;
    
    // Control value methods
    void setControlValue(double controlValue);
    double getControlValue() const;

    SimTK::Real getFiberWidth() const;
    SimTK::Real getSquareFiberWidth() const;

    SimTK::Real calcTendonForce(const SimTK::Real& tendonLength) const;
    SimTK::Real calcActiveForceLengthMultiplier(const SimTK::Real& normalizedFiberLength) const;
    SimTK::Real calcActiveForceVelocityMultiplier(const SimTK::Real& normalizedFiberVelocity) const;
    SimTK::Real calcPassiveForceLengthMultiplier(const SimTK::Real& normalizedFiberLength) const;
    SimTK::Real calcFiberForce(const SimTK::Real& activation,
                               const SimTK::Real& normalizedFiberLength,
                               const SimTK::Real& normalizedFiberVelocity) const;

    void implProduceForces(const SimTK::State& state,
                          ForceConsumer& forceConsumer) const override;

private:
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void computeStateVariableDerivatives(const SimTK::State& s) const override;
    void generateDecorations(bool fixed, 
                            const ModelDisplayHints& hints,
                            const SimTK::State& state,
                            SimTK::Array_<SimTK::DecorativeGeometry>& 
                            appendToThis) const override;
    void constructProperties();

    SimTK::Real m_maxIsometricForce;
    SimTK::Real m_optimalFiberLength;
    SimTK::Real m_tendonSlackLength;
    SimTK::Real m_pennationAngleAtOptimal;
    SimTK::Real m_vmax = 10.0;
    SimTK::Real m_fiberDamping = 0.1;
    SimTK::Real m_deactivationRate = 25.0;
    SimTK::Real m_activationRate = 100.0;

    // Tendon force-length curve.
    constexpr static SimTK::Real cT1 = 260.972;
    constexpr static SimTK::Real cT2 = 7.9706;
    // Fiber active force-length curve.
    constexpr static SimTK::Real cL1 = 1.5;
    constexpr static SimTK::Real cL2 = -2.75;
    constexpr static SimTK::Real r1 = 0.46899;
    constexpr static SimTK::Real r2 = 1.80528;
    // Fiber active force-velocity curve.
    constexpr static SimTK::Real cV1 = 0.227;
    constexpr static SimTK::Real cV2 = 0.110;
    constexpr static SimTK::Real Fvmax = 1.6;
    // Fiber passive force-length curve.
    constexpr static SimTK::Real cP1 = 1.08027;
    constexpr static SimTK::Real cP2 = 1.27368;
    // Activation dynamics.
    constexpr static SimTK::Real c1 = 75;
    constexpr static SimTK::Real c2 = 25;
};

} // namespace OpenSim

#endif // OPENSIM_POINTPATHMUSCLE_H 