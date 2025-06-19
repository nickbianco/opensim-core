#ifndef HYFYDY_MUSCLE_H
#define HYFYDY_MUSCLE_H
/* -------------------------------------------------------------------------- *
 *                        OpenSim:  HyfydyMuscle.h                            *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2025 Stanford University and the Authors                *
 * Author(s): Thomas Geijtenbeek                                              *
 * Contributors(s): Nicholas Bianco                                           *
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

#include <OpenSim/Actuators/osimActuatorsDLL.h>

#include <OpenSim/Common/DataTable.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>

namespace OpenSim {

/** This muscle model is based on the muscle model from the Hyfydy dynamics 
engine. */
class OSIMACTUATORS_API HyfydyMuscle : public Muscle {
    OpenSim_DECLARE_CONCRETE_OBJECT(HyfydyMuscle, Muscle);

public:
    OpenSim_DECLARE_PROPERTY(activation_rate, double,
            "Larger value means activation can increase more rapidly. "
            "Default: 100 [1/seconds]. Bounds: (0, inf]");
    OpenSim_DECLARE_PROPERTY(deactivation_rate, double,
            "Larger value means activation can decrease more rapidly. "
            "Default: 25 [1/seconds]. Bounds: (0, inf]");
    OpenSim_DECLARE_PROPERTY(default_activation, double,
            "Value of activation in the default state returned by "
            "initSystem(). Default: 0.5. Bounds: (0, inf]");
    OpenSim_DECLARE_PROPERTY(fiber_damping, double,
            "Use this property to define the linear damping force that is "
            "added to the total muscle fiber force. It is computed by "
            "multiplying this damping parameter by the normalized fiber "
            "velocity and the max isometric force. Default: 0. Bounds: [0, inf]");
    OpenSim_DECLARE_PROPERTY(ignore_passive_fiber_force, bool,
            "Make the passive fiber force 0. Default: false.");
    OpenSim_DECLARE_OUTPUT(passive_fiber_elastic_force, double,
            getPassiveFiberElasticForce, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(passive_fiber_elastic_force_along_tendon, double,
            getPassiveFiberElasticForceAlongTendon, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(passive_fiber_damping_force, double,
            getPassiveFiberDampingForce, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(passive_fiber_damping_force_along_tendon, double,
            getPassiveFiberDampingForceAlongTendon, SimTK::Stage::Dynamics);

    HyfydyMuscle() { constructProperties(); }


protected:
    //--------------------------------------------------------------------------
    // COMPONENT INTERFACE
    //--------------------------------------------------------------------------
    /// @name Component interface
    /// @{
    void extendFinalizeFromProperties() override;
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void extendInitStateFromProperties(SimTK::State& s) const override;
    void extendSetPropertiesFromState(const SimTK::State& s) override;
    void computeStateVariableDerivatives(const SimTK::State& s) const override;
    /// @}

    //--------------------------------------------------------------------------
    // ACTUATOR INTERFACE
    //--------------------------------------------------------------------------
    /// @name Actuator interface
    /// @{
    double computeActuation(const SimTK::State& s) const override;
    /// @}

public:
    //--------------------------------------------------------------------------
    // MUSCLE INTERFACE
    //--------------------------------------------------------------------------
    /// @name Muscle interface
    /// @{

    /// If ignore_activation_dynamics is true, this gets excitation instead.
    double getActivation(const SimTK::State& s) const override {
        // We override the Muscle's implementation because Muscle requires
        // realizing to Dynamics to access activation from MuscleDynamicsInfo,
        // which is unnecessary if the activation is a state.
        if (get_ignore_activation_dynamics()) {
            return getControl(s);
        } else {
            return getStateVariableValue(s, STATE_ACTIVATION_NAME);
        }
    }

    /// If ignore_activation_dynamics is true, this sets excitation instead.
    void setActivation(SimTK::State& s, double activation) const override {
        if (get_ignore_activation_dynamics()) {
            SimTK::Vector& controls(getModel().updControls(s));
            setControls(SimTK::Vector(1, activation), controls);
            getModel().setControls(s, controls);
        } else {
            setStateVariableValue(s, STATE_ACTIVATION_NAME, activation);
        }
        markCacheVariableInvalid(s, "velInfo");
        markCacheVariableInvalid(s, "dynamicsInfo");
    }

protected:
    double calcInextensibleTendonActiveFiberForce(
        SimTK::State&, double) const override;
    void calcMuscleLengthInfo(
            const SimTK::State& s, MuscleLengthInfo& mli) const override;
    void calcFiberVelocityInfo(
            const SimTK::State& s, FiberVelocityInfo& fvi) const override;
    void calcMuscleDynamicsInfo(
            const SimTK::State& s, MuscleDynamicsInfo& mdi) const override;
    void calcMusclePotentialEnergyInfo(const SimTK::State& s,
            MusclePotentialEnergyInfo& mpei) const override;

public:
    /// In this method, calcEquilibriumResidual() is used to find a value of the
    /// normalized tendon force state variable that produces muscle-tendon
    /// equilibrium. This relies on the implicit form of tendon compliance since
    /// the explicit form uses the normalized tendon force state variable
    /// directly to compute fiber force, which always produces a zero
    /// muscle-tendon equilibrium residual. The derivative of normalized tendon
    /// force is set to zero since a value is required for the implicit form of
    /// the model.  
    void computeInitialFiberEquilibrium(SimTK::State& s) const override;
    /// @}

    /// @name Get methods.
    /// @{

    /// Get the portion of the passive fiber force generated by the elastic
    /// element only (N).
    double getPassiveFiberElasticForce(const SimTK::State& s) const;
    /// Get the portion of the passive fiber force generated by the elastic
    /// element only, projected onto the tendon direction (N).
    double getPassiveFiberElasticForceAlongTendon(const SimTK::State& s) const;
    /// Get the portion of the passive fiber force generated by the damping
    /// element only (N).
    double getPassiveFiberDampingForce(const SimTK::State& s) const;
    /// Get the portion of the passive fiber force generated by the damping
    /// element only, projected onto the tendon direction (N).
    double getPassiveFiberDampingForceAlongTendon(const SimTK::State& s) const;

    static std::string getActivationStateName() {
        return STATE_ACTIVATION_NAME;
    }

    static double getMinNormalizedFiberLength() { return m_minNormFiberLength; }
    static double getMaxNormalizedFiberLength() { return m_maxNormFiberLength; }
    /// The first element of the Vec2 is the lower bound, and the second is the
    /// upper bound.
    /// Note that since fiber length is not used as a state variable, these
    /// bounds cannot be enforced directly. It is upon the user to ensure the
    /// muscle fiber is operating within the specified domain.
    SimTK::Vec2 getBoundsNormalizedFiberLength() const {
        return {getMinNormalizedFiberLength(), getMaxNormalizedFiberLength()};
    }
    /// @}

    /// @name Calculation methods.
    /// These functions compute the values of normalized/dimensionless curves
    /// their derivatives, and other quantities of the muscle.
    /// These do not depend on a SimTK::State.
    /// @{

    /// The active force-length curve is a polynomial function of the fiber 
    /// strain.
    /// Domain: [0, 1]
    /// Range: [0.46899, 1.80528]
    static SimTK::Real calcActiveForceLengthMultiplier(
            const SimTK::Real& normFiberLength) {
        if (normFiberLength > m_minNormFiberLength && 
                normFiberLength < m_maxNormFiberLength) {
            const SimTK::Real fiberStrain = normFiberLength - 1.0;
            const SimTK::Real fiberStrainSquared = fiberStrain * fiberStrain;
            const SimTK::Real fiberStrainCubed = 
                    fiberStrainSquared * fiberStrain;
            return cL1 * fiberStrainCubed + cL2 * fiberStrainSquared + 1.0;
        } else {
            return 0.0;
        } 
    }

    /// The derivative of the active force-length curve with respect to
    /// normalized fiber length.
    static SimTK::Real calcActiveForceLengthMultiplierDerivative(
            const SimTK::Real& normFiberLength) {
        if (normFiberLength > m_minNormFiberLength && 
                normFiberLength < m_maxNormFiberLength) {
            const SimTK::Real fiberStrain = normFiberLength - 1.0;
            const SimTK::Real fiberStrainSquared = fiberStrain * fiberStrain;
            return 3.0 * cL1 * fiberStrainSquared + 2.0 * cL2 * fiberStrain;
        } else {
            return 0.0;
        } 
    }


    /// The force velocity multiplier is a polynomial function of the normalized
    /// fiber velocity.
    /// Domain: [-1, 1]
    /// Range: [0, 1.6]
    static SimTK::Real calcForceVelocityMultiplier(
            const SimTK::Real& normFiberVelocity) {
        if (normFiberVelocity <= -1) {
            return 0.0;
        } else if (normFiberVelocity >= 0) {
            return (Fvmax*normFiberVelocity + cV2) / (cV2 + normFiberVelocity);
        } else {
            return (cV1*(normFiberVelocity + 1.0)) / (cV1 - normFiberVelocity);
        }
    }

    /// The passive force-length curve is a polynomial function of the
    /// normalized fiber length.
    SimTK::Real calcPassiveForceMultiplier(
            const SimTK::Real& normFiberLength) const {
        if (get_ignore_passive_fiber_force()) return 0;
        if (normFiberLength > 1.0) {
            const SimTK::Real fiberStrain = normFiberLength - 1.0;
            const SimTK::Real fiberStrainSquared = fiberStrain * fiberStrain;
            const SimTK::Real fiberStrainCubed = 
                    fiberStrainSquared * fiberStrain;
            return cP1 * fiberStrainCubed + cP2 * fiberStrainSquared;
        } else {
            return 0.0;
        }
    }

    /// This is the derivative of the passive force-length curve with respect to
    /// the normalized fiber length.
    SimTK::Real calcPassiveForceMultiplierDerivative(
            const SimTK::Real& normFiberLength) const {
        if (get_ignore_passive_fiber_force()) return 0;
        if (normFiberLength > 1.0) {
            const SimTK::Real fiberStrain = normFiberLength - 1.0;
            const SimTK::Real fiberStrainSquared = fiberStrain * fiberStrain;
            return 3.0 * cP1 * fiberStrainSquared + 2.0 * cP2 * fiberStrain;
        } else {
            return 0.0;
        }
    }

    /// The normalized tendon force as a polynomiaal function of normalized 
    /// tendon length.
    SimTK::Real calcTendonForceMultiplier(
            const SimTK::Real& normTendonLength) const {
        if (normTendonLength > 1.0) {
            const SimTK::Real tendonStrain = normTendonLength - 1.0;
            const SimTK::Real tendonStrainSquared = tendonStrain * tendonStrain;
            return cT1 * tendonStrainSquared + cT2 * tendonStrain;
        } else {
            return 0.0;
        }
    }

    /// This is the derivative of the tendon-force length curve with respect to
    /// normalized tendon length.
    SimTK::Real calcTendonForceMultiplierDerivative(
            const SimTK::Real& normTendonLength) const {
        if (normTendonLength > 1.0) {
            const SimTK::Real tendonStrain = normTendonLength - 1.0;
            return 2.0 * cT1 * tendonStrain + cT2;
        } else {
            return 0.0;
        }
    }

    /// This computes both the total fiber force and the individual components
    /// of fiber force (active, conservative passive, and non-conservative
    /// passive).
    /// @note based on Millard2012EquilibriumMuscle::calcFiberForce().
    void calcFiberForce(const SimTK::Real& activation,
            const SimTK::Real& activeForceLengthMultiplier,
            const SimTK::Real& forceVelocityMultiplier,
            const SimTK::Real& normPassiveFiberForce,
            const SimTK::Real& normFiberVelocity, 
            SimTK::Real& activeFiberForce,
            SimTK::Real& conPassiveFiberForce,
            SimTK::Real& nonConPassiveFiberForce,
            SimTK::Real& totalFiberForce) const {
        const auto& maxIsometricForce = get_max_isometric_force();
        // active force
        activeFiberForce =
                maxIsometricForce * (activation * activeForceLengthMultiplier *
                                            forceVelocityMultiplier);
        // conservative passive force
        conPassiveFiberForce = maxIsometricForce * normPassiveFiberForce;
        // non-conservative passive force
        nonConPassiveFiberForce =
                maxIsometricForce * get_fiber_damping() * normFiberVelocity;
        // total force
        totalFiberForce = activeFiberForce + conPassiveFiberForce +
                          nonConPassiveFiberForce;
    }

    /// The stiffness of the fiber in the direction of the fiber. This includes
    /// both active and passive force contributions to stiffness from the muscle
    /// fiber.
    /// @note based on Millard2012EquilibriumMuscle::calcFiberStiffness().
    SimTK::Real calcFiberStiffness(const SimTK::Real& activation,
            const SimTK::Real& normFiberLength,
            const SimTK::Real& fiberVelocityMultiplier) const {

        const SimTK::Real partialNormFiberLengthPartialFiberLength =
                1.0 / get_optimal_fiber_length();
        const SimTK::Real partialNormActiveForcePartialFiberLength =
                partialNormFiberLengthPartialFiberLength *
                calcActiveForceLengthMultiplierDerivative(normFiberLength);
        const SimTK::Real partialNormPassiveForcePartialFiberLength =
                partialNormFiberLengthPartialFiberLength *
                calcPassiveForceMultiplierDerivative(normFiberLength);

        // fiberStiffness = d_fiberForce / d_fiberLength
        return get_max_isometric_force() *
               (activation * partialNormActiveForcePartialFiberLength *
                               fiberVelocityMultiplier +
                       partialNormPassiveForcePartialFiberLength);
    }

    /// The stiffness of the tendon in the direction of the tendon.
    /// @note based on Millard2012EquilibriumMuscle.
    SimTK::Real calcTendonStiffness(const SimTK::Real& normTendonLength) const {

        if (get_ignore_tendon_compliance()) return SimTK::Infinity;
        return (get_max_isometric_force() / get_tendon_slack_length()) *
               calcTendonForceMultiplierDerivative(normTendonLength);
    }

    /// The stiffness of the whole musculotendon unit in the direction of the
    /// tendon.
    /// @note based on Millard2012EquilibriumMuscle.
    SimTK::Real calcMuscleStiffness(const SimTK::Real& tendonStiffness,
            const SimTK::Real& fiberStiffnessAlongTendon) const {

        if (get_ignore_tendon_compliance()) return fiberStiffnessAlongTendon;
        // TODO Millard2012EquilibriumMuscle includes additional checks that
        // the stiffness is non-negative and that the denominator is non-zero.
        // Checks are omitted here to preserve continuity and smoothness for
        // optimization (see #3685).
        return (fiberStiffnessAlongTendon * tendonStiffness) /
               (fiberStiffnessAlongTendon + tendonStiffness);
    }

    virtual double calcMuscleStiffness(const SimTK::State& s) const override
    {
        const MuscleDynamicsInfo& mdi = getMuscleDynamicsInfo(s);
        return calcMuscleStiffness(
                mdi.tendonStiffness,
                mdi.fiberStiffnessAlongTendon);
    }

    /// The derivative of pennation angle with respect to fiber length.
    /// @note based on
    /// MuscleFixedWidthPennationModel::calc_DPennationAngle_DFiberLength().
    SimTK::Real calcPartialPennationAnglePartialFiberLength(
            const SimTK::Real& fiberLength) const {

        using SimTK::square;
        // pennationAngle = asin(fiberWidth/fiberLength)
        // d_pennationAngle/d_fiberLength =
        //          d/d_fiberLength (asin(fiberWidth/fiberLength))
        return (-getFiberWidth() / square(fiberLength)) /
               sqrt(1.0 - square(getFiberWidth() / fiberLength));
    }

    /// The derivative of the fiber force along the tendon with respect to fiber
    /// length.
    /// @note based on
    /// Millard2012EquilibriumMuscle::calc_DFiberForceAT_DFiberLength().
    SimTK::Real calcPartialFiberForceAlongTendonPartialFiberLength(
            const SimTK::Real& fiberForce, const SimTK::Real& fiberStiffness,
            const SimTK::Real& sinPennationAngle,
            const SimTK::Real& cosPennationAngle,
            const SimTK::Real& partialPennationAnglePartialFiberLength) const {

        const SimTK::Real partialCosPennationAnglePartialFiberLength =
                -sinPennationAngle * partialPennationAnglePartialFiberLength;

        // The stiffness of the fiber along the direction of the tendon. For
        // small changes in length parallel to the fiber, this quantity is
        // d_fiberForceAlongTendon / d_fiberLength =
        //      d/d_fiberLength(fiberForce * cosPennationAngle)
        return fiberStiffness * cosPennationAngle +
               fiberForce * partialCosPennationAnglePartialFiberLength;
    }

    /// The derivative of the fiber force along the tendon with respect to the
    /// fiber length along the tendon.
    /// @note based on
    /// Millard2012EquilibriumMuscle::calc_DFiberForceAT_DFiberLengthAT.
    SimTK::Real calcFiberStiffnessAlongTendon(const SimTK::Real& fiberLength,
            const SimTK::Real& partialFiberForceAlongTendonPartialFiberLength,
            const SimTK::Real& sinPennationAngle,
            const SimTK::Real& cosPennationAngle,
            const SimTK::Real& partialPennationAnglePartialFiberLength) const {

        // The change in length of the fiber length along the tendon.
        // fiberLengthAlongTendon = fiberLength * cosPennationAngle
        const SimTK::Real partialFiberLengthAlongTendonPartialFiberLength =
                cosPennationAngle -
                fiberLength * sinPennationAngle *
                        partialPennationAnglePartialFiberLength;

        // fiberStiffnessAlongTendon
        //    = d_fiberForceAlongTendon / d_fiberLengthAlongTendon
        //    = (d_fiberForceAlongTendon / d_fiberLength) *
        //      (1 / (d_fiberLengthAlongTendon / d_fiberLength))
        return partialFiberForceAlongTendonPartialFiberLength *
               (1.0 / partialFiberLengthAlongTendonPartialFiberLength);
    }

    SimTK::Real calcPartialTendonLengthPartialFiberLength(
            const SimTK::Real& fiberLength,
            const SimTK::Real& sinPennationAngle,
            const SimTK::Real& cosPennationAngle,
            const SimTK::Real& partialPennationAnglePartialFiberLength) const {

        return fiberLength * sinPennationAngle *
                       partialPennationAnglePartialFiberLength -
               cosPennationAngle;
    }

    SimTK::Real calcPartialTendonForcePartialFiberLength(
            const SimTK::Real& tendonStiffness, const SimTK::Real& fiberLength, 
            const SimTK::Real& sinPennationAngle, 
            const SimTK::Real& cosPennationAngle) const {
        const SimTK::Real partialPennationAnglePartialFiberLength =
                calcPartialPennationAnglePartialFiberLength(fiberLength);

        const SimTK::Real partialTendonLengthPartialFiberLength =
                calcPartialTendonLengthPartialFiberLength(fiberLength, 
                    sinPennationAngle, cosPennationAngle,
                        partialPennationAnglePartialFiberLength);

        return tendonStiffness * partialTendonLengthPartialFiberLength;
    }
    /// @}

    /// @name Utilities
    /// @{

    /// Export the active force-length multiplier and passive force multiplier
    /// curves to a DataTable. If the normFiberLengths argument is omitted, we
    /// use createVectorLinspace(200, minNormFiberLength, maxNormFiberLength).
    DataTable exportFiberLengthCurvesToTable(
            const SimTK::Vector& normFiberLengths = SimTK::Vector()) const;
    /// Export the fiber force-velocity multiplier curve to a DataTable. If
    /// the normFiberVelocities argument is omitted, we use
    /// createVectorLinspace(200, -1.1, 1.1).
    DataTable exportFiberVelocityMultiplierToTable(
            const SimTK::Vector& normFiberVelocities = SimTK::Vector()) const;
    /// Export the fiber tendon force multiplier curve to a DataTable. If
    /// the normFiberVelocities argument is omitted, we use
    /// createVectorLinspace(200, 0.95, 1.6)
    DataTable exportTendonForceMultiplierToTable(
            const SimTK::Vector& normTendonLengths = SimTK::Vector()) const;
    /// Print the muscle curves to STO files. The files will be named as
    /// `<muscle-name>_<curve_type>.sto`.
    ///
    /// @param directory
    ///     The directory to which the data files should be written. Do NOT
    ///     include the filename. By default, the files are printed to the
    ///     current working directory.
    void printCurvesToSTOFiles(const std::string& directory = ".") const;

    /// Replace muscles of other types in the model with muscles of this type.
    /// Currently, only Millard2012EquilibriumMuscles and Thelen2003Muscles
    /// are replaced. For these two muscle classes, we copy property values into
    /// equivalent properties of the newly created DeGrooteFregly2016Muscle. 
    /// If the model has muscles of other types, an exception is
    /// thrown unless allowUnsupportedMuscles is true, in which a
    /// DeGrooteFregly2016Muscle is created using only the base Muscle class 
    /// property values. 
    /// Since the DeGrooteFregly2016Muscle implements tendon compliance dynamics
    /// with normalized tendon force as the state variable, this function
    /// ignores the 'default_fiber_length' property in replaced muscles.
    static void replaceMuscles(
            Model& model, bool allowUnsupportedMuscles = false);
    /// @}

    /// @name Scaling
    /// @{
    /// Adjust the properties of the muscle after the model has been scaled. The
    /// optimal fiber length and tendon slack length are each multiplied by the
    /// ratio of the current path length and the path length before scaling.
    void extendPostScale(
            const SimTK::State& s, const ScaleSet& scaleSet) override;
    /// @}

private:
    void constructProperties();

    void calcMuscleLengthInfoHelper(const SimTK::Real& muscleTendonLength, 
            MuscleLengthInfo& mli) const;
    void calcFiberVelocityInfoHelper(const SimTK::Real& muscleTendonVelocity, 
            const MuscleLengthInfo& mli, FiberVelocityInfo& fvi) const;
    void calcMuscleDynamicsInfoHelper(const SimTK::Real& activation,
            const MuscleLengthInfo& mli, const FiberVelocityInfo& fvi, 
            MuscleDynamicsInfo& mdi) const;
    void calcMusclePotentialEnergyInfoHelper(
            const MuscleLengthInfo& mli, MusclePotentialEnergyInfo& mpei) const;

    SimTK::Real getFiberWidth() const {
        const auto normFiberWidth = sin(get_pennation_angle_at_optimal());
        return get_optimal_fiber_length() * normFiberWidth;
    }
    SimTK::Real getSquareFiberWidth() const {
        return SimTK::square(getFiberWidth());
    }
    SimTK::Real getMaxContractionVelocityInMetersPerSecond() const {
        return get_max_contraction_velocity() * get_optimal_fiber_length();
    }

    // Curve parameters.
    // Notation comes from the Hyfydy documentation.

    // Parameters for the active fiber force-length curve.
    // ---------------------------------------------------
    constexpr static double cL1 = 1.5;
    constexpr static double cL2 = -2.75;
    constexpr static double m_minNormFiberLength = 0.46899;
    constexpr static double m_maxNormFiberLength = 1.80528;

    // Parameters for the passive fiber force-length curve.
    // ---------------------------------------------------
    constexpr static double cP1 = 1.08027;
    constexpr static double cP2 = 1.27368;

    // Parameters for the tendon force curve.
    // --------------------------------------
    constexpr static double cT1 = 260.972;
    constexpr static double cT2 = 7.9706;

    // Parameters for the force-velocity curve.
    // ----------------------------------------
    constexpr static double cV1 = 0.227;
    constexpr static double cV2 = 0.110;
    constexpr static double Fvmax = 1.6;

    static const std::string STATE_ACTIVATION_NAME;

    // Indices for MuscleDynamicsInfo::userDefinedDynamicsExtras.
    constexpr static int m_mdi_passiveFiberElasticForce = 0;
    constexpr static int m_mdi_passiveFiberDampingForce = 1;
    constexpr static int m_mdi_partialPennationAnglePartialFiberLength = 2;
    constexpr static int m_mdi_partialFiberForceAlongTendonPartialFiberLength =
            3;
    constexpr static int m_mdi_partialTendonForcePartialFiberLength = 4;
};

} // namespace OpenSim

#endif // HYFYDY_MUSCLE_H
