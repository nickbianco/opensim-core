/* -------------------------------------------------------------------------- *
 *                        OpenSim:  HyfydyMuscle.cpp                          *
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

#include "HyfydyMuscle.h"

#include <OpenSim/Actuators/Millard2012EquilibriumMuscle.h>
#include <OpenSim/Actuators/Thelen2003Muscle.h>
#include <OpenSim/Common/CommonUtilities.h>
#include <OpenSim/Simulation/Model/Model.h>
#include "OpenSim/Common/STOFileAdapter.h"

using namespace OpenSim;

const std::string HyfydyMuscle::STATE_ACTIVATION_NAME("activation");

// We must define these variables in some compilation unit (pre-C++17).
// https://stackoverflow.com/questions/40690260/undefined-reference-error-for-static-constexpr-member?noredirect=1&lq=1
constexpr double HyfydyMuscle::cL1;
constexpr double HyfydyMuscle::cL2;
constexpr double HyfydyMuscle::cP1;
constexpr double HyfydyMuscle::cP2;
constexpr double HyfydyMuscle::cT1;
constexpr double HyfydyMuscle::cV1;
constexpr double HyfydyMuscle::cV2;
constexpr double HyfydyMuscle::Fvmax;
constexpr double HyfydyMuscle::m_minNormFiberLength;
constexpr double HyfydyMuscle::m_maxNormFiberLength;
constexpr int HyfydyMuscle::m_mdi_passiveFiberElasticForce;
constexpr int HyfydyMuscle::m_mdi_passiveFiberDampingForce;
constexpr int HyfydyMuscle::m_mdi_partialPennationAnglePartialFiberLength;
constexpr int HyfydyMuscle::m_mdi_partialFiberForceAlongTendonPartialFiberLength;
constexpr int HyfydyMuscle::m_mdi_partialTendonForcePartialFiberLength;

void HyfydyMuscle::constructProperties() {
    constructProperty_activation_rate(100.0);
    constructProperty_deactivation_rate(25.0);
    constructProperty_default_activation(0.5);
    constructProperty_fiber_damping(0.0);
    constructProperty_ignore_passive_fiber_force(false);
}

void HyfydyMuscle::extendFinalizeFromProperties() {
    Super::extendFinalizeFromProperties();
    OPENSIM_THROW_IF_FRMOBJ(!getProperty_optimal_force().getValueIsDefault(),
            Exception,
            "The optimal_force property is ignored for this Force; "
            "use max_isometric_force instead.");

    OPENSIM_THROW_IF_FRMOBJ(!get_ignore_tendon_compliance(),
            Exception,
            "The ignore_tendon_compliance property must be set to 'true' for "
            "this Muscle.");

    SimTK_ERRCHK2_ALWAYS(get_activation_rate() > 0,
            "HyfydyMuscle::extendFinalizeFromProperties",
            "%s: activation_rate must be greater than zero, "
            "but it is %g.",
            getName().c_str(), get_activation_rate());

    SimTK_ERRCHK2_ALWAYS(get_deactivation_rate() > 0,
            "HyfydyMuscle::extendFinalizeFromProperties",
            "%s: deactivation_rate must be greater than zero, "
            "but it is %g.",
            getName().c_str(), get_deactivation_rate());

    SimTK_ERRCHK2_ALWAYS(get_default_activation() > 0,
            "HyfydyMuscle::extendFinalizeFromProperties",
            "%s: default_activation must be greater than zero, "
            "but it is %g.",
            getName().c_str(), get_default_activation());

    SimTK_ERRCHK2_ALWAYS(get_fiber_damping() >= 0,
            "HyfydyMuscle::extendFinalizeFromProperties",
            "%s: fiber_damping must be greater than or equal to zero, "
            "but it is %g.",
            getName().c_str(), get_fiber_damping());

    OPENSIM_THROW_IF_FRMOBJ(
            get_pennation_angle_at_optimal() < 0 ||
                    get_pennation_angle_at_optimal() >
                            SimTK::Pi / 2.0 - SimTK::SignificantReal,
            InvalidPropertyValue,
            getProperty_pennation_angle_at_optimal().getName(),
            "Pennation angle at optimal fiber length must be in the range [0, "
            "Pi/2).");
}

void HyfydyMuscle::extendAddToSystem(
        SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);
    if (!get_ignore_activation_dynamics()) {
        addStateVariable(STATE_ACTIVATION_NAME, SimTK::Stage::Dynamics);
    }
}

void HyfydyMuscle::extendInitStateFromProperties(SimTK::State& s) const {
    Super::extendInitStateFromProperties(s);
    if (!get_ignore_activation_dynamics()) {
        setActivation(s, get_default_activation());
    }
}

void HyfydyMuscle::extendSetPropertiesFromState(
        const SimTK::State& s) {
    Super::extendSetPropertiesFromState(s);
    if (!get_ignore_activation_dynamics()) {
        set_default_activation(getActivation(s));
    }
}

void HyfydyMuscle::computeStateVariableDerivatives(
        const SimTK::State& s) const {

    // Activation dynamics.
    // --------------------
    if (!get_ignore_activation_dynamics()) {
        const auto& activation = getActivation(s);
        const auto& excitation = getControl(s);
        static const double c2 = get_deactivation_rate();
        static const double c1 = get_activation_rate() - c2;
        const SimTK::Real derivative = 
                (excitation - activation) * (c1*excitation + c2);
        setStateVariableDerivativeValue(s, STATE_ACTIVATION_NAME, derivative);
    }
}

double HyfydyMuscle::computeActuation(const SimTK::State& s) const {
    const auto& mdi = getMuscleDynamicsInfo(s);
    setActuation(s, mdi.tendonForce);
    return mdi.tendonForce;
}

void HyfydyMuscle::calcMuscleLengthInfoHelper(
        const SimTK::Real& muscleTendonLength, MuscleLengthInfo& mli) const {

    // Tendon.
    // -------
    mli.normTendonLength = 1.0;
    mli.tendonStrain = mli.normTendonLength - 1.0;
    mli.tendonLength = get_tendon_slack_length() * mli.normTendonLength;

    // Fiber.
    // ------
    mli.fiberLengthAlongTendon = muscleTendonLength - mli.tendonLength;
    mli.fiberLength = sqrt(
            SimTK::square(mli.fiberLengthAlongTendon) + getSquareFiberWidth());
    mli.normFiberLength = mli.fiberLength / get_optimal_fiber_length();

    // Pennation.
    // ----------
    mli.cosPennationAngle = mli.fiberLengthAlongTendon / mli.fiberLength;
    mli.sinPennationAngle = getFiberWidth() / mli.fiberLength;
    mli.pennationAngle = asin(mli.sinPennationAngle);

    // Multipliers.
    // ------------
    mli.fiberPassiveForceLengthMultiplier =
            calcPassiveForceMultiplier(mli.normFiberLength);
    mli.fiberActiveForceLengthMultiplier =
            calcActiveForceLengthMultiplier(mli.normFiberLength);
}

void HyfydyMuscle::calcFiberVelocityInfoHelper(
        const SimTK::Real& muscleTendonVelocity, const MuscleLengthInfo& mli,
        FiberVelocityInfo& fvi) const {

    fvi.normTendonVelocity = 0.0;
    fvi.tendonVelocity = get_tendon_slack_length() * fvi.normTendonVelocity;
    fvi.fiberVelocityAlongTendon =
            muscleTendonVelocity - fvi.tendonVelocity;
    fvi.fiberVelocity =
            fvi.fiberVelocityAlongTendon * mli.cosPennationAngle;
    fvi.normFiberVelocity =
            fvi.fiberVelocity / getMaxContractionVelocityInMetersPerSecond();
    fvi.fiberForceVelocityMultiplier =
            calcForceVelocityMultiplier(fvi.normFiberVelocity);
    
    const SimTK::Real tanPennationAngle =
            getFiberWidth() / mli.fiberLengthAlongTendon;
    fvi.pennationAngularVelocity =
            -fvi.fiberVelocity / mli.fiberLength * tanPennationAngle;
}

void HyfydyMuscle::calcMuscleDynamicsInfoHelper(
        const SimTK::Real& activation, const MuscleLengthInfo& mli, 
        const FiberVelocityInfo& fvi,
        MuscleDynamicsInfo& mdi) const {

    mdi.activation = activation;

    SimTK::Real activeFiberForce;
    SimTK::Real conPassiveFiberForce;
    SimTK::Real nonConPassiveFiberForce;
    SimTK::Real totalFiberForce;
    calcFiberForce(mdi.activation, mli.fiberActiveForceLengthMultiplier,
            fvi.fiberForceVelocityMultiplier,
            mli.fiberPassiveForceLengthMultiplier, fvi.normFiberVelocity,
            activeFiberForce, conPassiveFiberForce, nonConPassiveFiberForce,
            totalFiberForce);

    SimTK::Real passiveFiberForce =
            conPassiveFiberForce + nonConPassiveFiberForce;

    // When using a rigid tendon, avoid generating compressive fiber forces by
    // saturating the damping force produced by the parallel element.
    // Based on Millard2012EquilibriumMuscle::calcMuscleDynamicsInfo().
    if (get_ignore_tendon_compliance()) {
       if (totalFiberForce < 0) {
           totalFiberForce = 0.0;
           nonConPassiveFiberForce = -activeFiberForce - conPassiveFiberForce; 
           passiveFiberForce = conPassiveFiberForce + nonConPassiveFiberForce;
       }
    }

    // Compute force entries.
    // ----------------------
    const auto maxIsometricForce = get_max_isometric_force();
    mdi.fiberForce = totalFiberForce;
    mdi.activeFiberForce = activeFiberForce;
    mdi.passiveFiberForce = passiveFiberForce;
    mdi.normFiberForce = mdi.fiberForce / maxIsometricForce;
    mdi.fiberForceAlongTendon = mdi.fiberForce * mli.cosPennationAngle;
    mdi.normTendonForce = mdi.normFiberForce * mli.cosPennationAngle;
    mdi.tendonForce = mdi.fiberForceAlongTendon;

    // Compute stiffness entries.
    // --------------------------
    mdi.fiberStiffness = calcFiberStiffness(mdi.activation, mli.normFiberLength,
            fvi.fiberForceVelocityMultiplier);
    const auto& partialPennationAnglePartialFiberLength =
            calcPartialPennationAnglePartialFiberLength(mli.fiberLength);
    const auto& partialFiberForceAlongTendonPartialFiberLength =
            calcPartialFiberForceAlongTendonPartialFiberLength(mdi.fiberForce,
                    mdi.fiberStiffness, mli.sinPennationAngle,
                    mli.cosPennationAngle,
                    partialPennationAnglePartialFiberLength);
    mdi.fiberStiffnessAlongTendon = calcFiberStiffnessAlongTendon(
            mli.fiberLength, partialFiberForceAlongTendonPartialFiberLength,
            mli.sinPennationAngle, mli.cosPennationAngle,
            partialPennationAnglePartialFiberLength);
    mdi.tendonStiffness = calcTendonStiffness(mli.normTendonLength);

    const auto& partialTendonForcePartialFiberLength =
            calcPartialTendonForcePartialFiberLength(mdi.tendonStiffness,
                    mli.fiberLength, mli.sinPennationAngle,
                    mli.cosPennationAngle);

    // Compute power entries.
    // ----------------------
    // In order for the fiberPassivePower to be zero work, the non-conservative
    // passive fiber force is lumped into active fiber power. This is based on
    // the implementation in Millard2012EquilibriumMuscle (and verified over
    // email with Matt Millard).
    mdi.fiberActivePower = -(mdi.activeFiberForce + nonConPassiveFiberForce) *
                           fvi.fiberVelocity;
    mdi.fiberPassivePower = -conPassiveFiberForce * fvi.fiberVelocity;
    mdi.tendonPower = -mdi.tendonForce * fvi.tendonVelocity;

    mdi.userDefinedDynamicsExtras.resize(5);
    mdi.userDefinedDynamicsExtras[m_mdi_passiveFiberElasticForce] =
            conPassiveFiberForce;
    mdi.userDefinedDynamicsExtras[m_mdi_passiveFiberDampingForce] =
            nonConPassiveFiberForce;
    mdi.userDefinedDynamicsExtras
            [m_mdi_partialPennationAnglePartialFiberLength] =
            partialPennationAnglePartialFiberLength;
    mdi.userDefinedDynamicsExtras
            [m_mdi_partialFiberForceAlongTendonPartialFiberLength] =
            partialFiberForceAlongTendonPartialFiberLength;
    mdi.userDefinedDynamicsExtras[m_mdi_partialTendonForcePartialFiberLength] =
            partialTendonForcePartialFiberLength;
}

void HyfydyMuscle::calcMusclePotentialEnergyInfoHelper(
        const MuscleLengthInfo& mli, MusclePotentialEnergyInfo& mpei) const {

    // Based on Millard2012EquilibriumMuscle::calcMusclePotentialEnergyInfo().

    // Fiber potential energy.
    // -----------------------
    // TODO
    mpei.fiberPotentialEnergy = 0.0;

    // Tendon potential energy.
    // ------------------------
    mpei.tendonPotentialEnergy = 0;

    // Total potential energy.
    // -----------------------
    mpei.musclePotentialEnergy =
            mpei.fiberPotentialEnergy + mpei.tendonPotentialEnergy;
}

void HyfydyMuscle::calcMuscleLengthInfo(
        const SimTK::State& s, MuscleLengthInfo& mli) const {

    const auto& muscleTendonLength = getLength(s);
    calcMuscleLengthInfoHelper(muscleTendonLength, mli);
}

void HyfydyMuscle::calcFiberVelocityInfo(
        const SimTK::State& s, FiberVelocityInfo& fvi) const {

    const auto& mli = getMuscleLengthInfo(s);
    const auto& muscleTendonVelocity = getLengtheningSpeed(s);
    const auto& activation = getActivation(s);

    calcFiberVelocityInfoHelper(muscleTendonVelocity, mli, fvi);

    if (fvi.normFiberVelocity < -1.0) {
        log_info("HyfydyMuscle '{}' is exceeding maximum "
                 "contraction velocity at time {} s.",
                getName(), s.getTime());
    }
}

void HyfydyMuscle::calcMuscleDynamicsInfo(
        const SimTK::State& s, MuscleDynamicsInfo& mdi) const {
    const auto& activation = getActivation(s);
    const auto& mli = getMuscleLengthInfo(s);
    const auto& fvi = getFiberVelocityInfo(s);

    calcMuscleDynamicsInfoHelper(activation, mli, fvi, mdi);
}

void HyfydyMuscle::calcMusclePotentialEnergyInfo(
        const SimTK::State& s, MusclePotentialEnergyInfo& mpei) const {
    const MuscleLengthInfo& mli = getMuscleLengthInfo(s);
    calcMusclePotentialEnergyInfoHelper(mli, mpei);
}

double
OpenSim::HyfydyMuscle::calcInextensibleTendonActiveFiberForce(
        SimTK::State& s, double activation) const {
    MuscleLengthInfo mli;
    FiberVelocityInfo fvi;
    MuscleDynamicsInfo mdi;
    const double muscleTendonLength = getLength(s);
    const double muscleTendonVelocity = getLengtheningSpeed(s);

    calcMuscleLengthInfoHelper(muscleTendonLength, mli);
    calcFiberVelocityInfoHelper(muscleTendonVelocity, mli, fvi);
    calcMuscleDynamicsInfoHelper(activation, mli, fvi, mdi);

    return mdi.activeFiberForce;
}

void HyfydyMuscle::computeInitialFiberEquilibrium(SimTK::State& s) const {
    return;
}

double HyfydyMuscle::getPassiveFiberElasticForce(
        const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s)
            .userDefinedDynamicsExtras[m_mdi_passiveFiberElasticForce];
}
double HyfydyMuscle::getPassiveFiberElasticForceAlongTendon(
        const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s)
                   .userDefinedDynamicsExtras[m_mdi_passiveFiberElasticForce] *
           getMuscleLengthInfo(s).cosPennationAngle;
}
double HyfydyMuscle::getPassiveFiberDampingForce(
        const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s)
            .userDefinedDynamicsExtras[m_mdi_passiveFiberDampingForce];
}
double HyfydyMuscle::getPassiveFiberDampingForceAlongTendon(
        const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s)
                   .userDefinedDynamicsExtras[m_mdi_passiveFiberDampingForce] *
           getMuscleLengthInfo(s).cosPennationAngle;
}

DataTable HyfydyMuscle::exportFiberLengthCurvesToTable(
        const SimTK::Vector& normFiberLengths) const {
    SimTK::Vector def;
    const SimTK::Vector* x = nullptr;
    if (normFiberLengths.nrow()) {
        x = &normFiberLengths;
    } else {
        def = createVectorLinspace(
                200, m_minNormFiberLength, m_maxNormFiberLength);
        x = &def;
    }

    DataTable table;
    table.setColumnLabels(
            {"active_force_length_multiplier", "passive_force_multiplier"});
    SimTK::RowVector row(2);
    for (int irow = 0; irow < x->nrow(); ++irow) {
        const auto& normFiberLength = x->get(irow);
        row[0] = calcActiveForceLengthMultiplier(normFiberLength);
        row[1] = calcPassiveForceMultiplier(normFiberLength);
        table.appendRow(normFiberLength, row);
    }
    return table;
}

DataTable HyfydyMuscle::exportTendonForceMultiplierToTable(
        const SimTK::Vector& normTendonLengths) const {
    SimTK::Vector def;
    const SimTK::Vector* x = nullptr;
    if (normTendonLengths.nrow()) {
        x = &normTendonLengths;
    } else {
        // Evaluate the inverse of the tendon curve at y = 1.
        def = createVectorLinspace(200, 0.95, 1.6);
        x = &def;
    }

    DataTable table;
    table.setColumnLabels({"tendon_force_multiplier"});
    SimTK::RowVector row(1);
    for (int irow = 0; irow < x->nrow(); ++irow) {
        const auto& normTendonLength = x->get(irow);
        row[0] = calcTendonForceMultiplier(normTendonLength);
        table.appendRow(normTendonLength, row);
    }
    return table;
}

DataTable HyfydyMuscle::exportFiberVelocityMultiplierToTable(
        const SimTK::Vector& normFiberVelocities) const {
    SimTK::Vector def;
    const SimTK::Vector* x = nullptr;
    if (normFiberVelocities.nrow()) {
        x = &normFiberVelocities;
    } else {
        def = createVectorLinspace(200, -1.1, 1.1);
        x = &def;
    }

    DataTable table;
    table.setColumnLabels({"force_velocity_multiplier"});
    SimTK::RowVector row(1);
    for (int irow = 0; irow < x->nrow(); ++irow) {
        const auto& normFiberVelocity = x->get(irow);
        row[0] = calcForceVelocityMultiplier(normFiberVelocity);
        table.appendRow(normFiberVelocity, row);
    }
    return table;
}

void HyfydyMuscle::printCurvesToSTOFiles(
        const std::string& directory) const {
    const std::string prefix =
            directory + SimTK::Pathname::getPathSeparator() + getName();
    STOFileAdapter::write(exportFiberLengthCurvesToTable(),
            prefix + "_fiber_length_curves.sto");
    STOFileAdapter::write(exportFiberVelocityMultiplierToTable(),
            prefix + "_fiber_velocity_multiplier.sto");
    STOFileAdapter::write(exportTendonForceMultiplierToTable(),
            prefix + "_tendon_force_multiplier.sto");
}

void HyfydyMuscle::extendPostScale(
        const SimTK::State& s, const ScaleSet& scaleSet) {
    Super::extendPostScale(s, scaleSet);

    AbstractGeometryPath& path = updPath();
    if (path.getPreScaleLength(s) > 0.0)
    {
        double scaleFactor = path.getLength(s) / path.getPreScaleLength(s);
        upd_optimal_fiber_length() *= scaleFactor;
        upd_tendon_slack_length() *= scaleFactor;

        // Clear the pre-scale length that was stored in the path.
        path.setPreScaleLength(s, 0.0);
    }
}
