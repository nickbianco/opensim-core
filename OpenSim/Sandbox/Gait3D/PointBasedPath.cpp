/* -------------------------------------------------------------------------- *
 *                        OpenSim:  PointBasedPath.cpp                        *
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

#include "PointBasedPath.h"
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PhysicalFrame.h>
#include <OpenSim/Simulation/Model/Station.h>
#include <OpenSim/Simulation/Model/ForceConsumer.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>
#include <OpenSim/Simulation/MomentArmSolver.h>
#include <SimTKsimbody.h>

using namespace OpenSim;

//=============================================================================
// CONSTRUCTORS
//=============================================================================

PointBasedPath::PointBasedPath() : AbstractGeometryPath() {
    setAuthors("Nicholas Bianco");
}

PointBasedPath::PointBasedPath(const std::string& name) : AbstractGeometryPath() {
    setAuthors("Nicholas Bianco");
    setName(name);
}

//=============================================================================
// PUBLIC METHODS
//=============================================================================

void PointBasedPath::addStation(const std::string& name, 
                              const PhysicalFrame& frame, 
                              const SimTK::Vec3& location) {
    Station* station = new Station(frame, location);
    station->setName(name);
    addComponent(station);
    appendSocketConnectee_stations(*station);
}

void PointBasedPath::invalidateGeometryCache(const SimTK::State& s) const {
    markCacheVariableInvalid(s, _lengthCV);
    markCacheVariableInvalid(s, _speedCV);
    markCacheVariableInvalid(s, _directionsCV);
    markCacheVariableInvalid(s, _pointsCV);
}

double PointBasedPath::getLength(const SimTK::State& s) const {
    if (isCacheVariableValid(s, _lengthCV)) {
        return getCacheVariableValue<double>(s, _lengthCV);
    }
    
    double length, speed;
    SimTK::Array_<SimTK::UnitVec3> directions;
    SimTK::Array_<SimTK::Vec3> points;
    computeGeometry(s, length, speed, directions, points);
    
    // Cache the results
    setCacheVariableValue(s, _lengthCV, length);
    setCacheVariableValue(s, _speedCV, speed);
    setCacheVariableValue(s, _directionsCV, directions);
    setCacheVariableValue(s, _pointsCV, points);
    markCacheVariableValid(s, _lengthCV);
    markCacheVariableValid(s, _speedCV);
    markCacheVariableValid(s, _directionsCV);
    markCacheVariableValid(s, _pointsCV);
    
    return length;
}

double PointBasedPath::getLengtheningSpeed(const SimTK::State& s) const {
    if (isCacheVariableValid(s, _speedCV)) {
        return getCacheVariableValue<double>(s, _speedCV);
    }
    
    // If length is cached but speed isn't, we need to recompute
    if (isCacheVariableValid(s, _lengthCV)) {
        double length, speed;
        SimTK::Array_<SimTK::UnitVec3> dirs;
        SimTK::Array_<SimTK::Vec3> s_G;
        computeGeometry(s, length, speed, dirs, s_G);
        
        // Update speed cache
        setCacheVariableValue(s, _speedCV, speed);
        markCacheVariableValid(s, _speedCV);
        
        return speed;
    }
    
    // Neither is cached, compute both
    double length, speed;
    SimTK::Array_<SimTK::UnitVec3> directions;
    SimTK::Array_<SimTK::Vec3> points;
    computeGeometry(s, length, speed, directions, points);
    
    // Cache the results
    setCacheVariableValue(s, _lengthCV, length);
    setCacheVariableValue(s, _speedCV, speed);
    setCacheVariableValue(s, _directionsCV, directions);
    setCacheVariableValue(s, _pointsCV, points);
    markCacheVariableValid(s, _lengthCV);
    markCacheVariableValid(s, _speedCV);
    markCacheVariableValid(s, _directionsCV);
    markCacheVariableValid(s, _pointsCV);
    
    return speed;
}

void PointBasedPath::produceForces(const SimTK::State& state,
                                   double tension,
                                   ForceConsumer& forceConsumer) const {
    if (tension <= 0.0) return;

    const auto& stations = getSocket<Station>("stations");
    int numPoints = stations.getNumConnectees();

    if (numPoints < 2) return;

    // Use cached geometry if available, otherwise compute
    double length, speed;
    SimTK::Array_<SimTK::UnitVec3> directions;
    SimTK::Array_<SimTK::Vec3> points;
    
    if (isCacheVariableValid(state, _lengthCV) && isCacheVariableValid(state, _speedCV)) {
        // Use cached values
        length = getCacheVariableValue<double>(state, _lengthCV);
        speed = getCacheVariableValue<double>(state, _speedCV);
        directions = getCacheVariableValue<SimTK::Array_<SimTK::UnitVec3>>(state, _directionsCV);
        points = getCacheVariableValue<SimTK::Array_<SimTK::Vec3>>(state, _pointsCV);
    } else {
        // Compute all geometry quantities
        computeGeometry(state, length, speed, directions, points);
        
        // Cache the results for future use
        setCacheVariableValue(state, _lengthCV, length);
        setCacheVariableValue(state, _speedCV, speed);
        setCacheVariableValue(state, _directionsCV, directions);
        setCacheVariableValue(state, _pointsCV, points);
        markCacheVariableValid(state, _lengthCV);
        markCacheVariableValid(state, _speedCV);
        markCacheVariableValid(state, _directionsCV);
        markCacheVariableValid(state, _pointsCV);
    }

    // Apply forces using the computed geometry
    for (int i = 0; i < static_cast<int>(directions.size()); ++i) {
        const SimTK::Vec3 f_G = tension * directions[i];
        SimTK::SpatialVec spatialVecLeft = 
            SimTK::SpatialVec(points[i] % f_G, f_G);
        SimTK::SpatialVec spatialVecRight = 
            -SimTK::SpatialVec(points[i+1] % f_G, f_G);

        forceConsumer.consumeBodySpatialVec(state, 
                stations.getConnectee(i).getParentFrame(), 
                spatialVecLeft);
        forceConsumer.consumeBodySpatialVec(state, 
                stations.getConnectee(i+1).getParentFrame(), 
                spatialVecRight);
    }
}

double PointBasedPath::computeMomentArm(const SimTK::State& s,
                                        const Coordinate& aCoord) const {
    // TODO: Implement moment arm calculation
    return 0.0;
}

//=============================================================================
// PRIVATE METHODS
//=============================================================================

void PointBasedPath::computeGeometry(const SimTK::State& s,
                                    double& length,
                                    double& speed,
                                    SimTK::Array_<SimTK::UnitVec3>& directions,
                                    SimTK::Array_<SimTK::Vec3>& points) const {
    const auto& stations = getSocket<Station>("stations");
    int numPoints = stations.getNumConnectees();

    // Initialize outputs
    length = 0.0;
    speed = 0.0;
    directions.clear();
    points.clear();

    if (numPoints < 2) {
        return;
    }

    // Pre-allocate arrays
    directions.resize(numPoints - 1);
    points.resize(numPoints);

    // Get the first point's location and velocity in ground
    const auto& station0 = stations.getConnectee(0);
    const SimTK::MobilizedBody& body0 = 
        station0.getParentFrame().getMobilizedBody();
    const SimTK::Transform& X_GB0 = body0.getBodyTransform(s);
    const SimTK::Vec3& p0 = station0.get_location();
    SimTK::Vec3 s0_G = X_GB0.R() * p0;
    SimTK::Vec3 p0_G = X_GB0.p() + s0_G;
    SimTK::Vec3 v0_G = body0.findStationVelocityInGround(s, p0);
    points[0] = s0_G;

    // Compute all geometric quantities in a single pass
    for (int i = 1; i < numPoints; ++i) {
        const auto& stationi = stations.getConnectee(i);
        const SimTK::MobilizedBody& bodyi = 
            stationi.getParentFrame().getMobilizedBody();
        const SimTK::Transform& X_GBi = bodyi.getBodyTransform(s);
        const SimTK::Vec3& pi = stationi.get_location();
        SimTK::Vec3 si_G = X_GBi.R() * pi;
        SimTK::Vec3 pi_G = X_GBi.p() + si_G;
        SimTK::Vec3 vi_G = bodyi.findStationVelocityInGround(s, pi);
        points[i] = si_G;

        const SimTK::Vec3& ri_G = pi_G - p0_G; // vector from point0 to pointi
        const SimTK::Real dist = ri_G.norm();  // distance between the points
        
        if (dist >= SimTK::SignificantReal) {
            // Compute direction vector
            const SimTK::UnitVec3 dir(ri_G);
            directions[i-1] = dir;
            
            // Accumulate length
            length += dist;
            
            // Accumulate lengthening speed
            speed += dot(vi_G - v0_G, dir);
        } else {
            // For very small distances, use a zero direction
            directions[i-1] = SimTK::UnitVec3(0, 0, 0);
        }

        // Update for next iteration
        s0_G = si_G;
        p0_G = pi_G;
        v0_G = vi_G;
    }
}

void PointBasedPath::extendAddToSystem(SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);
    
    // Add cache variables for geometry computation
    _lengthCV = addCacheVariable("length", 0.0, SimTK::Stage::Position);
    _speedCV = addCacheVariable("speed", 0.0, SimTK::Stage::Velocity);
    _directionsCV = addCacheVariable("directions", 
            SimTK::Array_<SimTK::UnitVec3>(), SimTK::Stage::Position);
    _pointsCV = addCacheVariable("points", 
            SimTK::Array_<SimTK::Vec3>(), SimTK::Stage::Position);
}

void PointBasedPath::generateDecorations(bool fixed, 
                                         const ModelDisplayHints& hints,
                                         const SimTK::State& state,
                                         SimTK::Array_<SimTK::DecorativeGeometry>& 
                                         appendToThis) const {
    // Don't generate fixed decorations
    if (fixed) return;

    const auto& stations = getSocket<Station>("stations");
    int numPoints = stations.getNumConnectees();

    if (numPoints < 2) return;

    // Get the first point's location in ground
    const auto& station0 = stations.getConnectee(0);
    const SimTK::MobilizedBody& body0 = 
        station0.getParentFrame().getMobilizedBody();
    const SimTK::Transform& X_GB0 = body0.getBodyTransform(state);
    const SimTK::Vec3& p0 = station0.get_location();
    SimTK::Vec3 p0_G = X_GB0.p() + X_GB0.R() * p0;

    // Draw line segments between consecutive points
    for (int i = 1; i < numPoints; ++i) {
        const auto& stationi = stations.getConnectee(i);
        const SimTK::MobilizedBody& bodyi = 
            stationi.getParentFrame().getMobilizedBody();
        const SimTK::Transform& X_GBi = bodyi.getBodyTransform(state);
        const SimTK::Vec3& pi = stationi.get_location();
        SimTK::Vec3 pi_G = X_GBi.p() + X_GBi.R() * pi;

        // Add line segment
        appendToThis.push_back(SimTK::DecorativeLine(p0_G, pi_G)
            .setLineThickness(4)
            .setColor(SimTK::Vec3(1, 0, 0))
            .setBodyId(0)
            .setIndexOnBody(i));

        // Update for next iteration
        p0_G = pi_G;
    }
} 