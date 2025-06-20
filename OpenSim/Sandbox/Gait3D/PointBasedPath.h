/* -------------------------------------------------------------------------- *
 *                        OpenSim:  PointBasedPath.h                          *
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

#ifndef OPENSIM_POINTBASEDPATH_H
#define OPENSIM_POINTBASEDPATH_H

#include <OpenSim/Simulation/Model/AbstractGeometryPath.h>
#include <OpenSim/Simulation/Model/Station.h>

namespace OpenSim {

/**
 * A concrete implementation of AbstractGeometryPath that uses a simple
 * point-based geometry similar to PointPathMuscle.
 * 
 * This class represents a path as a series of straight line segments
 * connecting points attached to PhysicalFrames. The path length is
 * computed as the sum of the distances between consecutive points,
 * and the lengthening speed is computed as the time derivative of
 * the path length.
 */
class PointBasedPath : public AbstractGeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(PointBasedPath, AbstractGeometryPath);

public:
    OpenSim_DECLARE_LIST_SOCKET(stations, Station, 
                                "The stations defining the path geometry.");

    PointBasedPath();
    PointBasedPath(const std::string& name);

    /**
     * Add a station to the path.
     * 
     * @param name     The name for the new station
     * @param frame    The PhysicalFrame to attach the station to
     * @param location The location of the station in the frame's coordinate system
     */
    void addStation(const std::string& name, 
                  const PhysicalFrame& frame, 
                  const SimTK::Vec3& location);

    /**
     * Invalidate the geometry cache. Call this when the path geometry changes.
     */
    void invalidateGeometryCache(const SimTK::State& s) const;

    // AbstractGeometryPath interface implementation
    double getLength(const SimTK::State& s) const override;
    double getLengtheningSpeed(const SimTK::State& s) const override;
    void produceForces(const SimTK::State& state,
                       double tension,
                       ForceConsumer& forceConsumer) const override;
    double computeMomentArm(const SimTK::State& s,
                           const Coordinate& aCoord) const override;
    bool isVisualPath() const override { return true; }

private:
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void generateDecorations(bool fixed, 
                            const ModelDisplayHints& hints,
                            const SimTK::State& state,
                            SimTK::Array_<SimTK::DecorativeGeometry>& 
                            appendToThis) const override;

    /**
     * Compute all geometric quantities efficiently in a single pass.
     * This method computes path length, lengthening speed, and direction
     * vectors for all segments, avoiding redundant computations.
     * 
     * @param s           The current state
     * @param length     Output: total path length
     * @param speed      Output: path lengthening speed
     * @param directions Output: direction vectors for each segment
     * @param points     Output: point vectors in ground frame for each station
     */
    void computeGeometry(const SimTK::State& s,
                        double& length,
                        double& speed,
                        SimTK::Array_<SimTK::UnitVec3>& directions,
                        SimTK::Array_<SimTK::Vec3>& points) const;

    // Cache variables to avoid recomputing geometry
    mutable CacheVariable<double> _lengthCV;
    mutable CacheVariable<double> _speedCV;
    mutable CacheVariable<SimTK::Array_<SimTK::UnitVec3>> _directionsCV;
    mutable CacheVariable<SimTK::Array_<SimTK::Vec3>> _pointsCV;
};

} // namespace OpenSim

#endif // OPENSIM_POINTBASEDPATH_H 