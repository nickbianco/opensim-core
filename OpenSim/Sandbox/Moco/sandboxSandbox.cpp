/* -------------------------------------------------------------------------- *
 * OpenSim Moco: sandboxSandbox.cpp                                           *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2019 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Christopher Dembia                                              *
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

// This file provides a way to easily prototype or test temporary snippets of
// code during development.

#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/Common/CommonUtilities.h>
#include <OpenSim/Simulation/Model/Scholz2015GeometryPath.h>

using namespace OpenSim;

std::vector<std::pair<double,double>>
computeDataPoints(const Model& model, std::string coordinatePath,
        std::string musclePath) {

    Model modelCopy(model);
    SimTK::State state = modelCopy.initSystem();
    const Coordinate& coordinate =
            modelCopy.getComponent<Coordinate>(coordinatePath);
    const Muscle& muscle = modelCopy.getComponent<Muscle>(musclePath);
    const double start = coordinate.getRangeMin();
    const double end = coordinate.getRangeMax();
    SimTK::Vector x = createVectorLinspace(202, start, end);

    std::vector<std::pair<double,double>> data_points;
    for (int i = 0; i < x.size(); ++i) {
        coordinate.setValue(state, x[i]);
        modelCopy.equilibrateMuscles(state);
        modelCopy.realizeReport(state);
        const double y = muscle.getPath().computeMomentArm(state, coordinate);
        data_points.emplace_back(x[i], y);
    }

    return data_points;
}

int main() {
    Model model("athlete_modified_right_ellipsoid.osim");
    model.initSystem();

    // BUG: If the default coordinate value is set to zero, the moment arm
    // values are negative for all coordinate values (which is the expected
    // behavior for a tricep muscle crossing the elbow). However, if the default
    // value is set to a larger value (e.g., Pi/1.25), then the moment arm
    // values are positive for all coordinate values.
    Coordinate& elbow_flex_r =
        model.updComponent<Coordinate>("/jointset/elbow_r/elbow_flex_r");
    elbow_flex_r.setDefaultValue(0);
    model.finalizeFromProperties();
    SimTK::State state = model.initSystem();

    model.realizePosition(state);
    const Muscle& muscle = model.getComponent<Muscle>("/forceset/triceps_long_r");
    const Scholz2015GeometryPath& path = muscle.getPath<Scholz2015GeometryPath>();
    double moment_arm = path.computeMomentArm(state, elbow_flex_r);
    std::cout << "state.getQ() = " << state.getQ() << "\n";
    std::cout << "moment arm at q_default: " << moment_arm << "\n";

    // Reset the cable warm-start to LiftedFromSurface so the solver finds the
    // correct topology at the new coordinate value, independent of what it
    // solved previously. setValue() will then invalidate the Position cache, and
    // the subsequent realizePosition() re-solves from that neutral start.
    path.resetWarmStart(state);

    elbow_flex_r.setValue(state, 0);
    model.realizePosition(state);
    moment_arm = path.computeMomentArm(state, elbow_flex_r);
    std::cout << "state.getQ() = " << state.getQ() << "\n";
    std::cout << "moment arm at q=0: " << moment_arm << "\n";

    // auto data_points = computeDataPoints(model,
    //         "/jointset/elbow_r/elbow_flex_r",
    //         "/forceset/triceps_long_r");

    // std::cout << "x,y\n";
    // for (const auto& [x, y] : data_points) {
    //     std::cout << x << "," << y << "\n";
    // }

    return EXIT_SUCCESS;
}
