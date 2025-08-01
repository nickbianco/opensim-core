/* -------------------------------------------------------------------------- *
 *                        OpenSim:  sandboxRajagopal.cpp                      *
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

#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Common/STOFileAdapter.h>
#include <OpenSim/Simulation/VisualizerUtilities.h>
#include <OpenSim/Actuators/ModelOperators.h>

#include "DiscreteController.h"

using namespace OpenSim;

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool visualize = false;
    double controlValue = 0.1;  // Default control value
    bool randomizeSpeeds = true;  // Default to randomizing speeds
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--visualize" || arg == "-v") {
            visualize = true;
        } else if (arg == "--control-value" || arg == "-c") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --control-value requires a value\n";
                return EXIT_FAILURE;
            }
            try {
                controlValue = std::stod(argv[++i]);
                if (controlValue < 0.0 || controlValue > 1.0) {
                    std::cerr << "Error: Control value must be between 0.0 and 1.0\n";
                    return EXIT_FAILURE;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid control value '" << argv[i] << "'\n";
                return EXIT_FAILURE;
            }
        } else if (arg == "--randomize-speeds" || arg == "-r") {
            randomizeSpeeds = true;
        } else if (arg == "--no-randomize-speeds" || arg == "-n") {
            randomizeSpeeds = false;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
            std::cout << "Options:\n";
            std::cout << "  --visualize, -v              Enable visualization\n";
            std::cout << "  --control-value VALUE, -c VALUE  Set control value\n";
            std::cout << "                                (0.0 to 1.0, default: 0.1)\n";
            std::cout << "  --randomize-speeds, -r        Randomize initial speeds (default)\n";
            std::cout << "  --no-randomize-speeds, -n     Use zero initial speeds\n";
            std::cout << "  --help, -h                   Show this help message\n";
            return EXIT_SUCCESS;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            std::cerr << "Use --help for usage information.\n";
            return EXIT_FAILURE;
        }
    }
    
    // Display control value
    std::cout << "Using control value: " << controlValue << "\n";
    
    // Display speed initialization setting
    if (randomizeSpeeds) {
        std::cout << "Using randomized initial speeds\n";
    } else {
        std::cout << "Using zero initial speeds\n";
    }

    // OpenSim::Logger::setLevel(Logger::Level::Error);

    // Load the model
    ModelProcessor modelProcessor = 
        ModelProcessor("RajagopalContact.osim") |
        // ModOpReplaceMusclesWithDeGrooteFregly2016() |
        ModOpIgnoreTendonCompliance(); 
    Model model = modelProcessor.process();

    DiscreteController* controller = new DiscreteController();
    controller->setName("controller");
    model.addController(controller);

    // Default state
    // -------------
    SimTK::State state = model.initSystem();
    model.getComponent<Coordinate>(
        "/jointset/ground_pelvis/pelvis_ty").setValue(state, 1.05);
    model.getComponent<Coordinate>(
        "/jointset/back/lumbar_extension").setValue(state, -SimTK::Pi/8);
    model.getComponent<Coordinate>(
        "/jointset/hip_l/hip_flexion_l").setValue(state,  
            SimTK::convertDegreesToRadians(30.0));
    model.getComponent<Coordinate>(
        "/jointset/hip_l/hip_adduction_l").setValue(state,  
            SimTK::convertDegreesToRadians(-10.0));
    model.getComponent<Coordinate>(
        "/jointset/hip_r/hip_flexion_r").setValue(state,  
            SimTK::convertDegreesToRadians(30.0));
    model.getComponent<Coordinate>(
        "/jointset/hip_r/hip_adduction_r").setValue(state,  
            SimTK::convertDegreesToRadians(-10.0));
    model.getComponent<Coordinate>(
        "/jointset/walker_knee_l/knee_angle_l").setValue(state,
            SimTK::convertDegreesToRadians(60.0));
    model.getComponent<Coordinate>(
        "/jointset/ankle_l/ankle_angle_l").setValue(state,
            SimTK::convertDegreesToRadians(-10.0));
    model.getComponent<Coordinate>(
        "/jointset/walker_knee_r/knee_angle_r").setValue(state,
            SimTK::convertDegreesToRadians(60.0));
    model.getComponent<Coordinate>(
        "/jointset/ankle_r/ankle_angle_r").setValue(state,
            SimTK::convertDegreesToRadians(-10.0));
    
    // Initialize speeds based on command line option
    if (randomizeSpeeds) {
        // Randomize the initial speeds.
        state.updU() = SimTK::Test::randVector(state.getNU());
    } else {
        // Use zero initial speeds.
        state.updU() = SimTK::Vector(state.getNU(), 0.0);
    }

    // Default controls
    // ----------------
    SimTK::Vector controls(model.getNumControls(), controlValue);
    controller->setDiscreteControls(state, controls);

    // Simulate
    // --------
    // model.equilibrateMuscles(state);
    SimTK::Real finalTime = 20.0;
    Manager manager(model);
    manager.setWriteToStorage(visualize);
    manager.setPerformAnalyses(visualize);
    manager.setIntegratorMethod(Manager::IntegratorMethod::CPodes);
    manager.setIntegratorAccuracy(1e-2);
    manager.setIntegratorMaximumStepSize(1e-2);
    manager.initialize(state);
    manager.integrate(finalTime);
    
    // Visualize    
    // --------
    if (visualize) {
        TimeSeriesTable table = manager.getStatesTable();
        VisualizerUtilities::showMotion(model, table);
    }


}