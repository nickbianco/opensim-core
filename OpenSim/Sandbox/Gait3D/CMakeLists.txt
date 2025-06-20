# CMakeLists.txt for Gait3D sandbox
# This file builds the sandboxGait3D executable with custom components

file(GLOB GAIT3D_FILES "PointPathMuscle.cpp" "MobilityLinearDamper.cpp" "MobilityLinearStop.cpp" "PointBasedPath.cpp" "DiscreteController.cpp")
file(GLOB TO_COMPILE "sandboxGait3D.cpp")

# These executables are *not* tests.

set(SANDBOX_FILES
    # Add any data files here if needed
)

foreach(dataFile ${SANDBOX_FILES})
    configure_file(${dataFile} ${CMAKE_CURRENT_BINARY_DIR}/${dataFile} COPYONLY)
endforeach()

get_filename_component(_target_name ${TO_COMPILE} NAME_WE)
add_executable(${_target_name} EXCLUDE_FROM_ALL ${TO_COMPILE} ${GAIT3D_FILES})
target_link_libraries(${_target_name} osimTools osimSimulation osimActuators)
set_target_properties(${_target_name} PROPERTIES
    FOLDER "Future sandbox"
)
