# CMakeLists.txt for Gait3D sandbox
# This file builds the sandboxGait3D executable with custom components

file(GLOB GAIT3D_FILES "Rajagopal22Muscles.cpp" "Rajagopal22MusclesContact.osim")
file(GLOB TO_COMPILE "sandboxGait3D.cpp" "sandboxRajagopal.cpp" "sandboxDoublePendulum.cpp")

# These executables are *not* tests.

set(SANDBOX_FILES
    # Add any data files here if needed
)

foreach(dataFile ${SANDBOX_FILES})
    configure_file(${dataFile} ${CMAKE_CURRENT_BINARY_DIR}/${dataFile} COPYONLY)
endforeach()

foreach(sourceFile ${TO_COMPILE})
    get_filename_component(_target_name ${sourceFile} NAME_WE)
    add_executable(${_target_name} EXCLUDE_FROM_ALL ${sourceFile} ${GAIT3D_FILES})
    target_link_libraries(${_target_name} osimTools osimSimulation osimActuators)
    set_target_properties(${_target_name} PROPERTIES
        FOLDER "Future sandbox"
    )
endforeach()

