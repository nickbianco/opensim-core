# Apply a patch to a dependency's source tree.
#
# Usage (from an ExternalProject PATCH_COMMAND):
#   ${CMAKE_COMMAND} -DPATCH_FILE=<patch> -DSOURCE_DIR=<SOURCE_DIR>
#                    -P ApplyPatch.cmake
#
# Tolerates an already-patched tree: ExternalProject re-runs the patch step
# whenever the download or update step re-runs, and a second `git apply` of the
# same patch would otherwise fail the build.

if(NOT PATCH_FILE OR NOT SOURCE_DIR)
    message(FATAL_ERROR
            "ApplyPatch.cmake requires -DPATCH_FILE= and -DSOURCE_DIR=")
endif()

find_package(Git REQUIRED)

# `git apply --reverse --check` succeeds only if the patch is already applied.
execute_process(
    COMMAND "${GIT_EXECUTABLE}" apply --reverse --check "${PATCH_FILE}"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE already_applied
    OUTPUT_QUIET ERROR_QUIET)
if(already_applied EQUAL 0)
    message(STATUS "Patch already applied: ${PATCH_FILE}")
    return()
endif()

execute_process(
    COMMAND "${GIT_EXECUTABLE}" apply "${PATCH_FILE}"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE result
    ERROR_VARIABLE error)
if(NOT result EQUAL 0)
    message(FATAL_ERROR
            "Failed to apply ${PATCH_FILE} in ${SOURCE_DIR}: ${error}")
endif()
message(STATUS "Applied patch: ${PATCH_FILE}")
