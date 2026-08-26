

set(MOCO_SOURCES
        osimMoco.h
        osimMocoDLL.h
        About.h
        About.cpp
        RegisterTypes_osimMoco.h
        RegisterTypes_osimMoco.cpp
        MocoUtilities.h
        MocoUtilities.cpp
        MocoStudy.h
        MocoStudy.cpp
        MocoBounds.h
        MocoBounds.cpp
        MocoVariableInfo.h
        MocoVariableInfo.cpp
        MocoProblem.h
        MocoProblem.cpp
        MocoProblemInfo.h
        MocoProblemRep.h
        MocoProblemRep.cpp
        MocoGoal/MocoGoal.h
        MocoGoal/MocoGoal.cpp
        MocoGoal/MocoMarkerFinalGoal.h
        MocoGoal/MocoMarkerFinalGoal.cpp
        MocoGoal/MocoMarkerTrackingGoal.h
        MocoGoal/MocoMarkerTrackingGoal.cpp
        MocoWeightSet.h
        MocoGoal/MocoStateTrackingGoal.h
        MocoGoal/MocoStateTrackingGoal.cpp
        MocoGoal/MocoSumSquaredStateGoal.h
        MocoGoal/MocoSumSquaredStateGoal.cpp
        MocoGoal/MocoContactTrackingGoal.cpp
        MocoGoal/MocoContactTrackingGoal.h
        MocoGoal/MocoContactImpulseTrackingGoal.cpp
        MocoGoal/MocoContactImpulseTrackingGoal.h
        MocoGoal/MocoControlGoal.h
        MocoGoal/MocoControlGoal.cpp
        MocoGoal/MocoControlTrackingGoal.h
        MocoGoal/MocoControlTrackingGoal.cpp
        MocoGoal/MocoExpressionBasedParameterGoal.h
        MocoGoal/MocoExpressionBasedParameterGoal.cpp
        MocoGoal/MocoJointReactionGoal.h
        MocoGoal/MocoJointReactionGoal.cpp
        MocoGoal/MocoOrientationTrackingGoal.h
        MocoGoal/MocoOrientationTrackingGoal.cpp
        MocoGoal/MocoOutputGoal.h
        MocoGoal/MocoOutputGoal.cpp
        MocoGoal/MocoTranslationTrackingGoal.h
        MocoGoal/MocoTranslationTrackingGoal.cpp
        MocoGoal/MocoAngularVelocityTrackingGoal.h
        MocoGoal/MocoAngularVelocityTrackingGoal.cpp
        MocoGoal/MocoAccelerationTrackingGoal.h
        MocoGoal/MocoAccelerationTrackingGoal.cpp
        MocoGoal/MocoInitialActivationGoal.cpp
        MocoGoal/MocoInitialActivationGoal.h
        MocoGoal/MocoInitialVelocityEquilibriumDGFGoal.h
        MocoGoal/MocoInitialVelocityEquilibriumDGFGoal.cpp
        MocoGoal/MocoInitialForceEquilibriumDGFGoal.h
        MocoGoal/MocoInitialForceEquilibriumDGFGoal.cpp
        MocoGoal/MocoPeriodicityGoal.h
        MocoGoal/MocoPeriodicityGoal.cpp
        MocoGoal/MocoStepTimeAsymmetryGoal.h
        MocoGoal/MocoStepTimeAsymmetryGoal.cpp
        MocoGoal/MocoStepLengthAsymmetryGoal.h
        MocoGoal/MocoStepLengthAsymmetryGoal.cpp
        MocoGoal/MocoGeneralizedForceTrackingGoal.h
        MocoGoal/MocoGeneralizedForceTrackingGoal.cpp
        MocoSolver.h
        MocoSolver.cpp
        MocoDirectCollocationSolver.h
        MocoDirectCollocationSolver.cpp
        MocoTrajectory.h
        MocoTrajectory.cpp
        MocoParameter.h
        MocoParameter.cpp
        MocoConstraint.h
        MocoConstraint.cpp
        MocoControlBoundConstraint.cpp
        MocoControlBoundConstraint.h
        MocoFrameDistanceConstraint.h
        MocoFrameDistanceConstraint.cpp
        MocoOutputConstraint.h
        MocoOutputConstraint.cpp
        Components/DiscreteForces.h
        Components/DiscreteForces.cpp
        Components/AccelerationMotion.h
        Components/AccelerationMotion.cpp
        Components/ControlDistributor.h
        Components/ControlDistributor.cpp
        Components/ActuatorInputController.h
        Components/ActuatorInputController.cpp
        MocoCasADiSolver/MocoCasADiSolver.h
        MocoCasADiSolver/MocoCasADiSolver.cpp
        MocoInverse.cpp
        MocoInverse.h
        MocoTrack.h
        MocoTrack.cpp
        ModelOperatorsDGF.h
        MocoTool.h
        MocoTool.cpp
        MocoConstraintInfo.h
        MocoConstraintInfo.cpp
        MocoStudyFactory.h
        MocoStudyFactory.cpp
        MocoScaleFactor.h
        MocoScaleFactor.cpp
        MocoOutputBoundConstraint.cpp
        MocoOutputBoundConstraint.h
        MocoStateBoundConstraint.cpp
        MocoStateBoundConstraint.h
        )
if(OPENSIM_WITH_CASADI)
    list(APPEND MOCO_SOURCES
            MocoCasADiSolver/CasOCProblem.h
            MocoCasADiSolver/CasOCProblem.cpp
            MocoCasADiSolver/CasOCSolver.h
            MocoCasADiSolver/CasOCSolver.cpp
            MocoCasADiSolver/CasOCFunction.h
            MocoCasADiSolver/CasOCFunction.cpp
            MocoCasADiSolver/CasOCTranscription.h
            MocoCasADiSolver/CasOCTranscription.cpp
            MocoCasADiSolver/CasOCTrapezoidal.h
            MocoCasADiSolver/CasOCTrapezoidal.cpp
            MocoCasADiSolver/CasOCHermiteSimpson.h
            MocoCasADiSolver/CasOCHermiteSimpson.cpp
            MocoCasADiSolver/CasOCLegendreGauss.h
            MocoCasADiSolver/CasOCLegendreGauss.cpp
            MocoCasADiSolver/CasOCLegendreGaussRadau.h
            MocoCasADiSolver/CasOCLegendreGaussRadau.cpp
            MocoCasADiSolver/CasOCIterate.h
            MocoCasADiSolver/MocoCasOCProblem.h
            MocoCasADiSolver/MocoCasOCProblem.cpp
            )
endif()

OpenSimAddLibrary(
    KIT Moco
    AUTHORS "Chris_Dembia-Nick_Bianco"
    LINKLIBS PUBLIC osimTools
    SOURCES ${MOCO_SOURCES}
    )

if (OPENSIM_WITH_CASADI)
    if (OPENSIM_CASADI_STATIC)
        get_target_property(casadi_type casadi TYPE)
        if (NOT casadi_type STREQUAL "STATIC_LIBRARY")
            message(FATAL_ERROR "OPENSIM_CASADI_STATIC requires a CasADi built "
                    "with ENABLE_STATIC=ON and ENABLE_SHARED=OFF; the CasADi "
                    "found at ${casadi_DIR} is a ${casadi_type}.")
        endif()
        # CasADi's plugins are archives rather than shared libraries in a static
        # build, so link the Ipopt plugin explicitly; osimMoco registers it at
        # runtime (see CasOCTranscription.cpp). Search only alongside the CasADi
        # we actually found, and only for a static archive: the dependencies
        # directory is on CMAKE_PREFIX_PATH and may hold a shared CasADi whose
        # plugins would pull a second CasADi back into the process.
        get_filename_component(casadi_lib_dir "${casadi_DIR}/../.." ABSOLUTE)
        set(casadi_saved_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
        find_library(CASADI_NLPSOL_IPOPT_LIBRARY casadi_nlpsol_ipopt
                PATHS "${casadi_lib_dir}" NO_DEFAULT_PATH)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${casadi_saved_suffixes})
        if (NOT CASADI_NLPSOL_IPOPT_LIBRARY)
            message(FATAL_ERROR "OPENSIM_CASADI_STATIC requires CasADi's static "
                    "Ipopt plugin (casadi_nlpsol_ipopt) in ${casadi_lib_dir}.")
        endif()
        # We bypass the imported CasADi target when linking so that we can hide
        # its symbols, so its usage requirements must be applied by hand.
        get_target_property(casadi_include_dirs casadi
                INTERFACE_INCLUDE_DIRECTORIES)
        target_include_directories(osimMoco SYSTEM PRIVATE
                ${casadi_include_dirs})
        get_target_property(casadi_compile_defs casadi
                INTERFACE_COMPILE_DEFINITIONS)
        if (casadi_compile_defs)
            target_compile_definitions(osimMoco PRIVATE ${casadi_compile_defs})
        endif()
        if (APPLE)
            set(casadi_link
                    "-Wl,-load_hidden,${CASADI_NLPSOL_IPOPT_LIBRARY}"
                    "-Wl,-load_hidden,$<TARGET_FILE:casadi>")
        elseif (MSVC)
            # Symbols are exported only via __declspec(dllexport) on Windows,
            # so a static CasADi leaks nothing into osimMoco's export table.
            set(casadi_link "${CASADI_NLPSOL_IPOPT_LIBRARY}" casadi)
        else()
            # Untested; --exclude-libs takes archive names, not paths.
            set(casadi_link "${CASADI_NLPSOL_IPOPT_LIBRARY}" casadi
                    "-Wl,--exclude-libs,libcasadi_nlpsol_ipopt.a"
                    "-Wl,--exclude-libs,libcasadi.a")
        endif()
        target_link_libraries(osimMoco PRIVATE ${casadi_link})
        # Ipopt stays a shared dependency; it exports no coalescable symbols.
        if (TARGET PkgConfig::IPOPT)
            target_link_libraries(osimMoco PRIVATE PkgConfig::IPOPT)
        endif()
        target_compile_definitions(osimMoco PRIVATE
                OPENSIM_CASADI_STATIC_PLUGINS)
    else()
        target_link_libraries(osimMoco PRIVATE casadi)
    endif()
endif()

target_include_directories(osimMoco INTERFACE
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if (OPENSIM_WITH_CASADI)
    target_compile_definitions(osimMoco PUBLIC OPENSIM_WITH_CASADI)
endif ()

# TODO: Remove the concept of a separate Moco version.
target_compile_definitions(osimMoco PRIVATE
    OPENSIM_MOCO_VERSION=${OPENSIM_QUALIFIED_VERSION}
    SWIG_PYTHON
    )

OpenSimAddInstallRPATHSelf(TARGET osimMoco LOADER)
OpenSimAddInstallRPATHSimbody(TARGET osimMoco LOADER
        FROM "${CMAKE_INSTALL_LIBDIR}")
install(TARGETS EXPORT OpenSimTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# TODO do not install InverseMuscleSolverMotionData.h b/c it uses Eigen.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenSim/Moco
        FILES_MATCHING PATTERN "*.h*"
        PATTERN "Archive" EXCLUDE
        PATTERN "doc" EXCLUDE
        PATTERN "Test" EXCLUDE
        )





if(BUILD_TESTING)
    add_subdirectory(Test)
endif()
