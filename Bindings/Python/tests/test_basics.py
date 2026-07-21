"""These are basic tests to make sure that C++ classes were wrapped properly.
There shouldn't be any python-specifc tests here.

"""

import os
import unittest

import opensim as osim

resources_dir = os.path.join(os.path.dirname(os.path.abspath(osim.__file__)),
                             'tests', 'resources')

# Silence warning messages if mesh (.vtp) files cannot be found.
osim.Model.setDebugLevel(0)

class TestBasics(unittest.TestCase):
    def test_version(self):
        print(osim.__version__)

    def test_set_mobilizer_frame_translations(self):
        # Setting many mobilizer frames at once must match setting them one at a
        # time through the Joint interface.
        import numpy as np

        def build():
            model = osim.Model()
            previous = model.getGround()
            for i in range(4):
                body = osim.Body(f'b{i}', 1.0, osim.Vec3(0),
                                 osim.Inertia(1, 1, 1, 0, 0, 0))
                model.addBody(body)
                joint = osim.PinJoint(f'j{i}', previous, osim.Vec3(0.1 * i, 0, 0),
                                      osim.Vec3(0), body, osim.Vec3(0, -0.4, 0),
                                      osim.Vec3(0))
                model.addJoint(joint)
                previous = body
            model.finalizeConnections()
            return model

        rng = np.random.default_rng(0)
        inboard = rng.uniform(-0.5, 0.5, (4, 3))
        outboard = rng.uniform(-0.5, 0.5, (4, 3))

        # Reference: one Joint at a time, preserving each frame's rotation.
        reference = build()
        state = reference.initSystem()
        reference.realizePosition(state)
        indexes = osim.SimTKArrayInt()
        for i in range(reference.getNumJoints()):
            joint = reference.getJointSet().get(i)
            indexes.push_back(int(joint.getChildFrame().getMobilizedBodyIndex()))
        for i in range(reference.getNumJoints()):
            joint = reference.getJointSet().get(i)
            X_PF = joint.getInboardFrame(state)
            joint.setInboardFrame(state, osim.Transform(
                X_PF.R(), osim.Vec3(*[float(v) for v in inboard[i]])))
        reference.realizePosition(state)
        for i in range(reference.getNumJoints()):
            joint = reference.getJointSet().get(i)
            X_BM = joint.getOutboardFrame(state)
            joint.setOutboardFrame(state, osim.Transform(
                X_BM.R(), osim.Vec3(*[float(v) for v in outboard[i]])))
        reference.realizePosition(state)
        expected = np.array(
            [reference.getBodySet().get(i).getPositionInGround(state).to_numpy()
             for i in range(reference.getNumBodies())])

        # Bulk: two calls, each taking every frame at once. The rotations are
        # supplied rather than read back, so they are captured from a pristine
        # model before anything is modified.
        model = build()
        bulk_state = model.initSystem()
        model.realizePosition(bulk_state)
        inboard_rotations = osim.SimTKArrayRotation()
        outboard_rotations = osim.SimTKArrayRotation()
        for i in range(model.getNumJoints()):
            joint = model.getJointSet().get(i)
            inboard_rotations.push_back(
                osim.Rotation(joint.getInboardFrame(bulk_state).R()))
            outboard_rotations.push_back(
                osim.Rotation(joint.getOutboardFrame(bulk_state).R()))
        model.setInboardFrames(
            bulk_state, indexes, inboard_rotations,
            osim.Vector.createFromMat(inboard.flatten()))
        model.setOutboardFrames(
            bulk_state, indexes, outboard_rotations,
            osim.Vector.createFromMat(outboard.flatten()))
        model.realizePosition(bulk_state)
        got = np.array(
            [model.getBodySet().get(i).getPositionInGround(bulk_state).to_numpy()
             for i in range(model.getNumBodies())])

        assert np.array_equal(got, expected), f'{got} != {expected}'

        # The call invalidates Stage::Instance but reads nothing, so calling it
        # twice in a row from an unrealized State has to succeed.
        model.setOutboardFrames(
            bulk_state, indexes, outboard_rotations,
            osim.Vector.createFromMat(outboard.flatten()))
        model.setOutboardFrames(
            bulk_state, indexes, outboard_rotations,
            osim.Vector.createFromMat(outboard.flatten()))
        model.realizePosition(bulk_state)
        np.testing.assert_allclose(
            np.array([model.getBodySet().get(i).getPositionInGround(
                bulk_state).to_numpy() for i in range(model.getNumBodies())]),
            expected, rtol=0, atol=0)

        # Mismatched translation count, and mismatched rotation count.
        with self.assertRaises(RuntimeError):
            model.setInboardFrames(bulk_state, indexes, inboard_rotations,
                                   osim.Vector.createFromMat(np.zeros(5)))
        with self.assertRaises(RuntimeError):
            model.setInboardFrames(bulk_state, indexes,
                                   osim.SimTKArrayRotation(),
                                   osim.Vector.createFromMat(inboard.flatten()))

    def test_muscle_helper_classes(self):
        # This test exists because some classes that Thelen2003Muscle used were
        # not accessibly in the bindings.
        muscle = osim.Thelen2003Muscle()

        fwpm = muscle.getPennationModel()
        fwpm.get_optimal_fiber_length()

        adm = muscle.getActivationModel()
        adm.get_activation_time_constant()

        muscle = osim.Millard2012EquilibriumMuscle()

        tendonFL = osim.TendonForceLengthCurve()
        muscle.setTendonForceLengthCurve(tendonFL)

    def test_SimTKArray(self):
        # Initally created to test the creation of a separate simbody module.
        ad = osim.SimTKArrayDouble()
        ad.push_back(1)

        av3 = osim.SimTKArrayVec3()
        av3.push_back(osim.Vec3(8))
        assert av3.at(0).get(0) == 8

    def test_ToolAndModel(self):
        # Test tools module.
        cmc = osim.CMCTool()
        model = osim.Model()
        model.setName('alphabet')
        cmc.setModel(model)
        assert cmc.getModel().getName() == 'alphabet'

    def test_AnalysisToolModel(self):
        # Test analyses module.
        cmc = osim.CMCTool()
        model = osim.Model()
        model.setName('eggplant')
        fr = osim.ForceReporter()
        fr.setName('strong')
        cmc.setModel(model)
        cmc.getAnalysisSet().adoptAndAppend(fr)

        assert cmc.getModel().getName() == 'eggplant'
        assert cmc.getAnalysisSet().get(0).getName() == 'strong'

    def test_ManagerConstructorCreatesIntegrator(self):
        # Make sure that the Manager is able to create a default integrator.
        # This tests a bug fix: previously, it was impossible to use the
        # Manager to integrate from MATLAB/Python, since it was not possible
        # to provide an Integrator to the Manager.
        model = osim.Model(os.path.join(resources_dir, "arm26.osim"))
        state = model.initSystem()

        manager = osim.Manager(model)
        state.setTime(0);
        manager.initialize(state);
        state = manager.integrate(0.00001)

    def test_WrapObject(self):
        # Make sure the WrapObjects are accessible.
        model = osim.Model()

        sphere = osim.WrapSphere()
        model.getGround().addWrapObject(sphere)

        cylinder = osim.WrapCylinder()
        cylinder.set_radius(0.5)
        model.getGround().addWrapObject(cylinder)

        torus = osim.WrapTorus()
        model.getGround().addWrapObject(torus)

        ellipsoid = osim.WrapEllipsoid()
        model.getGround().addWrapObject(ellipsoid)

    def test_ToyReflexController(self):
        controller = osim.ToyReflexController()
        
    def test_GCVSplineSet(self):
        splineset = osim.GCVSplineSet(os.path.join(resources_dir,
            'std_subject01_walk1_ik.mot'))
        splineset = osim.GCVSplineSet(
                osim.TimeSeriesTable(os.path.join(resources_dir,
                    'std_subject01_walk1_ik.mot')), [], 5, 0)

    def test_deserialize_tool_with_empty_model_file(self):
        # Ensure an exception is thrown when loading an AbstractTool (e.g.,
        # ForwardTool) setup file with an empty model_file. In particular, we
        # want to check the case where force_set_files is not empty.
        with self.assertRaises(RuntimeError):
            rra = osim.ForwardTool(os.path.join(resources_dir,
                'gait2392_setup_forward_empty_model.xml'))

        # No exception if we pass loadModel=False
        rra = osim.ForwardTool(
                os.path.join(resources_dir,
                    'gait2392_setup_forward_empty_model.xml'),
                True, # updateFromXMLNode
                False, # loadModel
                )

    def test_property_helper(self):
        muscle = osim.Thelen2003Muscle()
        # set max_isometric_force property using PropertyHelper
        # then retrieve using get_max_isometric_force native method
        property = muscle.getPropertyByName('max_isometric_force')
        osim.PropertyHelper.setValueDouble(200, property);
        assert muscle.get_max_isometric_force()==200

    def test_StatesTrajectory_and_StatesDocument(self):
        model = osim.ModelFactory.createDoublePendulum()
        state = model.initSystem()

        # Run a forward simulation and retreive the StatesTrajectory.
        manager = osim.Manager(model);
        manager.setRecordStatesTrajectory(True);
        manager.initialize(state);
        manager.integrate(5.0);
        statesTraj = manager.getStatesTrajectory();

        # Check initial time.
        initialState = statesTraj.get(0);
        self.assertEqual(initialState.getTime(), 0)
        self.assertEqual(initialState.getQ().size(), 2)
        self.assertEqual(initialState.getU().size(), 2)

        # Check final time.
        finalState = statesTraj.get(statesTraj.getSize()-1);
        self.assertEqual(finalState.getTime(), 5.0)
        self.assertEqual(finalState.getQ().size(), 2)
        self.assertEqual(finalState.getU().size(), 2)

        # Create a StatesDocument to serialize the trajectory.
        doc = statesTraj.exportToStatesDocument(model);
        doc.serialize('pendulum.ostates');

        # Deserialize and check equality.
        statesTrajDeserialized = osim.StatesTrajectory.createFromStatesDocument(
                model, 'pendulum.ostates')
        self.assertEqual(statesTraj.getSize(), statesTrajDeserialized.getSize())
        for i in range(statesTraj.getSize()):
            state = statesTraj.get(i)
            stateDeserialized = statesTrajDeserialized.get(i)
            self.assertEqual(state.getTime(), stateDeserialized.getTime())
