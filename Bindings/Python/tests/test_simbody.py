"""Test Simbody classes.

"""

import os
import unittest

import numpy as np

import opensim as osim

resources_dir = os.path.join(os.path.dirname(os.path.abspath(osim.__file__)),
                             'tests', 'resources')

# Silence warning messages if mesh (.vtp) files cannot be found.
osim.Model.setDebugLevel(0)

class TestSimbody(unittest.TestCase):

    def test_vec3_typemaps(self):
        npv = np.array([5, 3, 6])
        v1 = osim.Vec3(npv)
        v2 = v1.to_numpy()
        assert (npv == v2).all()

        # Incorrect size.
        with self.assertRaises(RuntimeError):
            osim.Vec3(np.array([]))
        with self.assertRaises(RuntimeError):
            osim.Vec3(np.array([5, 1]))
        with self.assertRaises(RuntimeError):
            osim.Vec3(np.array([5, 1, 6, 3]))

        # createFromMat()
        npv = np.array([1, 6, 8])
        v1 = osim.Vec3.createFromMat(npv)
        v2 = v1.to_numpy()
        assert (npv == v2).all()

        # Incorrect size.
        with self.assertRaises(RuntimeError):
            osim.Vec3.createFromMat(np.array([]))
        with self.assertRaises(RuntimeError):
            osim.Vec3.createFromMat(np.array([5, 1]))
        with self.assertRaises(RuntimeError):
            osim.Vec3.createFromMat(np.array([5, 1, 6, 3]))

        # Incorrect number of args for `Vec2`.
        osim.Vec2(1.0, 2.0)  # This is fine
        with self.assertRaises(TypeError):
            osim.Vec2(1.0, 2.0, 3.0)
        with self.assertRaises(TypeError):
            osim.Vec2(1.0, 2.0, 3.0, 4.0)

        # Incorrect number of args for `Vec3`.
        with self.assertRaises(TypeError):
            osim.Vec3(1.0, 2.0)
        osim.Vec3(1.0, 2.0, 3.0)  # This is fine
        with self.assertRaises(TypeError):
            osim.Vec3(1.0, 2.0, 3.0, 4.0)
        with self.assertRaises(TypeError):
            osim.Vec3(1.0, 2.0, 3.0, 4.0, 5.0)

        # Incorrect number of args for `Vec4`.
        with self.assertRaises(TypeError):
            osim.Vec4(1.0, 2.0)
        with self.assertRaises(TypeError):
            osim.Vec4(1.0, 2.0, 3.0)
        osim.Vec4(1.0, 2.0, 3.0, 4.0)  # This is fine
        with self.assertRaises(TypeError):
            osim.Vec4(1.0, 2.0, 3.0, 4.0, 5.0)
        with self.assertRaises(TypeError):
            osim.Vec4(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

    def test_vec3_operators(self):
        v1 = osim.Vec3(1, 2, 3)
        # Tests __getitem__().
        assert v1[0] == 1
        assert v1[1] == 2

        # Out of bounds.
        with self.assertRaises(RuntimeError):
            v1[-1]
        with self.assertRaises(RuntimeError):
            v1[3]
        with self.assertRaises(RuntimeError):
            v1[5]

        # Tests __setitem__().
        v1[0] = 5
        assert v1[0] == 5

        # Out of bounds.
        with self.assertRaises(RuntimeError):
            v1[-1] = 5
        with self.assertRaises(RuntimeError):
            v1[3] = 1.3

        # Add. TODO removed for now.
        v2 = osim.Vec3(5, 6, 7)
        #v3 = v1 + v2
        #assert v3[0] == 10
        #assert v3[1] == 8
        #assert v3[2] == 10

        # Length.
        assert len(v2) == 3

    def test_mat33_and_rotation_to_numpy(self):
        rotation = osim.Rotation(0.3, osim.Vec3(0, 0, 1))
        expected = np.array([[rotation.get(i, j) for j in range(3)]
                             for i in range(3)])

        # Rotation is extended directly, and also exposes its Mat33.
        np.testing.assert_allclose(rotation.to_numpy(), expected, rtol=0, atol=0)
        np.testing.assert_allclose(rotation.asMat33().to_numpy(), expected,
                                   rtol=0, atol=0)

    def test_transform_to_numpy(self):
        rotation = osim.Rotation(0.4, osim.Vec3(1, 0, 0))
        translation = osim.Vec3(1.0, -2.0, 3.5)
        transform = osim.Transform(rotation, translation)

        mat = transform.to_numpy()
        assert mat.shape == (3, 4)
        np.testing.assert_allclose(mat[:, :3], rotation.to_numpy(),
                                   rtol=0, atol=0)
        np.testing.assert_allclose(mat[:, 3], translation.to_numpy(),
                                   rtol=0, atol=0)

    def test_vector_vec3_typemaps(self):
        npv = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v1 = osim.VectorVec3.createFromMat(npv.flatten())
        assert v1.size() == 2
        np.testing.assert_allclose(v1.to_numpy(), npv, rtol=0, atol=0)

        # Round trip through an element accessor, to confirm the packing order.
        assert v1.get(1)[0] == 4.0
        assert v1.get(1)[2] == 6.0

        # updFromMat() overwrites in place.
        v1.updFromMat(np.array([[7.0, 8.0, 9.0], [0.5, 0.25, 0.125]]).flatten())
        np.testing.assert_allclose(
            v1.to_numpy(), np.array([[7.0, 8.0, 9.0], [0.5, 0.25, 0.125]]),
            rtol=0, atol=0)

        # Empty.
        v2 = osim.VectorVec3.createFromMat(np.array([]))
        assert v2.size() == 0
        assert v2.to_numpy().shape == (0, 3)

        # Sizes that are not a multiple of three, or do not match.
        with self.assertRaises(RuntimeError):
            osim.VectorVec3.createFromMat(np.array([1.0, 2.0]))
        with self.assertRaises(RuntimeError):
            v1.updFromMat(np.array([1.0, 2.0, 3.0]))

    def test_simtk_array_vec3_typemaps(self):
        array = osim.SimTKArrayVec3()
        array.push_back(osim.Vec3(1.0, 2.0, 3.0))
        array.push_back(osim.Vec3(4.0, 5.0, 6.0))
        expected = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        np.testing.assert_allclose(array.to_numpy(), expected, rtol=0, atol=0)

        array.updFromMat(np.array([[9.0, 8.0, 7.0], [6.0, 5.0, 4.0]]).flatten())
        np.testing.assert_allclose(
            array.to_numpy(), np.array([[9.0, 8.0, 7.0], [6.0, 5.0, 4.0]]),
            rtol=0, atol=0)
        assert array.getElt(0)[0] == 9.0

        with self.assertRaises(RuntimeError):
            array.updFromMat(np.array([1.0, 2.0, 3.0]))

    def test_vector_spatialvec_typemaps(self):
        npv = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                        [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]])
        v1 = osim.VectorOfSpatialVec.createFromMat(npv.flatten())
        assert v1.size() == 2
        np.testing.assert_allclose(v1.to_numpy(), npv, rtol=0, atol=0)

        # Columns 0-2 are the first Vec3, columns 3-5 the second. SpatialVec is
        # not subscriptable from Python, so its halves are read with get().
        assert v1.get(0).get(0)[0] == 1.0
        assert v1.get(0).get(1)[0] == 4.0
        assert v1.get(1).get(1)[2] == 12.0

        v1.updFromMat(np.zeros(12))
        np.testing.assert_allclose(v1.to_numpy(), np.zeros((2, 6)),
                                   rtol=0, atol=0)

        with self.assertRaises(RuntimeError):
            osim.VectorOfSpatialVec.createFromMat(np.array([1.0, 2.0]))

    def test_vector_typemaps(self):
        npv = np.array([5, 3, 6, 2, 9])
        v1 = osim.Vector.createFromMat(npv)
        v2 = v1.to_numpy()
        assert (npv == v2).all()

        npv = np.array([])
        v1 = osim.Vector.createFromMat(npv)
        v2 = v1.to_numpy()
        assert (npv == v2).all()

    def test_rowvector_typemaps(self):
        npv = np.array([5, 3, 6, 2, 9])
        v1 = osim.RowVector.createFromMat(npv)
        v2 = v1.to_numpy()
        assert (npv == v2).all()

        npv = np.array([])
        v1 = osim.RowVector.createFromMat(npv)
        v2 = v1.to_numpy()
        assert (npv == v2).all()

    def test_vectorview_typemaps(self):
        # Use a TimeSeriesTable to obtain VectorViews.
        table = osim.TimeSeriesTable()
        table.setColumnLabels(['a', 'b'])
        table.appendRow(0.0, osim.RowVector([1.5, 2.0]))
        table.appendRow(1.0, osim.RowVector([2.5, 3.0]))
        column = table.getDependentColumn('a').to_numpy()
        assert len(column) == 2
        assert column[0] == 1.5
        assert column[1] == 2.5
        row = table.getRowAtIndex(0).to_numpy()
        assert len(row) == 2
        assert row[0] == 1.5
        assert row[1] == 2.0
        row = table.getRowAtIndex(1).to_numpy()
        assert len(row) == 2
        assert row[0] == 2.5
        assert row[1] == 3.0

    def test_matrix_typemaps(self):
        npm = np.array([[5, 3], [3, 6], [8, 1]])
        m1 = osim.Matrix.createFromMat(npm)
        m2 = m1.to_numpy()
        assert (npm == m2).all()

        npm = np.array([[]])
        m1 = osim.Matrix.createFromMat(npm)
        m2 = m1.to_numpy()
        assert (npm == m2).all()

        npm = np.array([])
        with self.assertRaises(TypeError):
            osim.Matrix.createFromMat(npm)

    def test_vector_operators(self):
        v = osim.Vector(5, 3)

        # Tests __getitem__()
        assert v[0] == 3
        assert v[4] == 3

        # Out of bounds.
        with self.assertRaises(RuntimeError):
            v[-1]
        with self.assertRaises(RuntimeError):
            v[5]

        # Tests __setitem__()
        v[0] = 15
        assert v[0] == 15

        with self.assertRaises(RuntimeError):
            v[-1] = 12
        with self.assertRaises(RuntimeError):
            v[5] = 14
        with self.assertRaises(RuntimeError):
            v[9] = 18

        # Size.
        assert len(v) == 5

    def test_exceptions(self):
        with self.assertRaises(RuntimeError):
            osim.Model("NONEXISTANT_FILE_NAME")
        with self.assertRaises(RuntimeError):
            m = osim.Model()
            # Asking for the visualizer just after constructing a model
            # throws an exception.
            m.getVisualizer()

    # def test_typemaps(self):
    #     # TODO disabled for now
    #     m = osim.Model()
    #     m.setGravity(osim.Vec3(1, 2, 3))
    #     m.setGravity([1, 2, 3])

    #     with self.assertRaises(ValueError):
    #         m.setGravity(['a', 2, 3])
    #     with self.assertRaises(ValueError):
    #         m.setGravity([1, 2])
    #     with self.assertRaises(ValueError):
    #         m.setGravity([1, 2, 6, 3])

    def test_printing(self):
        v1 = osim.Vec3(1, 3, 2)
        assert v1.__str__() == "~[1,3,2]"

        v2 = osim.Vector(7, 3)
        assert v2.__str__() == "~[3 3 3 3 3 3 3]"

    def test_SimbodyMatterSubsystem(self):
        model = osim.Model(os.path.join(resources_dir,
            "gait10dof18musc_subject01.osim"))
        s = model.initSystem()
        smss = model.getMatterSubsystem()
        
        assert smss.calcSystemMass(s) == model.getTotalMass(s)
        assert (smss.calcSystemMassCenterLocationInGround(s)[0] == 
                model.calcMassCenterPosition(s)[0])
        assert (smss.calcSystemMassCenterLocationInGround(s)[1] == 
                model.calcMassCenterPosition(s)[1])
        assert (smss.calcSystemMassCenterLocationInGround(s)[2] == 
                model.calcMassCenterPosition(s)[2])

        coordNames = model.getCoordinateNamesInMultibodyTreeOrder();
        print('firstCoord', coordNames.getElt(0));
        
        J = osim.Matrix()
        smss.calcSystemJacobian(s, J)
        # 6 * number of mobilized bodies
        assert J.nrow() == 6 * (model.getBodySet().getSize() + 1)
        assert J.ncol() == model.getCoordinateSet().getSize()

        v = osim.Vec3()
        smss.calcStationJacobian(s, 2, v, J)
        assert J.nrow() == 3
        assert J.ncol() == model.getCoordinateSet().getSize()

        # Inverse dynamics from SimbodyMatterSubsystem.calcResidualForce().
        # For the given inputs, we will actually be computing the first column
        # of the mass matrix. We accomplish this by setting all inputs to 0
        # except for the acceleration of the first coordinate.
        #   f_residual = M udot + f_inertial + f_applied 
        #              = M ~[1, 0, ...] + 0 + 0
        model.realizeVelocity(s)
        appliedMobilityForces = osim.Vector()
        appliedBodyForces = osim.VectorOfSpatialVec()
        knownUdot = osim.Vector(s.getNU(), 0.0); knownUdot[0] = 1.0
        knownLambda = osim.Vector()
        residualMobilityForces = osim.Vector()
        smss.calcResidualForce(s, appliedMobilityForces, appliedBodyForces,
                          knownUdot, knownLambda, residualMobilityForces)
        assert residualMobilityForces.size() == s.getNU()

        # Explicitly compute the first column of the mass matrix, then copmare.
        massMatrixFirstColumn = osim.Vector() 
        smss.multiplyByM(s, knownUdot, massMatrixFirstColumn)
        assert massMatrixFirstColumn.size() == residualMobilityForces.size()
        for i in range(massMatrixFirstColumn.size()):
            self.assertAlmostEqual(massMatrixFirstColumn[i],
                                   residualMobilityForces[i])

        # InverseDynamicsSolver.
        # Using accelerations from forward dynamics should give 0 residual.
        model.realizeAcceleration(s)
        idsolver = osim.InverseDynamicsSolver(model)
        residual = idsolver.solve(s, s.getUDot())
        assert residual.size() == s.getNU()
        for i in range(residual.size()):
            assert abs(residual[i]) < 1e-10

    def test_Mat33(self):
        mat33 = osim.Mat33(1, 2, 3, 4, 5, 6, 7, 8, 9)
        assert mat33.nrow() == 3
        assert mat33.ncol() == 3
        assert mat33.get(0, 0) == 1
        assert mat33.get(0, 1) == 2
        assert mat33.get(0, 2) == 3
        assert mat33.get(1, 0) == 4
        assert mat33.get(1, 1) == 5
        assert mat33.get(1, 2) == 6
        assert mat33.get(2, 0) == 7
        assert mat33.get(2, 1) == 8
        assert mat33.get(2, 2) == 9
