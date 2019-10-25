import unittest
import numpy as np
from pdb2sql import pdb2sql
from pdb2sql import transform

class TestTools(unittest.TestCase):

    def setUp(self):
        self.db = pdb2sql('./pdb/dummy_transform.pdb')
        self.xyz = self.db.get('x,y,z')

    def tearDown(self):
        self.db.close()

    def test_get_xyz(self):
        """Verfify getting xyz from sql"""

        result = transform._get_xyz(self.db)
        target = np.array([[1.,  0.,  0.], [-1.,  0.,  0.],
                           [ 0., 1.,  0.], [ 0., -1.,  0.],
                           [ 0., 0.,  1.], [ 0.,  0., -1.]])
        np.testing.assert_equal(result, target)

    def test_translation(self):
        """Verify sql translation"""
        trans_vec = np.array([1,1,1])
        target = np.array([[2.,   1.,  1.], [ 0.,  1.,  1.],
                           [ 1.,  2.,  1.], [ 1.,  0.,  1.],
                           [ 1.,  1.,  2.], [ 1.,  1.,  0.]])
        transform.translation(self.db, trans_vec)
        result = self.db.get('x,y,z')
        np.testing.assert_almost_equal(result, target)

    def test_rot_axis(self):
        """Verify sql rotation using axis and angle"""
        # rotate pi around x axis
        angle = np.pi
        axis = (1., 0., 0.)
        target = np.array([[1.,   0.,  0.], [-1.,  0.,  0.],
                           [ 0., -1.,  0.], [ 0.,  1.,  0.],
                           [ 0.,  0., -1.], [ 0.,  0.,  1.]])
        transform.rot_axis(self.db, axis, angle)
        result = self.db.get('x,y,z')
        np.testing.assert_almost_equal(result, target)

    def test_rot_xyz_around_axis(self):
        """Verify xyz values rot ation using axis and angle"""
        # rotate pi around x, y and z axis
        angle = np.pi
        axes_xyz = [(1., 0., 0.),
                    (0., 1., 0.),
                    (0., 0., 1.)]
        targets = [np.array([[1.,   0.,  0.], [-1.,  0.,  0.],
                            [ 0., -1.,  0.], [ 0.,  1.,  0.],
                            [ 0.,  0., -1.], [ 0.,  0.,  1.]]),
                   np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                             [ 0.,  1.,  0.], [ 0., -1.,  0.],
                             [ 0.,  0., -1.], [ 0.,  0.,  1.]]),
                   np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                             [ 0., -1.,  0.], [ 0.,  1.,  0.],
                             [ 0.,  0.,  1.], [ 0.,  0., -1.]])]
        for axis, target in zip(axes_xyz, targets):
            with self.subTest(axis=axis, target=target):
                xyz_rot = transform.rot_xyz_around_axis(self.xyz, axis, angle)
                np.testing.assert_almost_equal(xyz_rot, target)

    def test_get_rot_axis_angle(self):
        """Verify generation of random axis and angle"""
        # number of repeats
        n = 1000
        for i in range(n):
            with self.subTest(i=i):
                axis, angle = transform.get_rot_axis_angle()
                # axis verctor must be unit vector
                result = axis[0]**2 + axis[1]**2 + axis[2]**2
                target = 1.
                np.testing.assert_almost_equal(result, target)
                # angle in the range [0, 2π）
                self.assertTrue(0. <= angle < 2 * np.pi)

    def test_get_rot_axis_angle_seed(self):
        """Verify specific random seed"""
        seed = 2019
        axis1, angle1 = transform.get_rot_axis_angle(seed)
        axis2, angle2 = transform.get_rot_axis_angle(seed)
        self.assertEqual(axis1, axis2)
        self.assertEqual(angle1, angle2)

    def test_rot_euler(self):
        """Verify sql rotation using Euler angles"""
        # rotate pi around z axis
        alpha, beta, gamma = 0, 0, np.pi
        target = np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                           [ 0., -1.,  0.], [ 0.,  1.,  0.],
                           [ 0.,  0.,  1.], [ 0.,  0., -1.]])
        transform.rot_euler(self.db, alpha, beta, gamma)
        result = self.db.get('x,y,z')
        np.testing.assert_almost_equal(result, target)

    def test_rotation_euler(self):
        """Verify xyz values rotation using Euler angles"""
        alpha, beta, gamma = 0, 0, np.pi
        # rotate pi around x, y and z axis
        angles = [(np.pi, 0., 0.),
                  (0., np.pi, 0.),
                  (0., 0., np.pi)]
        targets = [np.array([[1.,   0.,  0.], [-1.,  0.,  0.],
                            [ 0., -1.,  0.], [ 0.,  1.,  0.],
                            [ 0.,  0., -1.], [ 0.,  0.,  1.]]),
                   np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                             [ 0.,  1.,  0.], [ 0., -1.,  0.],
                             [ 0.,  0., -1.], [ 0.,  0.,  1.]]),
                   np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                             [ 0., -1.,  0.], [ 0.,  1.,  0.],
                             [ 0.,  0.,  1.], [ 0.,  0., -1.]])]
        for angle, target in zip(angles, targets):
            with self.subTest(angle=angle, target=target):
                result = transform.rotation_euler(
                    self.xyz, angle[0], angle[1], angle[2])
                np.testing.assert_almost_equal(result, target)

    def test_rot_mat(self):
        """Verify sql roation using rotation matrix"""
        # rotate pi around z-axis
        theta = np.pi
        cosa = np.cos(theta)
        sina = np.sin(theta)
        rot_mat = np.array([[cosa, -sina, 0],
                           [sina, cosa, 0],
                           [0, 0, 1]])
        target = np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                           [ 0., -1.,  0.], [ 0.,  1.,  0.],
                           [ 0.,  0.,  1.], [ 0.,  0., -1.]])
        transform.rot_mat(self.db, rot_mat)
        result = self.db.get('x,y,z')
        np.testing.assert_almost_equal(result, target)

    def test_rotation_matrix(self):
        """Verify xyz values roation using rotation matrix"""
        theta = np.pi
        cosa = np.cos(theta)
        sina = np.sin(theta)
        # rotate pi around x, y and z axis
        rot_mats = [np.array([[1, 0, 0], [0, cosa, -sina], [0, sina, cosa]]),
                    np.array([[cosa, 0, sina], [0, 1, 0], [-sina, 0, cosa]]),
                    np.array([[cosa, -sina, 0], [sina, cosa, 0], [0, 0, 1]])]
        targets = [np.array([[1.,   0.,  0.], [-1.,  0.,  0.],
                            [ 0., -1.,  0.], [ 0.,  1.,  0.],
                            [ 0.,  0., -1.], [ 0.,  0.,  1.]]),
                   np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                             [ 0.,  1.,  0.], [ 0., -1.,  0.],
                             [ 0.,  0., -1.], [ 0.,  0.,  1.]]),
                   np.array([[-1.,  0.,  0.], [ 1.,  0.,  0.],
                             [ 0., -1.,  0.], [ 0.,  1.,  0.],
                             [0., 0., 1.], [0., 0., - 1.]])]
        for mat, target in zip(rot_mats, targets):
            with self.subTest(mat=mat, target=target):
                result = transform.rotate(self.xyz, mat)
                np.testing.assert_almost_equal(result, target)

        def test_rotation_matrix_center(self):
            """Verify specific rotation center"""
            # rotate pi around z-axis with rotation center [1,1,1,]
            xyz = np.array([0., 0., 0.])
            rot_mat = np.array([[cosa, -sina, 0], [sina, cosa, 0], [0, 0, 1]])
            centers = [np.array([1., 1., 1.]), [1., 1., 1.]]
            for center in centers:
                with self.subTest(center=center):
                    result = transform.rotate(xyz, rot_mat)
                    target = np.array([2., 2., 0.])
                    np.testing.assert_almost_equal(result, target)

if __name__ == "__main__":
    unittest.main()