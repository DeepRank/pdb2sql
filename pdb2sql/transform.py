import numpy as np
from time import time
from .pdb2sqlcore import pdb2sql

'''
This file contains several transformations of the
molecular coordinate that might be usefull during the
definition of the data set.
'''

########################################################################
# Translation
########################################################################
def translation(db, vect, **kwargs):
    xyz = _get_xyz(db, **kwargs)
    xyz += vect
    _update(db, xyz, **kwargs)

########################################################################
# Rotation using axis–angle presentation
# see https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
########################################################################
def rot_axis(db, axis, angle, **kwargs):
    xyz = _get_xyz(db, **kwargs)
    xyz = rot_xyz_around_axis(xyz, axis, angle)
    _update(db, xyz, **kwargs)

def get_rot_axis_angle(seed=None):
    """Get the rotation angle and axis.

    Args:
        seed(int): random seed for numpy

    Returns:
        list(float): axis of rotation
        float: angle of rotation
    """
    if seed is not None:
        np.random.seed(seed)

    # define the rotation axis
    # uniform distribution on a sphere
    # eq1,2 in http://mathworld.wolfram.com/SpherePointPicking.html
    u1, u2 = np.random.rand(), np.random.rand()
    theta = 2 * np.pi * u1      # [0, 2*pi)
    phi = np.arccos(2 * u2 - 1) # [0, pi]
    # eq19 in http://mathworld.wolfram.com/SphericalCoordinates.html
    axis = [np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi)]

    # define the rotation angle
    angle =  2 * np.pi * np.random.rand()

    return axis, angle

def rot_xyz_around_axis(xyz, axis, angle, center=None):
    """Get the rotated xyz.

    Args:
        xyz(np.array): original xyz coordinates
        axis (list(float)): axis of rotation
        angle (float): angle of rotation
        center (list(float)): center of rotation,
            defaults to the mean of input xyz.

    Returns:
        np.array: rotated xyz coordinates
    """
    # get the data
    ct, st = np.cos(angle), np.sin(angle)
    ux, uy, uz = axis

    # definition of the rotation matrix
    rot_mat = np.array([[ct + ux**2 * (1 - ct),
                         ux * uy * (1 - ct) - uz * st,
                         ux * uz * (1 - ct) + uy * st],
                        [uy * ux * (1 - ct) + uz * st,
                         ct + uy**2 * (1 - ct),
                         uy * uz * (1 - ct) - ux * st],
                        [uz * ux * (1 - ct) - uy * st,
                         uz * uy * (1 - ct) + ux * st,
                         ct + uz**2 * (1 - ct)]])

    # apply the rotation
    return rotate(xyz, rot_mat, center)

########################################################################
# Rotation using Euler anlges
# see https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations
########################################################################

def rot_euler(db, alpha, beta, gamma, **kwargs):
    """Rotate molecule from Euler rotation axis.

    Args:
        alpha (float): angle of rotation around the x axis
        beta (float): angle of rotation around the y axis
        gamma (float): angle of rotation around the z axis
        **kwargs: keyword argument to select the atoms.
            See pdb2sql.get()
    """
    xyz = _get_xyz(db, **kwargs)
    xyz = rotation_euler(xyz, alpha, beta, gamma)
    _update(db, xyz, **kwargs)

def rotation_euler(xyz, alpha, beta, gamma, center=None):

    # precomte the trig
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)

    # rotation matrices
    rx = np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]])
    ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
    rz = np.array([[cg, -sg, 0], [sg, cg, 0], [0, 0, 1]])

    # get rotation matrix
    rot_mat = np.dot(rz, np.dot(ry, rx))

    # apply the rotation
    return rotate(xyz, rot_mat, center)

########################################################################
# Rotation using provided rotation matrix
########################################################################

def rot_mat(db, mat, **kwargs):
    """Rotate molecule from a rotation matrix.

    Args:
        mat (np.array): 3x3 rotation matrix
        **kwargs: keyword argument to select the atoms.
            See pdb2sql.get()
    """
    xyz = _get_xyz(db, **kwargs)
    xyz = rotate(xyz, mat)
    _update(db, xyz, **kwargs)

def rotate(xyz, rot_mat, center=None):
    """[summary]

    Args:
        xyz(np.ndarray): x,y,z coordinates
        rot_mat(np.ndarray): rotation matrix
        center (list or np.ndarray, optional): rotation center.
            Defaults to None, i.e. using molecule center as rotation
            center.

    Raises:
        TypeError: Rotation center must be list or 1D np.ndarray.

    Returns:
        np.ndarray: x,y,z coordinates after rotation
    """
    # the default rotation center is the center of molecule itself.
    if center is None:
        center = np.mean(xyz, 0)

    if not isinstance(center, (list, np.ndarray)):
        raise TypeError("Rotation center must be list or 1D np.ndarray")

    return np.dot(rot_mat, (xyz - center).T).T + center

########################################################################
# helper functions
########################################################################
def _get_xyz(db, **kwargs):
    return np.array(db.get('x,y,z', **kwargs))

def _update(db, xyz, **kwargs):
    db.update('x,y,z', xyz, **kwargs)