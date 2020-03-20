
import os
import numpy as np
from .pdb2sqlcore import pdb2sql
from .transform import rotate


def superpose(pdb1, pdb2, method='svd', only_backbone=True, **kwargs):
    """superpose two complexes

    Arguments:
        pdb1 {str or pdb2sql} -- name or sqldb of the first pdb
        pdb2 {str or pdb2sql} -- name or sqldb of the second pdb

    Keyword Arguments:
        method {str} -- method used to superpose the complex (default: {'svd'})
        only_backbone {bool} -- use only backbone atos to align (default: True)
        **kwargs -- keyword arguments used to select a portion of the pdb

    Example:
        >> pdb1 = '1AK4_5w.pdb'
        >> pdb2 = '1AK4_10w.pdb'
        >> superpose(pdb1, pdb2, chainID='A')

    """

    backbone_atoms = ['CA', 'C', 'N', 'O']

    if not isinstance(pdb1, pdb2sql):
        sql1 = pdb2sql(pdb1)
    else:
        sql1 = pdb1

    if not isinstance(pdb2, pdb2sql):
        sql2 = pdb2sql(pdb2)
    else:
        sql2 = pdb2

    if only_backbone:
        if 'name' not in kwargs:
            kwargs['name'] = backbone_atoms
        else:
            raise ValueError('Atom type specified but only_backbone == True')

    # xyz of the chains selected
    chain_xyz1 = np.array(sql1.get("x,y,z", **kwargs))
    chain_xyz2 = np.array(sql2.get("x,y,z", **kwargs))

    # translation vector
    tr1 = get_trans_vect(chain_xyz1)
    tr2 = get_trans_vect(chain_xyz2)

    # rotation matrix
    chain_xyz1 += tr1
    chain_xyz2 += tr2
    center = np.mean(chain_xyz2, 0)
    rmat = get_rotation_matrix(chain_xyz2, chain_xyz1, method=method)

    # transform the coordinate of second pdb
    xyz2 = np.array(sql2.get("x,y,z"))
    xyz2 += tr2
    xyz2 = rotate(xyz2, rmat, center=center)
    xyz2 -= tr1

    # update the second sql
    sql2.update('x,y,z', xyz2)
    pdb1_name = os.path.basename(pdb1)
    fname = pdb2.rstrip('.pdb') + '_superposed_on_' + pdb1_name.rstrip('.pdb') + '.pdb'
    sql2.exportpdb(fname)

# compute the translation vector to center a set of points


def get_trans_vect(P):
    """Get the translationv vector to the origin.

    Args:
        P (np.array(nx3)): position of the points in the molecule

    Returns:
        float: minus mean value of the xyz columns
    """
    return -np.mean(P, 0)

# main switch for the rotation matrix
# add new methods here if necessary


def get_rotation_matrix(P, Q, method='svd'):

    # get the matrix with Kabsh method
    if method.lower() == 'svd':
        mat = get_rotation_matrix_Kabsh(P, Q)

    # or with the quaternion method
    elif method.lower() == 'quaternion':
        mat = get_rotation_matrix_quaternion(P, Q)

    else:
        raise ValueError(
            f'{method} is not a valid method for rmsd alignement. '
            f'Options are svd or quaternions')

    return mat

# get the rotation matrix via a SVD
# decomposition of the correlation matrix


def get_rotation_matrix_Kabsh(P, Q):
    """Get the rotation matrix to aligh two point clouds.

    The method is based on th Kabsh approach
    https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures

    Args:
        P (np.array): xyz of the first point cloud
        Q (np.array): xyz of the second point cloud

    Returns:
        np.array: rotation matrix

    Raises:
        ValueError: matrix have different sizes
    """
    pshape = P.shape
    qshape = Q.shape

    if pshape[0] == qshape[0]:
        npts = pshape[0]
    else:
        raise ValueError("Matrix don't have the same number of points",
                         P.shape, Q.shape)

    p0, q0 = np.abs(np.mean(P, 0)), np.abs(np.mean(Q, 0))
    eps = 1E-6
    if any(p0 > eps) or any(q0 > eps):
        raise ValueError('You must center the fragment first', p0, q0)

    # form the covariance matrix
    A = np.dot(P.T, Q) / npts

    # SVD the matrix
    V, _, W = np.linalg.svd(A)

    # the W matrix returned here is
    # already its transpose
    # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.linalg.svd.html
    W = W.T

    # determinant
    d = np.linalg.det(np.dot(W, V.T))

    # form the U matrix
    Id = np.eye(3)
    if d < 0:
        Id[2, 2] = -1

    U = np.dot(W, np.dot(Id, V.T))

    return U

# get the rotation amtrix via the quaternion approach
# doesn't work great so far


def get_rotation_matrix_quaternion(P, Q):
    """Get the rotation matrix to aligh two point clouds.

    The method is based on the quaternion approach
    http://www.ams.stonybrook.edu/~coutsias/papers/rmsd17.pdf

    Args:
        P (np.array): xyz of the first point cloud
        Q (np.array): xyz of the second point cloud

    Returns:
        np.array: rotation matrix

    Raises:
        ValueError: matrix have different sizes
    """
    pshape = P.shape
    qshape = Q.shape

    if pshape[0] != qshape[0]:
        raise ValueError("Matrix don't have the same number of points",
                         P.shape, Q.shape)

    p0, q0 = np.abs(np.mean(P, 0)), np.abs(np.mean(Q, 0))
    eps = 1E-6
    if any(p0 > eps) or any(q0 > eps):
        raise ValueError('You must center the fragment first', p0, q0)

    # form the correlation matrix
    R = np.dot(P.T, Q)

    # form the F matrix (eq. 10 of ref[1])
    F = np.zeros((4, 4))

    F[0, 0] = np.trace(R)
    F[0, 1] = R[1, 2] - R[2, 1]
    F[0, 2] = R[2, 0] - R[0, 2]
    F[0, 3] = R[0, 1] - R[1, 0]

    F[1, 0] = R[1, 2] - R[2, 1]
    F[1, 1] = R[0, 0] - R[1, 1] - R[2, 2]
    F[1, 2] = R[0, 1] + R[1, 0]
    F[1, 3] = R[0, 2] + R[2, 0]

    F[2, 0] = R[2, 0] - R[0, 2]
    F[2, 1] = R[0, 1] + R[1, 0]
    F[2, 2] = -R[0, 0] + R[1, 1] - R[2, 2]
    F[2, 3] = R[1, 2] + R[2, 1]

    F[3, 0] = R[0, 1] - R[1, 0]
    F[3, 1] = R[0, 2] + R[2, 0]
    F[3, 2] = R[1, 2] + R[2, 1]
    F[3, 3] = -R[0, 0] - R[1, 1] + R[2, 2]

    # diagonalize it
    l, U = np.linalg.eig(F)

    # extract the eigenvect of the highest eigenvalues
    indmax = np.argmax(l)
    q0, q1, q2, q3 = U[:, indmax]

    # form the rotation matrix (eq. 33 ref[1])
    U = np.zeros((3, 3))

    U[0, 0] = q0**2 + q1**2 - q2**2 - q3**2
    U[0, 1] = 2 * (q1 * q2 - q0 * q3)
    U[0, 2] = 2 * (q1 * q3 + q0 * q2)
    U[1, 0] = 2 * (q1 * q2 + q0 * q3)
    U[1, 1] = q0**2 - q1**2 + q2**2 - q3**2
    U[1, 2] = 2 * (q2 * q3 - q0 * q1)
    U[2, 0] = 2 * (q1 * q3 - q0 * q2)
    U[2, 1] = 2 * (q2 * q3 + q0 * q1)
    U[2, 2] = q0**2 - q1**2 - q2**2 + q3**2

    return U
