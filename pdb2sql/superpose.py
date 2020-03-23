
import os
import numpy as np
from .pdb2sqlcore import pdb2sql
from .transform import rotate


def superpose(mobile, target, method='svd', only_backbone=True, export = True, **kwargs):
    """superpose two complexes

    Arguments:
        mobile {str or pdb2sql} -- name or sqldb of the mobile pdb
        target {str or pdb2sql} -- name or sqldb of the target pdb

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

    if not isinstance(mobile, pdb2sql):
        sql_mobile = pdb2sql(mobile)
    else:
        sql_mobile = mobile

    if not isinstance(target, pdb2sql):
        sql_target = pdb2sql(target)
    else:
        sql_target = target

    if only_backbone:
        if 'name' not in kwargs:
            kwargs['name'] = backbone_atoms
        else:
            raise ValueError('Atom type specified but only_backbone == True')

    # selections of some atoms
    selection_mobile = np.array(sql_mobile.get("x,y,z", **kwargs))
    selection_target = np.array(sql_target.get("x,y,z", **kwargs))

    # superpose the complexes
    sql_mobile = superpose_selection(sql_mobile,sql_target,
                                     selection_mobile,selection_target,
                                     method)
    # export a pdb file
    if export:
        target_name = os.path.basename(sql_target.pdbfile).rstip('pdb')
        mobile_name = os.path.basename(sql_mobile.pdbfile).rstip('pdb')
        fname = mobile_name + '_superposed_on_' + \
            target_name + '.pdb'
        sql_mobile.exportpdb(fname)

    return sql_mobile

def superpose_selection(sql_mobile, sql_target, 
                        selection_mobile, selection_target,
                        method):
    """superpose mobile on target using the selected atoms in the selection
    
    Arguments:
        sql_mobile {pdb2sql} -- sqldb of the mobile pdb
        sql_target {pdb2sql} -- sqldb of the target pdb
        selection_mobile {np.ndarray} -- coordinates of the atoms 
                                         in the mobile used to define 
                                         the rotation
        selection_target {np.ndarray} -- coordinates of the atoms 
                                         in the target used to define 
                                         the rotation
    
    Keyword Arguments:
        method {str} -- method used to define the rotation svd or quaternion 
                        (default: {'svd'})
    """

    # translation vector
    tr_mobile = get_trans_vect(selection_mobile)
    tr_target = get_trans_vect(selection_target)

    # rotation matrix
    selection_target += tr_target
    selection_mobile += tr_mobile
    center = np.mean(selection_mobile, 0)
    rmat = get_rotation_matrix(selection_mobile, selection_target, method=method)

    # transform the coordinate of second pdb
    xyz_mobile = np.array(sql_mobile.get("x,y,z"))
    xyz_mobile += tr_mobile
    xyz_mobile = rotate(xyz_mobile, rmat, center=center)
    xyz_mobile -= tr_target

    # update the second sql
    sql_mobile.update('x,y,z', xyz_mobile)

    return sql_mobile

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
