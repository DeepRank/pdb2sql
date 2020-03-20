import numpy as np
from .pdb2sqlcore import pdb2sql
from .transform import rot_xyz_around_axis


def align(pdb, axis='z', export=True, **kwargs):
    """Align the max principal component of a structure along one of the cartesian axis

    Arguments:
        pdb {str, pdb2sql} -- the pdbfile or the sql database of the complex

    Keyword Arguments:
        axis {str} -- cartesian axis for alignement (default: {'z'})
        export {bool} -- export the aligned structure to file
        **kwargs {dict} -- option to select subpart of the structure for alignement

    Returns:
        pd2sql -- sql databse of the aligned structure

    Example:
        >>> pdb = '1AK4'
        >>> sql = align(pdb,chainID='A')
    """

    if not isinstance(pdb, pdb2sql):
        sql = pdb2sql(pdb)
    else:
        sql = pdb

    # extract coordinate
    xyz = np.array(sql.get('x,y,z', **kwargs))

    # perform pca
    u, v = pca(xyz)

    # extract max eigenvector
    vmax = v[:, np.argmax(u)]
    x, y, z = vmax
    r = np.linalg.norm(vmax)

    # rotation angle
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)

    # complete coordinate
    xyz = np.array(sql.get('x,y,z'))

    # align along preferred axis
    if axis == 'x':
        xyz = rot_xyz_around_axis(xyz, np.array([0, 0, 1]), -phi)
        xyz = rot_xyz_around_axis(xyz, np.array([0, 1, 0]), np.pi/2 - theta)

    if axis == 'y':
        xyz = rot_xyz_around_axis(xyz, np.array([0, 0, 1]), np.pi/2-phi)
        xyz = rot_xyz_around_axis(xyz, np.array([0, 1, 0]), np.pi/2 - theta)

    if axis == 'z':
        xyz = rot_xyz_around_axis(xyz, np.array([0, 0, 1]), -phi)
        xyz = rot_xyz_around_axis(xyz, np.array([0, 1, 0]), -theta)

    # update the sql
    sql.update('x,y,z', xyz)

    # export the pdbfile
    if export:
        fname = sql.pdbfile.rstrip('.pdb') + '_aligned.pdb'
        sql.exportpdb(fname)

    return sql


def pca(A):
    """computes the principal component analysis of the points A

    Arguments:
        A {numpy.ndarray} -- matrix of points [npoints x ndim]

    Returns:
        tuple -- eigenvalues, eigenvectors, score
    """
    scat = (A-np.mean(A.T, axis=1)).T
    u, v = np.linalg.eig(np.cov(scat))
    return u, v
