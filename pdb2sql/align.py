import numpy as np
from .pdb2sqlcore import pdb2sql
from .interface import interface
from .transform import rot_xyz_around_axis


def align(pdb, axis=None, export=True, **kwargs):
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

    # get the pca eigenvect we want to align
    vect = get_max_pca_vect(xyz)

    # align the sql
    sql = align_pca_vect(sql, vect, axis)

    # export the pdbfile
    if export:
        export_aligned(sql)

    return sql


def align_interface(ppi, plane='xy', export=True, **kwargs):
    """align the interface of a complex in a given plane

    Arguments:
        ppi {interface} -- sql interface or pdb file
        plane {str} -- plane for alignement

    Keyword Arguments:
        export {bool} -- write a pdb file (default: {True})
        kwargs {dict} -- keywaord argument from interface.get_contact_atoms method
    """

    if not isinstance(ppi, interface):
        sql = interface(ppi)
    else:
        sql = ppi

    index_contact = sql.get_contact_atoms(**kwargs)
    row_id = []
    for _, v in index_contact.items():
        row_id += v
    xyz = np.array(sql.get('x,y,z', rowID=row_id))

    # get the pca eigenvect we want to align
    vect = get_min_pca_vect(xyz)

    # align the sql database
    dict_plane = {'xy': 'z', 'xz': 'y', 'yz': 'x'}
    sql = align_pca_vect(sql, vect, dict_plane[plane])

    # export the pdbfile
    if export:
        export_aligned(sql)

    return sql


def align_pca_vect(sql, vect, axis):
    """Align the pca vect of the sql along th axis

    Arguments:
        sql {pdb2sql} -- sqldb of the complex
        vect {np.ndarray} -- pca eigenvect
        axis {str} -- axis along which to align vect

    Returns:
        pdb2sql -- aligned sqldb
    """

    # rotation angles
    phi, theta = get_rotation_angle(vect)

    # complete coordinate
    xyz = np.array(sql.get('x,y,z'))

    # align them
    xyz = _align_along_axis(xyz, axis, phi, theta)

    # update the sql
    sql.update('x,y,z', xyz)

    return sql


def export_aligned(sql):
    """export a pdb file of the aligned pdb

    Arguments:
        sql {pdb2sql} -- aligned sqldb
    """
    fname = sql.pdbfile.rstrip('.pdb') + '_aligned.pdb'
    sql.exportpdb(fname)


def get_rotation_angle(vmax):
    """Extracts the rotation angles from the PCA

    Arguments:
        u {np.array} -- eigenvalues of the PCA
        V {np.array} -- eigenvectors of the PCA
    """

    # extract max eigenvector

    x, y, z = vmax
    r = np.linalg.norm(vmax)

    # rotation angle
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)

    return phi, theta


def get_max_pca_vect(xyz):
    """Get the max eigenvector of th pca

    Arguments:
        xyz {numpy.ndarray} -- matrix of the atoms coordinates
    """
    u, v = pca(xyz)
    return v[:, np.argmax(u)]


def get_min_pca_vect(xyz):
    """Get the min eigenvector of th pca

    Arguments:
        xyz {numpy.ndarray} -- matrix of the atoms coordinates
    """
    u, v = pca(xyz)
    return v[:, np.argmin(u)]

def pca(mat):
    """computes the principal component analysis of the points A

    Arguments:
        A {numpy.ndarray} -- matrix of points [npoints x ndim]

    Returns:
        tuple -- eigenvalues, eigenvectors, score
    """
    scat = (mat-np.mean(mat.T, axis=1)).T
    u, v = np.linalg.eig(np.cov(scat))
    return u, v


def _align_along_axis(xyz, axis, phi, theta):
    """align the xyz coordinates along the given axi

    Arguments:
        xyz {numpy.ndarray} -- coordinates of the atoms
        axis {str} -- axis to align
        phi {float} -- azimuthal angle
        theta {float} -- the other angles

    Raises:
        ValueError: axis should be x y or z

    Returns:
        nd.array -- rotated coordinates
    """

    # align along preferred axis
    if axis == 'x':
        xyz = rot_xyz_around_axis(xyz, np.array([0, 0, 1]), -phi)
        xyz = rot_xyz_around_axis(xyz, np.array([0, 1, 0]), np.pi/2 - theta)

    elif axis == 'y':
        xyz = rot_xyz_around_axis(xyz, np.array([0, 0, 1]), np.pi/2 - phi)
        xyz = rot_xyz_around_axis(xyz, np.array([0, 1, 0]), np.pi/2 - theta)

    elif axis == 'z':
        xyz = rot_xyz_around_axis(xyz, np.array([0, 0, 1]), -phi)
        xyz = rot_xyz_around_axis(xyz, np.array([0, 1, 0]), -theta)
    else:
        raise ValueError('axis should be x, y ,or z')

    return xyz
