import re
import os
import urllib.request
from urllib.error import HTTPError

def fetch(pdbid, outdir='.'):
    """Download PDB file from PDB website https://www.rcsb.org/.

    Args:
        pdbid(str): PDB ID
        outdir (str, optional): Output path. Defaults to '.'.

    Raises:
        ValueError: PDB ID is not valid
        ValueError: PDB ID is valid but does not exist on PDB website

    Examples:
        >>> from pdb2sql import fetch
        >>> fetch('1cbh')
    """
    # defaults
    hosturl = 'http://files.rcsb.org/download'

    # check pdbid
    p = re.compile('[0-9a-z]{4,4}$', re.IGNORECASE)
    if not p.match(pdbid):
        raise ValueError(f'Invalid PDB ID: {pdbid}.')
    pdb = pdbid + '.pdb'
    cif = pdbid + '.cif'

    # build downloading url
    url_pdb = os.path.join(hosturl, pdb)
    url_cif = os.path.join(hosturl, cif)
    fout = os.path.join(outdir, pdb)

    # get url content
    try:
        pdbdata = urllib.request.urlopen(url_pdb)
    except HTTPError:
        try:
            cifdata = urllib.request.urlopen(url_cif)
            raise ValueError(f'The PDB ID given is only represented in '
                f'mmCIF format and pdb2sql does not handle mmCIF format.')
        except HTTPError:
            raise ValueError(f'PDB ID not exist: {pdbid}')

    # write to file
    with open(fout, 'wb') as f:
        f.write(pdbdata.read())