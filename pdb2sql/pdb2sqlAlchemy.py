from sqlalchemy import *
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.schema import Sequence
from sqlalchemy.orm import sessionmaker
from .pdb2sql_base import pdb2sql_base
from time import time
import re

Base = declarative_base()


class ATOM(Base):
    '''SQLAlchemy object for atoms.'''

    __tablename__ = 'ATOM'
    rowID = Column(Integer, primary_key=True)
    serial = Column(Integer, nullable=False)
    name = Column(String(6), nullable=False)
    altLoc = Column(String(1), nullable=False)
    resName = Column(String(3), nullable=False)
    chainID = Column(String(1), nullable=False)
    resSeq = Column(Integer, nullable=False)
    iCode = Column(String(1), nullable=False)
    x = Column(Float, nullable=False)
    y = Column(Float, nullable=False)
    z = Column(Float, nullable=False)
    occ = Column(Float, nullable=False)
    temp = Column(Float, nullable=False)
    element = Column(String(2), nullable=False)
    model = Column(Integer, nullable=False)


class pdb2sql_alchemy(pdb2sql_base):

    def __init__(
            self,
            pdbfile,
            sqlfile=None,
            fix_chainID=False,
            verbose=False):
        '''Use sqlAlchemy to load the database.'''
        super().__init__(pdbfile, sqlfile, fix_chainID, verbose)
        self._create_sql()

    def _create_sql(self):

        # sqlalchemy init
        self.engine = create_engine('sqlite:///:memory:')
        Base.metadata.create_all(self.engine)
        DBSession = sessionmaker(bind=self.engine)
        self.session = DBSession()

        # read the pdb file a pure python way
        # RMK we go through the data twice here. Once to read the ATOM line and once to parse the data ...
        # we could do better than that. But the most time consuming step seems to be the CREATE TABLE query
        # with open(pdbfile,'r') as fi:
        #   data = [line.split('\n')[0] for line in fi if line.startswith('ATOM')]

        fi = open(self.pdbfile, 'r')
        self.nModel = 0
        _check_format = True

        for line in fi:

            if line.startswith('ATOM'):
                line = line.split('\n')[0]

            elif line.startswith('ENDMDL'):
                self.nModel += 1
                continue

            else:
                continue

            if _check_format:
                # old format chain ID fix
                del_copy = self.delimiter.copy()
                if line[del_copy['chainID'][0]] == ' ':
                    del_copy['chainID'] = [72, 73]
                if line[del_copy['chainID'][0]] == ' ':
                    raise ValueError('chainID not found sorry')
                _check_format_ = False

            # browse all attribute of each atom
            at = {}
            for ik, (colname, coltype) in enumerate(self.col.items()):

                # get the piece of data
                if colname in del_copy.keys():
                    data = line[del_copy[colname][0]:del_copy[colname][1]].strip()

                    # convert it if necessary
                    if coltype == 'INT':
                        data = int(data)
                    elif coltype == 'REAL':
                        data = float(data)


                    # append to dict
                    at[colname] = data

            # create a new ATOM
            newat = ATOM(
                serial=at['serial'],
                name=at['name'],
                altLoc=at['altLoc'],
                resName=at['resName'],
                chainID=at['chainID'],
                resSeq=at['resSeq'],
                iCode=at['iCode'],
                x=at['x'],
                y=at['y'],
                z=at['z'],
                occ=at['occ'],
                temp=at['temp'],
                element=at['element'],
                model=self.nModel)

            # add the atom to the data base
            self.session.add(newat)

        self.session.commit()
        fi.close()

    def get(self, attribute=None, **kwargs):
        '''Exectute a simple SQL query that extracts values of attributes for certain condition.

        Args :
            attribute (str) : attribute to retreive eg : ['x','y,'z'], 'xyz', 'resSeq'
                              if None all the attributes are returned

            **kwargs : argument to select atoms eg : name = ['CA','O'], chainID = 'A', no_name = ['H']

        Returns:
            data : array containing the value of the attributes

        Example :
        >>> db.get('x,y,z',chainID='A',no_name=['H']")
        '''

        # parse the commas in the attribute
        if attribute is not None:
            if ',' in attribute:
                attribute = attribute.split(',')

        if 'model' not in kwargs.keys() and self.nModel > 0:
            model_data = []
            for iModel in range(self.nModel):
                kwargs['model'] = iModel
                model_data.append(self.get(attribute, **kwargs))
            return model_data

        # no selection specified
        if len(kwargs) == 0:

            # no attribute specfied
            # we return all the ATOM objects
            if attribute is None:
                return self.session.query(ATOM).all()

            # if we specify attributes
            else:

                # extract data
                # as a list of tuples if we have several attributes
                if isinstance(attribute, list):
                    return [tuple(atom.__dict__[at] - 1 if at == 'rowID' else atom.__dict__[at]
                                  for at in attribute) for atom in self.session.query(ATOM)]

                # or as a list if we have onlye one attribute
                else:
                    return [
                        atom.__dict__[attribute] -
                        1 if attribute == 'rowID' else atom.__dict__[attribute] for atom in self.session.query(ATOM)]

        # if we have a condition
        else:

            # create the atom list
            atom_list = self.session.query(ATOM)

            # loop through the kwargs
            for key, values in kwargs.items():

                # make sure the values are array for extraction
                if not isinstance(values, list):
                    values = [values]

                # fix indexing python v SQL
                # 0 <--> 1
                if 'rowID' in key:
                    values = [v + 1 for v in values]

                #  filter the atom list
                if 'no_' in key:
                    atom_list = atom_list.filter(
                        ~ATOM.__dict__[re.sub('no_', '', key)].in_(values))
                else:
                    atom_list = atom_list.filter(
                        ATOM.__dict__[key].in_(values))

            # if no attribute specfied we return an array of ATOM objet
            if attribute is None:
                return atom_list.all()

            else:

                # or we extract data
                # as a list of tuples if we have several attributes
                if isinstance(attribute, list):
                    return [
                        tuple(
                            atom.__dict__[at] -
                            1 if at == 'rowID' else atom.__dict__[at] for at in attribute) for atom in atom_list]

                # or as a list if we have onlye one attribute
                else:
                    return [
                        atom.__dict__[attribute] -
                        1 if attribute == 'rowID' else atom.__dict__[attribute] for atom in atom_list]

    def update(self, attribute, values, **kwargs):
        '''Update the database.

        Args:
            attribute (str) : string of attribute names eg. ['x','y,'z'], 'xyz', 'resSeq'

            values (np.ndarray) : an array of values that corresponds
                                  to the number of attributes and atoms selected

            **kwargs : selection arguments eg : name = ['CA','O'], chainID = 'A', no_name = ['H']
        '''

        # parse the commas in the attribute
        if ',' in attribute:
            attribute = attribute.split(',')

        if not isinstance(attribute, list):
            attribute = [attribute]

        # size of the values
        nrow = len(values)
        ncol = len(values[0])

        nat = len(attribute)
        if nat != ncol:
            raise ValueError(
                'Values and Attribute have incompatible sizes', nat, ncol)

        # handle the multi model cases
        if 'model' not in kwargs.keys() and self.nModel > 0:
            for iModel in range(self.nModel):
                kwargs['model'] = iModel
                self.update(attribute, values, **kwargs)
            return

        # no selection
        if len(kwargs) == 0:

            # check size
            if nrow != len(self.session.query(ATOM).all()):
                raise ValueError(
                    'Wrong number of values for the ATOM selection')

            # goes through all the ros
            for irow in range(nrow):

                # create  a dict of values
                dict_values = {}
                for icol, at in enumerate(attribute):
                    dict_values[at] = values[irow][icol]

                # update
                self.session.query(ATOM).filter(
                    ATOM.rowID == irow).update(dict_values)

        # if there is a selection
        else:

            # create the atom list
            atom_list = self.session.query(ATOM)

            # loop through the kwargs
            for key, val in kwargs.items():

                # make sure the values are array for extraction
                if not isinstance(val, list):
                    val = [val]

                # fix indexing python v SQL
                # 0 <--> 1
                if 'rowID' in key:
                    val = [v + 1 for v in val]

                #  filter the atom list
                if 'no_' in key:
                    atom_list = atom_list.filter(
                        ~ATOM.__dict__[re.sub('no_', '', key)].in_(val))
                else:
                    atom_list = atom_list.filter(ATOM.__dict__[key].in_(val))

            # check size
            if nrow != len(atom_list.all()):
                raise ValueError(
                    'Wrong number of values for the ATOM selection')

            # get the indexes
            indexes = self.get('rowID', **kwargs)

            # here the conversion to 0 starting index is annoying
            indexes = [i + 1 for i in indexes]

            # goes through all the ros
            for ival, ind in enumerate(indexes):

                # create  a dict of values
                dict_values = {}
                for icol, at in enumerate(attribute):
                    dict_values[at] = values[ival][icol]

                # update
                self.session.query(ATOM).filter(
                    ATOM.rowID == ind).update(dict_values)

        self.session.commit()
