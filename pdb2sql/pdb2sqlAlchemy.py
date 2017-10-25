from sqlalchemy import *
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.schema import Sequence
from sqlalchemy.orm import sessionmaker

from time import time
import re 

Base = declarative_base()


class ATOM(Base):

	__tablename__ = 'ATOM'
	rowID = Column(Integer,primary_key=True)
	serial = Column(Integer,nullable=False)
	name = Column(String(5),nullable=False)
	altLoc = Column(String(5),nullable=False)
	resName = Column(String(5),nullable=False)
	chainID = Column(String(5),nullable=False)
	resSeq = Column(Integer,nullable=False)
	iCode = Column(String(5),nullable=False)
	x = Column(Float,nullable=False)
	y = Column(Float,nullable=False)
	z = Column(Float,nullable=False)
	occ = Column(Float,nullable=False)
	temp = Column(Float,nullable=False)


class pdb2sql_alchemy(object):

	def __init__(self,pdbfile,sqlfile=None,fix_chainID=False,verbose=False):

		self.pdbfile = pdbfile
		self.sqlfile = sqlfile
		self.verbose = verbose

		self._create_sql()

	############################################################################################
	#
	#	CREATE THE DB
	#
	############################################################################################

	def _create_sql(self):

		pdbfile = self.pdbfile
		sqlfile = self.sqlfile


		# column names and types
		self.col = {'serial' : 'INT',
		       'name'   : 'TEXT',
		       'altLoc' : 'TEXT',
		       'resName' : 'TEXT',
		       'chainID' : 'TEXT',
			   'resSeq'  : 'INT',
			   'iCode'   : 'TEXT',
			   'x'       : 'REAL',
			   'y'       : 'REAL',
			   'z'       : 'REAL',
			   'occ'     : 'REAL',
			   'temp'    : 'REAL'}

	    # delimtier of the column format
	    # taken from http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
		self.delimiter = {
					'serial' : [6,11],
					'name'   : [12,16],
					'altLoc' : [16,17],
					'resName' :[17,20],
					'chainID' :[21,22],
					'resSeq'  :[22,26],
					'iCode'   :[26,26],
					'x'       :[30,38],
					'y'       :[38,46],
					'z'       :[46,54],
					'occ'     :[54,60],
					'temp'    :[60,66]}

		# sqlalchemy init
		self.engine = create_engine('sqlite:///:memory:')
		Base.metadata.create_all(self.engine)
		DBSession = sessionmaker(bind=self.engine)
		self.session = DBSession()


		# read the pdb file a pure python way
		# RMK we go through the data twice here. Once to read the ATOM line and once to parse the data ... 
		# we could do better than that. But the most time consuming step seems to be the CREATE TABLE query
		with open(pdbfile,'r') as fi:
			data = [line.split('\n')[0] for line in fi if line.startswith('ATOM')]


		# old format chain ID fix
		del_copy = self.delimiter.copy()
		if data[0][del_copy['chainID'][0]] == ' ':
			del_copy['chainID'] = [72,73]

		for iatom,atom in enumerate(data):

			# sometimes we still have an empty line somewhere
			if len(atom) == 0:
				continue

			# browse all attribute of each atom
			at = {}
			for ik,(colname,coltype) in enumerate(self.col.items()):

				# get the piece of data
				data = atom[del_copy[colname][0]:del_copy[colname][1]].strip()

				# convert it if necessary
				if coltype == 'INT':
					data = int(data)
				elif coltype == 'REAL':
					data = float(data)

				# append to dict
				at[colname] = data

			# create a new ATOM
			newat = ATOM(serial=at['serial'],name=at['name'],altLoc=at['altLoc'],resName=at['resName'],
				         chainID=at['chainID'],resSeq=at['resSeq'],iCode=at['iCode'],
				         x=at['x'],y=at['y'],z=at['z'],occ=at['occ'],temp=at['temp'])

			# add the atom to the data base
			self.session.add(newat)

		self.session.commit()


	############################################################################################
	#
	#		GET FUNCTIONS
	#
	#			get(self,attribute,selection) -> return the atribute(s) value(s) for the given selection 
	#				
	#			attribute : example ['x','y,'z'] --> return the x y z in an array of tuple
	#						        'x,y,z'		 --> return the x y z in an array of tuple
	#                               'resSeq'     --> return the resSeq in an array
	#                                None        --> return a list of ATOM objects
	#		
	#			**kwargs  : example name = ['CA','O'], return only the atom whose name are CA or O
	#                               chainID = 'A'      return only the atom of chainID A
	#                               no_name = ['H']    return all atoms except the H
	#
	#		RMK : The rowID starts at 0 in the SLQ database. The rowID are etherefore converted 
	#             to start at 0 both in the attribute and the conditions
	###############################################################################################

	def get(self,attribute=None,**kwargs):


		# parse the commas in the attribute
		if attribute is not None:
			if ',' in attribute:
				attribute = attribute.split(',')


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
				if isinstance(attribute,list):
					return [  tuple( atom.__dict__[at]-1 if at=='rowID' else atom.__dict__[at] for at in attribute  ) for atom in self.session.query(ATOM) ]

				# or as a list if we have onlye one attribute
				else:
					return [ atom.__dict__[attribute]-1 if attribute=='rowID' else atom.__dict__[attribute] for atom in self.session.query(ATOM) ]

		# if we have a condition
		else:

			# create the atom list
			atom_list = self.session.query(ATOM)

			# loop through the kwargs
			for key,values in kwargs.items():

				# make sure the values are array for extraction
				if not isinstance(values,list):
					values = [values]

				# fix indexing python v SQL
				# 0 <--> 1
				if 'rowID' in key:
					values = [v+1 for v in values]

				#  filter the atom list
				if 'no_' in key:
					atom_list = atom_list.filter(~ATOM.__dict__[re.sub('no_','',key)].in_(values))
				else:
					atom_list = atom_list.filter(ATOM.__dict__[key].in_(values))

			# if no attribute specfied we return an array of ATOM objet
			if attribute is None:
				return atom_list.all()

			else:

				# or we extract data
				# as a list of tuples if we have several attributes
				if isinstance(attribute,list):
					return [  tuple( atom.__dict__[at]-1 if at=='rowID' else atom.__dict__[at]  for at in attribute  ) for atom in atom_list ]

				# or as a list if we have onlye one attribute
				else:
					return [ atom.__dict__[attribute]-1 if attribute=='rowID' else atom.__dict__[attribute]   for atom in atom_list ]

	############################################################################################
	#
	#		UPDATE FUNCTIONS
	#
	#			update(self,attribute,values,selection) -> update the atribute(s) value(s) for the given selection 
	#				
	#			attribute : example ['x','y,'z'] --> update the x y z in an array of tuple
	#						        'x,y,z'		 --> update the x y z in an array of tuple
	#                               'resSeq'     --> update the resSeq in an array
	# 
	#			values    : an array of values that corresponds to the number of 
	#                       attributes and number of atom seleted 
	#                       format : nATOM x nAttribute 
 	#		
	#			**kwargs  : example name = ['CA','O'], return only the atom whose name are CA or O
	#                               chainID = 'A'      return only the atom of chainID A
	#                               no_name = ['H']    return all atoms except the H
	#
	#
	#
	#		RMK : The rowID starts at 0 in the SLQ database. The rowID are etherefore converted 
	#             to start at 0 both in the attribute and the conditions
	###############################################################################################

	def update(self,attribute,values,**kwargs):

		# parse the commas in the attribute
		if ',' in attribute:
			attribute = attribute.split(',')

		if not isinstance(attribute,list):
			attribute = [attribute]

		# size of the values		
		nrow = len(values)
		ncol = len(values[0])
	

		nat = len(attribute)
		if nat != ncol:
			raise ValueError('Values and Attribute have incompatible sizes',nat,ncol)

		# no selection
		if len(kwargs) == 0:

			# check size
			if nrow != len(self.session.query(ATOM).all()):
				raise ValueError('Wrong number of values for the ATOM selection')

			# goes through all the ros
			for irow in range(nvalues):

				# create  a dict of values
				dict_values = {}
				for icol,at in enumerate(attribute):
					dict_values[at] = values[irow][icol]

				# update
				self.session.query(ATOM).filter(ATOM.rowID == irow).update(dict_values)

		# if there is a selection
		else:

			# create the atom list
			atom_list = self.session.query(ATOM)

			# loop through the kwargs
			for key,val in kwargs.items():

				# make sure the values are array for extraction
				if not isinstance(val,list):
					val = [val]

				# fix indexing python v SQL
				# 0 <--> 1
				if 'rowID' in key:
					val = [v+1 for v in val]

				#  filter the atom list
				if 'no_' in key:
					atom_list = atom_list.filter(~ATOM.__dict__[re.sub('no_','',key)].in_(val))
				else:
					atom_list = atom_list.filter(ATOM.__dict__[key].in_(val))

			# check size
			if nrow != len(atom_list.all()):
				raise ValueError('Wrong number of values for the ATOM selection')

			# get the indexes
			indexes = self.get('rowID',**kwargs)

			#here the conversion to 0 starting index is annoying
			indexes = [i+1 for i in indexes]

			# goes through all the ros
			for ival,ind in enumerate(indexes):

				# create  a dict of values
				dict_values = {}
				for icol,at in enumerate(attribute):
					dict_values[at] = values[ival][icol]

				# update
				self.session.query(ATOM).filter(ATOM.rowID == ind).update(dict_values)
			
		self.session.commit()

	###################################################################################################
	#
	#	Export a pdb file	
	#
	###################################################################################################

	def exportpdb(self,fname,**kwargs):

		'''
		Export a PDB file with kwargs selection
		not pretty so far but functional
		'''

		# get the data
		data = self.get('*',**kwargs)

		# write each line
		# the PDB format is pretty strict
		# http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
		f = open(fname,'w')
		for d in data:
			line = 'ATOM  '
			line += '{:>5}'.format(d[0])	# serial
			line += ' '
			line += '{:^4}'.format(d[1])	# name
			line += '{:>1}'.format(d[2])	# altLoc
			line += '{:>3}'.format(d[3])	#resname
			line += ' '
			line += '{:>1}'.format(d[4])	# chainID
			line += '{:>4}'.format(d[5])	# resSeq
			line += '{:>1}'.format(d[6])	# iCODE
			line += '   '
			line += '{: 8.3f}'.format(d[7])	#x
			line += '{: 8.3f}'.format(d[8])	#y
			line += '{: 8.3f}'.format(d[9])	#z
			line += '{: 6.2f}'.format(d[10])	# occ
			line += '{: 6.2f}'.format(d[11])	# temp
			line += '\n'

			f.write(line)

		# close
		f.close()




if __name__ == "__main__":

	import numpy as np

	# create the data base
	t0 = time()
	db = pdb2sql_alchemy('5hvd.pdb')
	print('ALCH %f' %(time()-t0))

	# extract the xyz position of all VAL and LEU resiues of chain A but not the H atoms
	xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],no_name=['H'])

	# put the data back 
	db.update('x,y,z',xyz,chainID='A',resName=['VAL','LEU'],no_name=['H'])