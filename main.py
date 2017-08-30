#import matpotlib as plt
import numpy as np
import pandas as pd
import os
import sys


# important things to note: positive residues will start
# from ph - pk and gradually ascend to 0 (when the x axis,
# ph value reaches that of pk. It will stay that way for the
# rest of the ph values on the x axis. For negative residues,
# the energy will start at 0, and then stay that way until ph
# (the x axis) matches the value of pk. After this event, the
# energy will decrease constantly until it reaches its final
# value of pk - ph. (last ph, of course, is 14)

class pk_shifts(object):
	def __init__(self):
		self.pksol_dict = dict(ARG = 12.5, HIS = 6.98, LYS = 10.4,
        		        ASP = 4.75, GLU = 4.75, TYR = 10.2)
		self.res_charges = dict(POSITIVE = ['ARG', 'HIS', 'LYS'],
				NEGATIVE = ['ASP', 'GLU', 'TYR'])
		self.x_axis = np.linspace(0,14,15)	

	def gather_data(self, pk_filepath):
        	if not(os.path.exists(pk_filepath)):
			sys.exit('SPECIFIED PATH FOR pK.out NOT FOUND')
	       	raw_data = pd.read_csv(pk_filepath, sep = '\s+', header = None,
					usecols = [0,1], names = ['Res','pK'])
	        raw_data.set_index('Res', inplace = True)

		# using pd to filter out residues for which we have default
		# pkSol values. These residues are in the self.res_charges
		# dictionary.
		self.data = raw_data[ [raw_data.index.values[i][0:3] 
				in np.ravel(self.res_charges.values()) for i in 
				range(len(raw_data.index))] ]
		# filter out amino acids that have pK values that are >14 
		# or that are <0. 
		self.data = self.data[ [self.data['pK'].values[i] not in 
				['<0.0','>14.0'] for i in 
				range(len(self.data['pK']))] ] 

		self.positive_res = self.data[ [self.data.index.values[i][0:3]
				in self.res_charges['POSITIVE'] for i in 
				range(len(self.data.index))] ]

		self.negative_res = self.data[ [self.data.index.values[i][0:3]
				in self.res_charges['NEGATIVE'] for i in 
				range(len(self.data.index))] ]

	def split_by_charge(self):
		'''Reads in the data and splits into two dataframes, one 
		including only positive residues, and one including only 
		negative residues'''

						
		
		
if(__name__ == '__main__'):
	p = pk_shifts()
	p.gather_data('pK.out')		
	
