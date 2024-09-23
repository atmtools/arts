#! /usr/bin/env python3

import numpy as np


class SphereDataPhysical:

	def __init__(self, filename = None):

		if filename != None:
			self.read_file(filename)

		pass

	def read_file(self, filename):
		"""
		Load sphere data stored in physical space from .csv file

		Return:
		-------
			Tuple with (data, longitude angles, latitude angles)
		"""
		print("Loading file: "+filename)

		try:
			data = np.loadtxt(filename, skiprows=0)
		except Exception as e:
			raise e

		# First row and col are longitude and latitude coordinates
		self.labelsx = data[0,0:]
		self.labelsy = data[0:,0]
		self.data = data[1:,1:]

