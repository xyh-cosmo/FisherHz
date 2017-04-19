#Base class for Hz Cosmo

import numpy as np
import scipy as sp
import matplotlib.pylab as plt

class HzCosmo:
	"""
	Simple cosmology module.  Interplolate over {zi,Ei} to compute luminosity distance.
	"""
	def __init__(self,z_interp=[],E_fid=[],zmax=1.5,z_size=10,OmegaM=0.3,OmegaK=0.0,H0=70.0,cosmo_array_size=100):
		self.clight = 2.9979e5	# speed of light, in unit of km/s
		self.rMB =  # reduce absolute magnitude in B-band
		self.Om	= OmegaM
		self.Ok	= OmegaK
		self.Ol = 1.0 - self.Om - self.Ok
		self.H0	= H0
		self.interplnE = None
		self.array_size = cosmo_array_size
		self.rhoTot = np.zeros(self.array_size)
		self.dL_z = np.zeros(self.array_size)
		self.dL = np.zeros(self.array_size)

		if len(z_interp) == z_size:
			self.z_interp = z_interp
			if len(E_fid) == z_size:
				self.E_bin = E_fid
				self.E_fid = E_fid
		else:
			def Efun(z):
				return (self.Om*(1+z)**3 + self.Ol)**0.5
			ztemp = np.linspace(0.,zmax,1+z_size)
			self.z_interp = ztemp
			self.E_fid = Efun(ztemp)
			self.E_interp = self.E_fid
			self.interplnE= interp1d(self.z_interp,log(self.E_interp))

	def UpdateEfun(self,E_idx,dE=0.01):
		if E_idx >=1 and E_idx <= len(self.z_interp):
			print '--> updating %2d -th Ei'
			self.E_interp[E_idx] += dE
			self.interplnE= interp1d(self.z_interp,log(self.E_interp))
			self.E_interp[E_idx] += self.E_fid[E_idx] # reset to the fiducial value

	def UpdateHzCosmo(self):


	def Efun(self,z):
		return exp(self.interplnE(z))

	def run_test(self):
		plt.plot(self.z_interp,self.E_fid,'-ro',label=r'test E(z)')
		plt.show()

	def update_H(self,


# TEST part

HC = HzCosmo()

HC.run_test()
