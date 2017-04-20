#Base class for Hz Cosmo
import sys,os
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
# import scipy as sp
import matplotlib.pylab as plt

class HzCosmo:
	"""
	Simple cosmology module.  Interplolate over {zi,Ei} to compute luminosity distance.
	"""
	def __init__(self,z_interp=[],zmax=1.5,OmegaM=0.3,OmegaK=0.0,H0=70.0,cosmo_array_size=100,dL_zmin=0.005):
		self.clight = 2.9979e5	# speed of light, in unit of km/s
		self.Om	= OmegaM
		self.Ok	= OmegaK
		self.Ol = 1.0 - self.Om - self.Ok
		self.H0	= H0
		self.MB = -19.3
		self.rMB = 5.*np.log10(self.clight/self.H0) + 25.0 + self.MB # reduce absolute magnitude in B-band

		self.interplnE = None
		self.interldL = None
		self.array_size = cosmo_array_size
		self.dL_z = np.linspace(dL_zmin,zmax,self.array_size)
		self.dL = np.zeros(self.array_size)
		self.mB = np.zeros(self.array_size)
		self.dmBdE = []  # store derivatives of mB wrt to Ei

		if len(z_interp) < 5:
			print('length of z_interp is too small, pls increase it!')
			sys.exit(0)

		if max(z_interp) != zmax:
			print('max(z_interp) does not match zmax = %g'%zmax)
			sys.exit(0)

		z_size = len(z_interp)
		if min(z_interp) == 0.0:
			z_size += 1
			self.z_interp = np.array(z_interp)
		else:
			ztemp = [0]
			for i in range(len(z_interp)):
				ztemp.append(z_interp[i])
			print 'ztemp = ',ztemp
			self.z_interp = np.array(ztemp)

		def Efun(z):
			return (self.Om*(1+z)**3 + self.Ol)**0.5
		self.E_fid = Efun(self.z_interp)
		self.E_interp = self.E_fid
		self.interplnE= interp1d(self.z_interp,np.log(self.E_interp))

		self.UpdateHzCosmo()
		self.PrintParams()

	def PrintParams(self):
		print('OmegaM = %g'%self.Om)
		print('OmegaK = %g'%self.Ok)
		print('OmegaL = %g'%self.Ol)
		print('H0     = %g'%self.H0)
		print('MB     = %g'%self.MB)
		print('rMB    = %g'%self.rMB)

	def UpdateEfun(self,E_idx,dE=0.01):
		if E_idx >=1 and E_idx <= len(self.z_interp):
			print '--> updating %2d -th Ei'
			self.E_interp[E_idx] += dE
			self.interplnE= interp1d(self.z_interp,log(self.E_interp))
			self.E_interp[E_idx] = self.E_fid[E_idx] # reset to the fiducial value

	def UpdateHzCosmo(self):
		def Einv(z):
			return 1./self.Efun(z)

		for i in range(len(self.dL_z)):
			self.dL[i] = quad(Einv,0.,self.dL_z[i])[0]
			self.mB[i] = 5.0*np.log10((1.0+self.dL_z[i])*self.dL[i]) + self.rMB


	def Efun(self,z):
		return np.exp(self.interplnE(z))

	def run_test(self):
		plt.subplot(1,2,1)
		plt.plot(self.z_interp,self.E_fid,'-ro',label=r'test E(z)')
		plt.legend()

		plt.subplot(1,2,2)
		plt.plot(self.dL_z,self.mB,label=r'apparent magnitude')
		plt.legend()

		plt.show()

	# def update_H(self,


# TEST part

z_interp = [0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5]
HC = HzCosmo(z_interp=z_interp)

HC.run_test()
