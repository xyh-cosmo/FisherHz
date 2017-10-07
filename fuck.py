import os,sys
sys.path.append('python')
from pylab import *
from HzCosmoBase import *

#mock_sne = loadtxt('wfirst_LCDM.txt')
mock_sne = loadtxt('wfirst_LCDM_int_0.09_0.0_meas_0.08_lens_0.07_sys_0.02.txt')
# mock_sne = loadtxt('wfirst_LCDM_sigma_int_0.1.txt')

sne_z = mock_sne[:,0]
sne_err = mock_sne[:,2]
sne_err2 = sne_err**2
#sne_err2 = zeros(len(sne_err))+0.1**2

zmin = 0.125
zmax = 1.7

z_interp = linspace(zmin,zmax,10)
HC0 = HzCosmo(z_interp=z_interp,cosmo_array_size=200,zmax=1.7)

def get_cov_H(z_interp):
    HC = HzCosmo(z_interp=z_interp,cosmo_array_size=200,zmax=1.7,kind='linear')
    HC.ComputedmB()

    dim = len(z_interp)+1
    fish = zeros((dim,dim))

    for i in range(dim-1):
        for j in range(dim-1):
            fish[i,j] = sum(HC.dmBdE[i](sne_z)*HC.dmBdE[j](sne_z)/sne_err2)

    for i in range(dim-1):
        fish[dim-1,i] = sum(HC.dmBdE[i](sne_z)/sne_err2)
        fish[i,dim-1] = fish[dim-1,i]

    fish[dim-1,dim-1] = sum(1./sne_err2)

    cov_H_MB = inv(fish)
    cov_H = cov_H_MB[0:dim-1,0:dim-1]
    return cov_H

z6 = array([0.17 ,  0.476,  0.782,  1.088,  1.394,  1.7])
z8 = array([0.17,  0.38857143,  0.60714286,  0.82571429,  1.04428571, 1.26285714,  1.48142857,  1.7])
z10 = array([0.17,  0.34,  0.51,  0.68,  0.85,  1.02,  1.19,  1.36, 1.53,  1.7])

zz = []
zz.append(z6)
zz.append(z8)
zz.append(z10)

for i in range(len(zz)):
    cov_H = get_cov_H(zz[i])
#    plot(z,HC0.Efun(z),'-o',label=r'# zbin = '+str(n))
    plot(zz[i],diag(cov_H)**0.5/HC0.Efun(zz[i]),'-o',label=r'# zbin = '+str(len(zz[i])))
#    plot(z,diag(cov_H)**0.5/HC0.Efun(z)*(20./n)**0.5,'-o',label=r'# zbin = '+str(n))
#    plot(z,diag(cov_H)**0.5/HC0.Efun(z)*(20./n)**2,'-o',label=r'# zbin = '+str(n))
#    plot(z,diag(cov_H)**0.5*(20./n)**0,'-o',label=r'# zbin = '+str(n))

#ylim(0,0.5)
legend()
show()    

