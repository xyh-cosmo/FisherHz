import os,sys
sys.path.append('python')
from pylab import *
from HzCosmoBase import *
from LogHzCosmoBase import *

#mock_sne = loadtxt('wfirst_LCDM.txt')
mock_sne = loadtxt('wfirst_LCDM_int_0.09_0.0_meas_0.08_lens_0.07_sys_0.02.txt')
# mock_sne = loadtxt('wfirst_LCDM_sigma_int_0.1.txt')

sne_z = mock_sne[:,0]
sne_err = mock_sne[:,2]
sne_err2 = sne_err**2
# sne_err2 = zeros(len(sne_err))+0.05**2

z_interp = linspace(0.05,1.7,10)
HC0 = HzCosmo(z_interp=z_interp,cosmo_array_size=200,zmax=1.7)
# HC0 = HzCosmo(z_interp=z_interp,cosmo_array_size=200,zmax=1.7)

def get_cov_H(z_interp,dE=0.01):
    HC = HzCosmo(z_interp=z_interp,cosmo_array_size=200,zmax=1.7,kind='linear',dE=dE)
    HC.ComputedmB()

    dim = len(z_interp)+1
    fish = zeros((dim,dim))

    for a in range(dim-1):
        for b in range(dim-1):
            fish[a,b] = sum( HC.dmBdE[a](sne_z)*HC.dmBdE[b](sne_z) / sne_err2 )

    for i in range(dim-1):
        fish[dim-1,i] = sum(HC.dmBdE[i](sne_z)/sne_err2)
        fish[i,dim-1] = fish[dim-1,i]

    fish[dim-1,dim-1] = sum(1./sne_err2)

    cov_H_MB = inv(fish)
    cov_H = cov_H_MB[0:dim-1,0:dim-1]
    return cov_H

def get_cov_LogH(z_interp,dlogE=0.01):
    HC = LogHzCosmo(z_interp=z_interp,cosmo_array_size=200,zmax=1.7,kind='linear',dlogE=dlogE)
    HC.ComputedmB()

    dim = len(z_interp)+1
    fish = zeros((dim,dim))

    for i in range(dim-1):
        for j in range(dim-1):
            fish[i,j] = sum(HC.dmBdlogE[i](sne_z)*HC.dmBdlogE[j](sne_z)/sne_err2)

    for i in range(dim-1):
        fish[dim-1,i] = sum(HC.dmBdlogE[i](sne_z)/sne_err2)
        fish[i,dim-1] = fish[dim-1,i]

    fish[dim-1,dim-1] = sum(1./sne_err2)

    cov_LogH_MB = inv(fish)
    cov_LogH = cov_LogH_MB[0:dim-1,0:dim-1]
    return cov_LogH


z6 = array([0.17 ,  0.476,  0.782,  1.088,  1.394,  1.7])
z8 = array([0.17,  0.38857143,  0.60714286,  0.82571429,  1.04428571, 1.26285714,  1.48142857,  1.7])
z10 = array([0.17,  0.34,  0.51,  0.68,  0.85,  1.02,  1.19,  1.36, 1.53,  1.7])


# zA = linspace(0.1,1.7,6)
# zB = linspace(0.1,1.7,8)
zC = linspace(0.1,1.7,15)
zD = linspace(0.1,1.7,100)
zE = linspace(0.1,1.7,150)

zz = []
# zz.append(zA)
# zz.append(zB)
zz.append(zC)
zz.append(zD)
zz.append(zE)

scale = []
scale.append( 1.0 / float(len(zz[0]))**0.5 )
scale.append( 1.0 / float(len(zz[1]))**0.5 )
scale.append( 1.0 / float(len(zz[2]))**0.5 )
# scale.append( 1.0 / float(len(zz[3]))**0.5 )
# scale.append( 1.0 / float(len(zz[4]))**0.5 )

fig = figure(figsize=(7,5))
ax = fig.add_subplot(111)

Allvecs = []
Allvals = []

for i in range(len(zz)):
    # print 'scale[%d] = %g'%(i,scale[i])
    cov_H       = get_cov_H(zz[i],dE=1E-8)
    vals,vecs = eig(cov_H)
    vals = (sort(vals)**0.5)
    idx = argsort(vals)
    # Allvecs.append(vecs[:,idx])
    Allvals.append(vals[idx])
    # semilogy(vals*scale[i],'--s',label=r'# zbin = '+str(len(zz[i])))
    plot(vals*scale[i],'--s',label=r'# zbin = '+str(len(zz[i])))
    # cov_LogH    = get_cov_LogH(zz[i],dlogE=1E-5)
    # plot(zz[i],diag(cov_H)**0.5/HC0.Efun(zz[i])*scale[i],'--s',label=r'# zbin = '+str(len(zz[i])))
    # plot(zz[i],diag(cov_LogH)**0.5*scale[i],'-o',label=r'# zbin = '+str(len(zz[i])))


# for i in range(len(Allvecs)):
#     for j in range(9):
#         subplot(3,3,j+1)
#         plot(Allvecs[i][:,j],'-')

xlim(0,25)
# ylim(1e-3,1)
grid()
legend()

savefig('fuck.pdf')

show()    

