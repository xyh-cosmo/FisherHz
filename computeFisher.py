import os,sys
sys.path.append('python')
from pylab import *
from HzCosmoBase import *

mock_sne = loadtxt('wfirst_LCDM.txt')

sne_z = mock_sne[:,0]
sne_err = mock_sne[:,2]
sne_err2 = sne_err**2

z_interp = linspace(0.15,1.7,10)
# z_interp = [0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7]
HC = HzCosmo(z_interp=z_interp,cosmo_array_size=100,zmax=1.7)
HC.ComputedmB()

dim = len(z_interp)+1
fish = zeros((dim,dim))
fishxx = zeros((dim-1,dim-1))

for i in range(dim-1):
    for j in range(dim-1):
        fish[i,j] = sum(HC.dmBdE[i](sne_z)*HC.dmBdE[j](sne_z)/sne_err2)
        fishxx[i,j] = sum(HC.dmBdE[i](sne_z)*HC.dmBdE[j](sne_z)/sne_err2)

for i in range(dim-1):
    fish[dim-1,i] = sum(HC.dmBdE[i](sne_z)/sne_err2)
    fish[i,dim-1] = fish[dim-1,i]

fish[dim-1,dim-1] = sum(1./sne_err2)

full_cov = inv(fish)
full_covxx = inv(fishxx)

Hz_cov = full_cov[0:(dim-1),0:(dim-1)]

# matshow(fish)
# matshow(Hz_cov)
# colorbar()

# semilogy(z_interp,diag(Hz_cov)**0.5,'-o',label=r'float rMB')
# semilogy(z_interp,diag(full_covxx)**0.5,'--s',label=r'fix rMB')
# legend(loc='best',frameon=False)

evals,evecs = eig(full_covxx)

IDs = argsort(evals)

subplot(1,2,1)
for i in range(5):
    plot(evecs[:,IDs[i]],'-.',label=r'sigma^2='+str(evals[IDs[i]]))
legend()

subplot(1,2,2)
for i in range(5,10):
    plot(evecs[:,IDs[i]],'-.',label=r'sigma^2='+str(evals[IDs[i]]))

legend()

show()
