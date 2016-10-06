"""
Reads hdf5 file 
Converts halo masses into numpy array
Plots halo mass function, which is basically number density of haloes binned by mass. 

Semi-analytic functions (see Press-Schecter Model, Sheth-Tormen model) and
observational data that are usually compared with halos in cosmological simulations
via such distributions.
"""

import numpy as np  # Library for numerical methods    
import matplotlib.pylab as plt # Plotting routine
import h5py  # Reads HDF5 

def errorfill(x, y, yerr, color=None, alpha_fill=0.25, ax=None):
    # Error bars
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color, lw = 1.0, label = strLabel)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)
    
def ReverseCumulativeSum(x):
    return np.cumsum(x[::-1])[::-1] 


#m_particle = 9.94e+09   # 
nbins = 50


planck = h5py.File('planck_p100_1000trees.hdf5','r')
HaloMass = np.array(planck.get('forestHalos/nodeMass')) # Reading masses of each halo into an array
planck.close()
L = 1.02   #Box length, Full box= 256 Mpc according to Eve

xlim1 = HaloMass.min()
xlim2= HaloMass.max()

print  xlim1, '< Halo mass (Msol) <', xlim2 

#plt.clf()
plt.figure()


y,binEdges = np.histogram( HaloMass
, bins = np.logspace(np.log10(xlim1), np.log10(xlim2), nbins), density= False)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])


strLabel = r"FOF-haloes"
errorfill(bincenters, ReverseCumulativeSum(y/(L**3.0)), np.sqrt(ReverseCumulativeSum(y))/(L**3.0), color='darkred')


plt.yscale('log')
plt.xscale('log')
plt.ylabel(r"$n(> M)(h^3/Mpc^3)$")
plt.xlabel("M $(h^{-1} M_\odot)$")
plt.legend(loc = "upper right")


# Some pretty plotting commands below
plt.xlim(xlim1, xlim2)
plt.minorticks_on()
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)  
plt.rc('font',size=18)
plt.rc('axes',labelsize= 30)
plt.tick_params(axis='both', which='major', labelsize=25)

plt.savefig('HMF.pdf', bbox_inches='tight')
plt.show()

#===========================================================================================


#def neglog(array1d):
#    neglog = np.zeros_like(array1d)
#    
#    for x in range(array1d.size):        
#        if ( array1d[x] == 0): neglog[x] = 0
#        elif ( array1d[x] > 0): neglog[x] = np.log10(array1d[x])
#        else: neglog[x] = -np.log10(-array1d[x])
#
#    return neglog
#
#
#def dydx(y, x):
#    dx = np.gradient(x)
#    dydx = np.gradient(y, dx, edge_order=2)
#    return dydx
#    
#
#plt.figure(71)
#
#xlim1 = npartHalo_ahf[:,1].min()
#xlim2 = npartHalo_ahf[:,1].max()
#nbins = 15
#nbins = np.logspace(np.log10(xlim1*m_particle), np.log10(xlim2*m_particle), nbins)
#
#y,binEdges = np.histogram( npartHalo_l3[:,1]*m_particle
#, bins = nbins, density= False)
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#
#M = bincenters
#n = (y/(L**3.0))
#
##plt.plot( np.log10(M), neglog( dydx(n, np.log10(M)) ) , lw = 1.5, color='k', label = r"$\lambda_3$-haloes")
#
#plt.plot( np.log10(M)[:-1], neglog( n )[:-1] , lw = 1.5, color='k', label = r"$\lambda_3$-haloes")
#
#
#plt.minorticks_on()
##plt.yscale('symlog')
##plt.xscale('symlog')
##plt.xlim(2.5e11,2.5e14)
#plt.ylabel(r"$dn/dlogM [h^3/Mpc^3]$")
#plt.xlabel("log(M $[h^{-1} M_\odot])$")
#plt.legend(loc = "upper right")
#plt.savefig('plots/HMF_3x.pdf', bbox_inches='tight')
#
#
#
#plt.show()