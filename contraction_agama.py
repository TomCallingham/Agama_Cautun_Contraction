import numpy as np
import agama
from scipy import integrate
agama.setUnits( mass=1., length=1., velocity=1.)  # Msun, kpc, km/s

def contract_agama_potential(dm_pot, baryon_pot, fbar=0.157,rmax=500, rmin=0.1):
    '''Given spherical DM and axisymmetric baryon agama potentials,
    creates a contracted DM agama potential using procedure found in Cautun 20'''
    r_space =  np.geomspace(rmin,rmax,501)
    xyz = np.stack((r_space,np.zeros_like(r_space),np.zeros_like(r_space))).transpose()
    
    dens_dm = dm_pot.density(xyz)
    try:
        Mcum_dm = np.array([dm_pot.totalMass(r) for r in r_space])
        Mcum_bar = np.array([baryon_pot.totalMass(r) for r in r_space])
    except Exception: # Assume Old agama, no totalMass 
        Mcum_dm = Mcum_from_sph_dens(r_space, dm_pot.density)
        Mcum_bar = Mcum_from_axi_dens(r_space, baryon_pot.density)
        
    dens_bar = density_from_enclosed_mass(r_space,Mcum_bar,r_space) #note spherical average
    dens_contracted = contract_density_fit(dens_dm,dens_bar,Mcum_dm,Mcum_bar,fbar)
    def contracted_dens_func(xyz):
        r = np.linalg.norm(xyz,axis=1)
        return np.interp(r, r_space, dens_contracted)

    contracted_pot = agama.Potential(type="Multipole", density=contracted_dens_func,
                                     symmetry="spherical", rmin=1e-3, rmax=1e3)
    return contracted_pot

### Functions below lifted from Cautun Github
def contract_density_fit( density_DM, density_bar, mass_DM, mass_bar, f_bar=0.157 ):
    """ Returns the contracted DM density profile given the 'uncontracted' density and that of the baryonic distribution.
    It uses the differential (d/dr) form of Eq. (11) from Cautun et al (2020).
    
    WARNING: This function is based on a fit to the iterative calculated result of the 'contract_enclosed_mass' and works best when the
    baryonic fraction is close to the Planck (2015) value, f_bar=0.157. Do not use it if that is not the case.
   
   Args:
      density_DM    : array of DM densities. 
                          It corresponds to '(1-baryon_fraction) * density in
                          DMO (dark matter only) simulations'.
      density_bar   : array of baryonic densities.
      mass_DM       : enclosed mass in the DM component in the absence of baryons. 
                          It corresponds to '(1-baryon_fraction) * enclosed mass in
                          DMO (dark matter only) simulations'.
      mass_bar      : enclosed baryonic mass for which to calculate the DM profile.
      f_bar         : optional cosmic baryonic fraction.
   Returns:
      Array of 'contracted' DM densities.
   """
        
    eta_bar = mass_bar / mass_DM * (1.-f_bar) / f_bar  # the last two terms account for transforming the DM mass into the corresponding baryonic mass in DMO simulations
    first_factor = 0.45 + 0.41 * (eta_bar + 0.98)**0.53
    temp         = density_bar - eta_bar * density_DM * f_bar / (1.-f_bar)
    const_term   = 0.41 * 0.53 * (eta_bar + 0.98)**(0.53-1.) * (1.-f_bar) / f_bar * temp
    
    return density_DM * first_factor + const_term

def density_from_enclosed_mass( r_bins, enclosed_mass, out_r_bins ):
    """ Converts an array of enclosed masses to 3D densities.
   
   Args:
      r_bins             : the radial bins at which the enclosed mass is defined.
      enclosed_mass      : array of enclosed masses.
      out_r_bins         : array of radial distances at which the density will be interpolated.
   Returns:
      Array of densities at the location of the output radial distances.
   """
    bin_mass = enclosed_mass[1:] - enclosed_mass[:-1]
    shell_vol= 4.*np.pi/3. * (r_bins[1:]**3 - r_bins[:-1]**3)
    bin_dens = bin_mass / shell_vol
    r_vals = np.sqrt( r_bins[1:] * r_bins[:-1] )
    return np.interp( out_r_bins, r_vals, bin_dens )




def Mcum_from_sph_dens( r_space, density_xyz ):
    xyz = np.stack((r_space,np.zeros_like(r_space),np.zeros_like(r_space))).transpose()
    dens = density_xyz(xyz)
    # calculate the bin edges
    r_edges = np.zeros( len(dens)+1 )
    r_edges[0]    = 0.
    r_edges[1:-1] = 0.5 * (r_space[:-1] + r_space[1:])
    r_edges[-1]   = r_space[-1]
    
    # calculate the volume of each bin
    delta_V = 4./3. * np.pi * (r_edges[1:]**3 - r_edges[:-1]**3)
    
    return np.interp( r_space, r_edges[1:], (delta_V * dens).cumsum() )

def Mcum_from_axi_dens(rbins, density_xyz, N=500 ):  #Gives the enclosed mass in r for (R,z) func
    r1 = np.logspace( np.log10(rbins[0]*1.e-3), np.log10(rbins[-1]*1.1), N+1 )  # the bins used for the enclosed mass calculation
    r  = np.sqrt( r1[1:] * r1[:-1] )
    dr = r1[1:] - r1[:-1]
    
    def Int(theta,r_sph):
        xyz = np.stack((r_sph*np.cos(theta),np.zeros_like(r_sph), r_sph*np.sin(theta))).transpose()
        I = density_xyz(xyz) * 4*np.pi * r_sph**2 * np.cos(theta) 
        return I
    
    shellMass = np.zeros( r.shape[0] )
    for i in range( r.shape[0] ):
        shellMass[i] = integrate.quad( Int, 0., np.pi/2, args=( r[i], ) )[0] * dr[i]
    return np.interp( rbins, r1[1:], shellMass.cumsum() )
