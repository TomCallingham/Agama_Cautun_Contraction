import numpy as np
import agama
agama.setUnits( mass=1., length=1., velocity=1.)  # Msun, kpc, km/s

def contract_agama_potential(dm_pot, baryon_pot, fbar=0.157,rmax=500, rmin=0.1):
    '''Given spherical DM and axisymmetric baryon agama potentials,
    creates a contracted DM agama potential using procedure found in Cautun 20'''
    r_space =  np.geomspace(rmin,rmax,501)
    Mcum_dm = np.array([dm_pot.totalMass(r) for r in r_space])
    xyz = np.stack((r_space,np.zeros_like(r_space),np.zeros_like(r_space))).transpose()
    dens_dm = dm_pot.density(xyz)
    Mcum_bar = np.array([baryon_pot.totalMass(r) for r in r_space])
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
