import agama
agama.setUnits(mass=1., length=1., velocity=1.)  # Msun, kpc, km/s
# Units 1 Msun, 1 Kpc, 1 km/s
# Cautun DM halo
fb = 4.825 / 30.7  # Planck 1 baryon fraction
m200 = 0.969e12  # the DM halo mass
conc = 8.76
NFW_rho0 = 3486926.735447284
NFW_rs = 25.20684733101539
# Note subtletly in Cautun20 NFW definition, scaled overdensity changes from paper value!

# Cautun Bulge
r0_bulge = 0.075
rcut_bulge = 2.1
rho0_bulge = 103
# Cautun Stellar Discs
zd_thin = 0.3  # kcp
Rd_thin = 2.63  # kpc
Sigma0_thin = 731. * (1e6)  # Msun / kpc
zd_thick = 0.9  # kpc
Rd_thick = 3.80  # kpc
Sigma0_thick = 101. * (1e6)
# Cautun Gas Discs
Rd_HI = 7.  # kpc
Rm_HI = 4.
zd_HI = 0.085  # kpc
Sigma0_HI = 53 * (1e6)
Rd_H2 = 1.5  # kpc
Rm_H2 = 12.
zd_H2 = 0.045  # kpc
Sigma0_H2 = 2200 * (1e6)
# Cautun CGM
A = 0.19
Beta = 1.46
critz0 = 127.5e-9
R200 = 219  # kpc
cgm_amp = 200 * critz0 * A * fb

uncontracted_DM_pot = agama.Potential(type="Spheroid",
                                      densityNorm=NFW_rho0,
                                      axisRatioZ=1,
                                      gamma=1,
                                      beta=3,
                                      scaleRadius=NFW_rs)

thin_disc_pot = agama.Potential(type="disk",
                                surfaceDensity=Sigma0_thin,
                                scaleRadius=Rd_thin,
                                scaleHeight=zd_thin)

thick_disc_pot = agama.Potential(type="disk",
                                 surfaceDensity=Sigma0_thick,
                                 scaleRadius=Rd_thick,
                                 scaleHeight=zd_thick)

HI_disc_pot = agama.Potential(type="disk",
                              surfaceDensity=Sigma0_HI,
                              scaleRadius=Rd_HI,
                              scaleHeight=zd_HI,
                              innerCutoffRadius=4)

H2_disc_pot = agama.Potential(type="disk",
                              surfaceDensity=Sigma0_H2,
                              scaleRadius=Rd_H2,
                              scaleHeight=zd_H2,
                              innerCutoffRadius=12)

bulge_pot = agama.Potential(type="Spheroid",
                            densityNorm=10.3e10,
                            axisRatioZ=0.5,
                            gamma=0,
                            beta=1.8,
                            scaleRadius=r0_bulge,
                            outerCutoffRadius=rcut_bulge)


cgm_pot = agama.Potential(type="Spheroid",
    densityNorm = 7.615e2,
    gamma = Beta,
    beta = Beta,
    scaleRadius = R200,
    cutoffStrength = 2,
    outerCutoffRadius = 2*R200)


disc_pot = agama.Potential(thin_disc_pot, thick_disc_pot, HI_disc_pot, H2_disc_pot)
C20_Baryon_pot = agama.Potential(disc_pot, bulge_pot, cgm_pot)

data_file = 'Contracted_Cautun20' + '.coef_mul'
try:
    print('Load Pot')
    C20_contracted_pot = agama.Potential(data_file)
except Exception:
    print("Can't load potential, creating new instance...")
    from contraction_agama import contract_agama_potential
    C20_contracted_pot = contract_agama_potential(uncontracted_DM_pot, C20_Baryon_pot)
    print('Saving new potential...')
    C20_contracted_pot.export(data_file)
    print('Saved!')

# C20_pot = agama.Potential(C20_contracted_pot, disc_pot, bulge_pot, cgm_pot) #More convenient split
C20_pot = agama.Potential(thin_disc_pot, thick_disc_pot, HI_disc_pot, H2_disc_pot, 
                          bulge_pot, cgm_pot, C20_contracted_pot)
