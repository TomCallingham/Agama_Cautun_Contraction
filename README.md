# Agama_Cautun_Contraction
Python scripts to load Cautun20 potential into Agama, and to contract DM potentials using the Cautun20 procedure.

Load the potential in python with either agama.Potential('Cautun20.ini') or from Cautun20 import C20_pot.

To contract an Agama potential:

from contraction_agama import contract_agama_potential

contracted_dm_pot = contract_agama_potential(uncontracted_dm_pot, baryon_pot)
