import sys
from run_adapt import run_adapt_vqe

# define output file
file_path = 'h4_gsd.30.321g.act4.out'
print("Running ADAPT-VQE, "
      "results will be saved in " + file_path)
sys.stdout = open(file_path, "w")


unit = 'Angstrom'
r = 3.0
geometry = [('H', (0,0,0)), ('H', (0,0,r)), ('H', (0,0,2*r)), ('H', (0,0,3*r))]
charge = 0
spin = 0
basis = '3-21g'
#basis = 'sto-3g'
n_act = 4

# Secondary pool
basis2 = None
n_act2 = None
initial_ind = [1]

# pool type: SD, GSD, sc_GSD
pool_type = 'GSD'

run_adapt_vqe(unit,geometry,charge,spin,basis,
              pool_type=pool_type,
              basis2=basis2,
              n_act=n_act,
              n_act2=n_act2,
              initial_ind=initial_ind)
