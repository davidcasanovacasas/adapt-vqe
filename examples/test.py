import sys
from run_adapt import run_adapt_vqe

# define output file
file_path = 'test.out'
print("Running ADAPT-VQE, "
      "results will be saved in " + file_path)
#sys.stdout = open(file_path, "w")

# molecule
unit = 'Angstrom'
r = 3.0
geometry = [('H', (0,0,0)), ('H', (0,0,r))]
charge = 0
spin = 0
basis = '3-21g'
#basis = 'sto-3g'
n_act = 4

# pool type: SD, GSD, sc_GSD
pool_type = 'SD'

# Secondary pool
basis2 = None
n_act2 = None

# adapt-vqe parameters
initial_ind = [1]
adapt_conver = 'norm'  # default = 'norm'
adapt_thresh = 1e-5    # default = 1e-3
theta_thresh = 1e-7    # default = 1e-7
adapt_maxiter = 200    # default = 200

# launch calculation
run_adapt_vqe(unit,geometry,charge,spin,basis,
              pool_type=pool_type,
              basis2=basis2,
              n_act=n_act,
              n_act2=n_act2,
              initial_ind=initial_ind,
              adapt_conver = adapt_conver,
              adapt_thresh = adapt_thresh,
              theta_thresh = theta_thresh,
              adapt_maxiter = adapt_maxiter)
