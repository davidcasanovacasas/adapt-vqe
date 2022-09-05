import sys
from run_adapt import run_adapt_vqe

# define output file
file_path = 'h4_gsd.60.321g.act4.out'
print("Running ADAPT-VQE, "
      "results will be saved in " + file_path)
sys.stdout = open(file_path, "w")

unit = 'Angstrom'
r = 6.0
geometry = [('H', (0,0,0)), ('H', (0,0,r)), ('H', (0,0,2*r)), ('H', (0,0,3*r))]
charge = 0
spin = 0
basis = '3-21g'
#basis = 'sto-3g'
n_act = 4
reference = 'rhf'  # rhf, uno, uhf (not ready. we should add rohf)

# pool type: SD, GSD, sc_GSD
pool_type = 'GSD'

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
              reference = reference,
              pool_type=pool_type,
              basis2=basis2,
              n_act=n_act,
              n_act2=n_act2,
              initial_ind=initial_ind,
              adapt_conver = adapt_conver,
              adapt_thresh = adapt_thresh,
              theta_thresh = theta_thresh,
              adapt_maxiter = adapt_maxiter)
