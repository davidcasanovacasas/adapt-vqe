import scipy
import vqe_methods
import operator_pools
import pyscf_helper
from run_adapt import run_adapt_vqe

import pyscf
from pyscf import lib
#from pyscf import gto, scf, mcscf, fci, ao2mo, lo, molden, cc
from pyscf import gto, scf, mcscf, fci, ao2mo, lo, tools, cc
from pyscf.cc import ccsd

import openfermion
from openfermion import *
from tVQE import *

import sys

# define output file
file_path = 'test.out'
print("Running ADAPT-VQE, "
      "results will be saved in " + file_path)
#sys.stdout = open(file_path, "w")


unit = 'Angstrom'
r = 3.0
geometry = [('H', (0,0,0)), ('H', (0,0,r))]
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
pool_type = 'SD'

run_adapt_vqe(unit,geometry,charge,spin,basis,
              pool_type=pool_type,
              basis2=basis2,
              n_act=n_act,
              n_act2=n_act2,
              initial_ind=initial_ind)
