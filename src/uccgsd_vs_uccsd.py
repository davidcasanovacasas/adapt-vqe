import scipy
import openfermion
import openfermionpsi4
import os
import numpy as np
import copy
import random 
import sys

import pyscf
from pyscf import lib
from pyscf import gto, scf, mcscf, fci, ao2mo, lo, molden, cc
from pyscf.cc import ccsd

import operator_pools
import vqe_methods
from tVQE import *

from openfermion import *


import pyscf_helper

r = 1.5
#geometry = [('H', (0,0,1*r)), ('H', (0,0,2*r)), ('H', (0,0,3*r)), ('H', (0,0,4*r))]
geometry = [('H',  (0, 0, 0)), 
            ('Li', (0, 0, r*2.39))]
#geometry = [('H', (0,0,1*r)), ('H', (0,0,2*r)), ('H', (0,0,3*r)), ('H', (0,0,4*r)), ('H', (0,0,5*r)), ('H', (0,0,6*r))]


charge = 0
spin = 0
basis = 'sto-3g'

[n_orb, n_a, n_b, h, g, mol, E_nuc, E_scf, C, S] = pyscf_helper.init(geometry,charge,spin,basis)

print(" n_orb: %4i" %n_orb)
print(" n_a  : %4i" %n_a)
print(" n_b  : %4i" %n_b)

sq_ham = pyscf_helper.SQ_Hamiltonian()
sq_ham.init(h, g, C, S)
print(" HF Energy: %12.8f" %(E_nuc + sq_ham.energy_of_determinant(range(n_a),range(n_b))))

fermi_ham  = sq_ham.export_FermionOperator()

hamiltonian = openfermion.transforms.get_sparse_operator(fermi_ham)

s2 = Make_S2(n_orb)

#build reference configuration
occupied_list = []
for i in range(n_a):
    occupied_list.append(i*2)
for i in range(n_b):
    occupied_list.append(i*2+1)

print(" Build reference state with %4i alpha and %4i beta electrons" %(n_a,n_b), occupied_list)
reference_ket = scipy.sparse.csc_matrix(openfermion.jw_configuration_state(occupied_list, 2*n_orb)).transpose()

[e,v] = scipy.sparse.linalg.eigsh(hamiltonian.real,1,which='SA',v0=reference_ket.todense())
for ei in range(len(e)):
    S2 = v[:,ei].conj().T.dot(s2.dot(v[:,ei]))
    print(" State %4i: %12.8f au  <S2>: %12.8f" %(ei,e[ei]+E_nuc,S2))
fermi_ham += FermionOperator((),E_nuc)
pyscf.molden.from_mo(mol, "full.molden", sq_ham.C)

pool = operator_pools.singlet_GSD()
pool.init(n_orb)

[e,v,params] = vqe_methods.ucc(fermi_ham, pool, reference_ket, adapt_thresh=1e-6, theta_thresh=1e-9)

print(" Final energy: %12.8f" %e)
print(" <S^2> of final state  : %12.8f" %(v.conj().T.dot(s2.dot(v))[0,0].real))

