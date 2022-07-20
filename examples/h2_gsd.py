import scipy
import vqe_methods
import operator_pools
import pyscf_helper

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
file_path = 'h2_gsd.basis2.out'
print("Running ADAPT-VQE"
      "results will be saved in " + file_path)
sys.stdout = open(file_path, "w")

def test():
    r = 3.5
    geometry = [('H', (0,0,0)), ('H', (0,0,r))]
    charge = 0
    spin = 0
    basis  = '3-21g'
    #basis2 = None
    basis2 = 'sto-3g'
    initial_ind = [2]

    n_orb, n_a, n_b, h, g, mol, E_nuc, E_scf, C, S = pyscf_helper.init(geometry,charge,spin,basis,reference='rhf')

    print(" n_orb: %4i" %n_orb)
    print(" n_a  : %4i" %n_a)
    print(" n_b  : %4i" %n_b)

    sq_ham = pyscf_helper.SQ_Hamiltonian()
    sq_ham.init(h, g, C, S)
    print(" HF Energy: %12.8f" %(E_nuc + sq_ham.energy_of_determinant(range(n_a),range(n_b))))

    # convert from spation orbitals to spin orbitals
    fermi_ham  = sq_ham.export_FermionOperator()
    #print("Fermionic Hamiltonian")
    #print(fermi_ham)

    #hamiltonian = openfermion.transforms.get_sparse_operator(fermi_ham)
    hamiltonian = openfermion.linalg.get_sparse_operator(fermi_ham)
    #print("Sparse Hamiltonian")
    #print(hamiltonian)

    # build S^2 matrix
    s2 = vqe_methods.Make_S2(n_orb)

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
    pyscf.tools.molden.from_mo(mol, "full.molden", sq_ham.C)
    #   pyscf.molden.from_mo(mol, "full.molden", sq_ham.C)

    #  Get operator pool (you can change it to: singlet_GSD() ...)
    pool = operator_pools.singlet_GSD()
    pool.init(n_orb, n_occ_a=n_a, n_occ_b=n_b, n_vir_a=n_orb-n_a, n_vir_b=n_orb-n_b)

    print("List of operators")
    for oi in range(pool.n_ops):
        orbitals = pool.op_index[oi]
        print(pool.op_index.index(orbitals), orbitals)

    print("List of terms")
    for oi in range(pool.n_ops):
        opstring = pool.get_string_for_term(pool.fermi_ops[oi])
        print(oi,opstring)

    if basis2:
        print("Performing ADAPT-VQE with basis2")
        # set n_orb2
        mol2 = mol
        mol2.basis = basis2
        mol2.build()
        n_orb2 = mol2.nao_nr()
        print(" # orbitals in basis2: %s = %4i" %(basis2, n_orb2))
        # Get basis2 operator pool
        pool_basis2 = operator_pools.singlet_GSD()
        pool_basis2.init(n_orb2, n_occ_a=n_a, n_occ_b=n_b, n_vir_a=n_orb2-n_a, n_vir_b=n_orb2-n_b)
        # Do the calculation
        [e,v,params] = vqe_methods.adapt_vqe_basis2(fermi_ham, pool, pool_basis2, reference_ket,
                                                theta_thresh=1e-9,
                                                trial_indices = initial_ind)
    else:
        print("Performing ADAPT-VQE")
        [e,v,params] = vqe_methods.adapt_vqe(fermi_ham, pool, reference_ket, theta_thresh=1e-9)

    print(" Final ADAPT-VQE energy: %12.8f" %e)
    print(" <S^2> of final state  : %12.8f" %(v.conj().T.dot(s2.dot(v))[0,0].real))


if __name__== "__main__":
    test()
