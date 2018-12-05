#!/usr/bin/python
from __future__ import print_function
from scipy.optimize import minimize
import numpy as np
import sys
import os
from multiprocessing import Pool
if __package__ == None: sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from zetafold.parameters import get_params
from zetafold.partition import partition
from zetafold.score_structure import score_structure
from zetafold.data.training_examples import *

def calc_dG_gap( training_example ):
    sequence, structure = training_example.sequence, training_example.structure
    dG_structure = score_structure( sequence, structure, params = training_example.params )

    force_base_pairs = training_example.force_base_pairs
    p = partition( sequence, params = training_example.params, suppress_all_output = True, mfe = True, force_base_pairs = force_base_pairs )
    dG = p.dG

    dG_gap = dG_structure - dG # will be a positive number, best case zero.
    print(p.struct_MFE, dG_gap)
    return dG_gap

def calc_dG_gap_deriv( training_example ):
    sequence, structure = training_example.sequence, training_example.structure
    (dG_structure, log_derivs_structure ) = score_structure( sequence, structure, params = training_example.params, deriv_params = train_parameters )

    force_base_pairs = training_example.force_base_pairs
    p = partition( sequence, params = training_example.params, suppress_all_output = True, mfe = True, force_base_pairs = force_base_pairs, deriv_params = train_parameters )
    log_derivs = p.log_derivs

    return np.array( log_derivs ) - np.array( log_derivs_structure )

def free_energy_gap( x ):
    for n,param_tag in enumerate(train_parameters): params.set_parameter( param_tag, np.exp(x[n]))
    print('\n',np.exp(x))
    all_dG_gap = pool.map( calc_dG_gap, training_examples )
    return sum( all_dG_gap )

def free_energy_gap_deriv( x ):
    for n,param_tag in enumerate(train_parameters):
        assert( param_tag in params.parameter_tags )
        params.set_parameter( param_tag, np.exp(x[n]))
    all_dG_gap_deriv = pool.map( calc_dG_gap_deriv, training_examples )
    return sum( all_dG_gap_deriv )

# set up starter parameters
params = get_params( suppress_all_output = True )
params.set_parameter( 'K_coax', 0.0 )

# set up training examples
training_examples = [ tRNA, P4P6_outerjunction ]
train_parameters  = ['C_init','C_eff_stack_WC_WC','C_eff_stack_WC_GU','C_eff_stack_WC_UG','C_eff_stack_GU_GU','C_eff_stack_UG_GU','C_eff_stack_GU_UG']
x0      = np.array( [     0.5,                5.3,                4.1,                5.6,                4.9,                4,                    4] ) * np.log(10.0)

#x0 = np.array( [0.5, 5.3, 4.1, 5.6, 4.9, 4, 4, 0.5] )
#apply_params_func = apply_params_Cinit_CeffSix_Kcoax
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P5abc_sequence, P5abc_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure, P4P6_outerjunction_force_bps), (add_sequence, add_structure) ]
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_sequence, P4P6_structure) ]

for training_example in training_examples: training_example.params = params

pool = Pool( 4 )
result = minimize( free_energy_gap, x0, method = 'L-BFGS-B') #, jac = free_energy_gap_deriv )
final_loss = free_energy_gap( result.x )
print( 'Deriv: ', free_energy_gap_deriv( result.x ) )

print(result)
print('Final parameters:', result.x, 'Loss:',final_loss)
