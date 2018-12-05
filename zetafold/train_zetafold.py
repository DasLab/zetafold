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
    all_dG_gap_deriv = map( calc_dG_gap_deriv, training_examples )
    return sum( all_dG_gap_deriv )

def apply_params_Ceff_Cinit_KdAU_KdGU( params, x ):
    q = 10.0**x
    params.C_eff_stacked_pair = q[0]
    params.initialize_C_eff_stack()
    params.C_init = q[1]
    params.base_pair_types[2].Kd = q[2]
    params.base_pair_types[3].Kd = q[2] # A-U
    params.base_pair_types[4].Kd = q[3]
    params.base_pair_types[5].Kd = q[3] # G-U
    params.K_coax = 0.0

def apply_params_Cinit_CeffSix( params, x ):
    q = 10.0**x
    params.C_init = q[0]
    bpts_WC = params.base_pair_types[0:4]
    bpt_GU  = params.base_pair_types[4]
    bpt_UG  = params.base_pair_types[5]
    for bpt1 in bpts_WC:
        for bpt2 in bpts_WC:
            params.C_eff_stack[bpt1][bpt2] = q[1]
    for bpt in bpts_WC:
        params.C_eff_stack[bpt][bpt_GU] = q[2]
        params.C_eff_stack[bpt_UG][bpt] = q[2]
        params.C_eff_stack[bpt][bpt_UG] = q[3]
        params.C_eff_stack[bpt_GU][bpt] = q[3]
    params.C_eff_stack[bpt_GU][bpt_GU] = q[4]
    params.C_eff_stack[bpt_UG][bpt_UG] = q[4]
    params.C_eff_stack[bpt_GU][bpt_UG] = q[5]
    params.C_eff_stack[bpt_UG][bpt_GU] = q[6]
    params.K_coax = 0.0

def apply_params_Cinit_CeffSix_Kcoax( params, x ):
    apply_params_Cinit_CeffSix( params, x[:-1] )
    params.K_coax = 10.0**x[-1]

def apply_params_Cinit_l_lBP_CeffSix( params, x ):
    apply_params_Cinit_CeffSix( params, np.array([x[i] for i in [0,3,4,5,6,7,8] ]) )
    q = 10.0**x
    params.l = q[1]
    params.l_BP = q[2]

def apply_params_Ceff_Cinit_KdAU_KdGU_Kcoax( params, x ):
    q = 10.0**x
    params.C_eff_stacked_pair = q[0]
    params.initialize_C_eff_stack()
    params.C_init = q[1]
    params.base_pair_types[2].Kd = q[2]
    params.base_pair_types[3].Kd = q[2] # A-U
    params.base_pair_types[4].Kd = q[3]
    params.base_pair_types[5].Kd = q[3] # G-U
    params.K_coax = q[4]

# set up starter parameters
params = get_params( suppress_all_output = True )
params.set_parameter( 'K_coax', 0.0 )

# set up training examples
training_examples = [ tRNA ]
train_parameters = ['Kd_CG']
x0 = np.array( [5] )

#x0 = np.array( [5, 1, 3, 3] )
#apply_params_func = apply_params_Ceff_Cinit_KdAU_KdGU
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [5, 4, 4, 3, 3, 3, 3] )
#apply_params_func = apply_params_Cinit_CeffSix
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [1.56, 5.4, 5, 4, 4, 4, 4] )
#apply_params_func = apply_params_Cinit_CeffSix
##sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_sequence, P4P6_structure) ]
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P5abc_sequence, P5abc_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure, P4P6_outerjunction_force_bps) ]

#x0 = np.array( [1.56, 0, 0, 5, 5, 4, 4, 4, 4] )
#apply_params_func = apply_params_Cinit_l_lBP_CeffSix
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [5, 1, 3, 3, 1] )
#apply_params_func = apply_params_Ceff_Cinit_KdAU_KdGU_Kcoax
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]
#sequence_structure_pairs  = [ (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [0.5, 5.3, 4.1, 5.6, 4.9, 4, 4, 0.5] )
#apply_params_func = apply_params_Cinit_CeffSix_Kcoax
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P5abc_sequence, P5abc_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure, P4P6_outerjunction_force_bps), (add_sequence, add_structure) ]
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_sequence, P4P6_structure) ]

for training_example in training_examples: training_example.params = params

pool = Pool( 4 )
result = minimize( free_energy_gap, x0, method = 'L-BFGS-B', jac = free_energy_gap_deriv )
final_loss = free_energy_gap( result.x )
print( 'Deriv: ', free_energy_gap_deriv( result.x ) )

print(result)
print('Final parameters:', result.x, 'Loss:',final_loss)
