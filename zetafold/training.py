#!/usr/bin/python
from __future__ import print_function
import numpy as np
import sys
import os
from .parameters import get_params
from .partition import partition
from .score_structure import score_structure
from .util.constants import KT_IN_KCAL
from scipy.optimize import check_grad

def calc_dG_gap( training_example ):
    ( sequence, structure, force_base_pairs, params, train_parameters, allow_extra_base_pairs ) = ( training_example.sequence, training_example.structure, training_example.force_base_pairs, training_example.params, training_example.train_parameters, training_example.allow_extra_base_pairs )
    dG_structure = score_structure( sequence, structure, params = params, allow_extra_base_pairs = allow_extra_base_pairs  )
    p = partition( sequence, params = params, suppress_all_output = True, mfe = True, structure = force_base_pairs, allow_extra_base_pairs = allow_extra_base_pairs )
    dG = p.dG
    dG_gap = dG_structure - dG # will be a positive number, best case zero.
    print(p.struct_MFE, training_example.name, dG_gap)
    return dG_gap

def calc_dG_gap_deriv( training_example ):
    ( sequence, structure, force_base_pairs, params, train_parameters, allow_extra_base_pairs ) = ( training_example.sequence, training_example.structure, training_example.force_base_pairs, training_example.params, training_example.train_parameters, training_example.allow_extra_base_pairs )
    (dG_structure, log_derivs_structure ) = score_structure( sequence, structure, params = params, deriv_params = train_parameters, allow_extra_base_pairs = allow_extra_base_pairs )
    p = partition( sequence, params = params, suppress_all_output = True, mfe = True, structure = force_base_pairs, allow_extra_base_pairs = allow_extra_base_pairs, deriv_params = train_parameters )
    log_derivs = p.log_derivs
    dG_gap = dG_structure - p.dG
    print(p.struct_MFE, training_example.name, dG_gap, ' in deriv' )
    return KT_IN_KCAL * ( np.array( log_derivs ) - np.array( log_derivs_structure ) )

def pack_variables( x, params, train_parameters, training_examples = None, allow_extra_base_pairs = False):
    for n,param_tag in enumerate(train_parameters):
        assert( param_tag in params.parameter_tags )
        params.set_parameter( param_tag, np.exp(x[n]))
    if not training_examples: return
    for training_example in training_examples:
        training_example.params = params
        training_example.train_parameters = train_parameters
        training_example.allow_extra_base_pairs = allow_extra_base_pairs

def free_energy_gap( x, params, train_parameters, training_examples, allow_extra_base_pairs, pool, outfile ):
    pack_variables( x, params, train_parameters, training_examples, allow_extra_base_pairs )
    params.output_to_file( 'current.params' )
    print('\n',np.exp(x))
    all_dG_gap = pool.map( calc_dG_gap, training_examples )
    sum_dG_gap = sum( all_dG_gap )
    output_info( outfile, x, sum_dG_gap )
    return sum_dG_gap

def free_energy_gap_deriv( x, params, train_parameters, training_examples, allow_extra_base_pairs, pool ):
    pack_variables( x, params, train_parameters, training_examples, allow_extra_base_pairs )
    all_dG_gap_deriv = pool.map( calc_dG_gap_deriv, training_examples )
    return sum( all_dG_gap_deriv )

def output_info( outfile, x, sum_dG_gap ):
    if outfile == None: return
    fid = open( outfile, 'a' )
    fid.write( '%12.6f' % sum_dG_gap )
    for val in x: fid.write( '%25.6f' % val )
    fid.write( '\n' )

def create_outfile( outfile, params, training_params ):
    if outfile == None: return
    fid = open( outfile, 'w' )
    fid.write( '%12s' % 'Loss' )
    for param_tag in training_params: fid.write( '%25s' % param_tag )
    fid.write( '\n' )


def train_deriv_check( x0, loss, grad, train_parameters ):
    # Not enough output from scipy.check_grad, so I wrote my own deriv_check
    #print( 'Deriv error: ', check_grad( loss, grad, x0 ) )
    loss_val = loss( x0 )
    analytic_grad_val = grad( x0 )
    numerical_grad_val = []
    epsilon = 1.0e-8
    for n in range( len(x0) ):
        x = np.array( x0 ) # need to make an actual copy
        x[ n ] += epsilon
        numerical_grad_val.append( ( loss( x ) - loss_val ) / epsilon )
    print()
    print( '%20s %25s' % ('','-dG/d(log parameter)' ) )
    print( '%20s %25s %25s %25s' % ('parameter','analytic','numerical','diff' ) )
    for i,parameter in enumerate(train_parameters):
           print( '%20s %25.12f %25.12f %25.12f' % (parameter, analytic_grad_val[i], numerical_grad_val[i], analytic_grad_val[i] - numerical_grad_val[i]  ) )
    print()

    exit()

def get_bounds( train_parameters ):
    '''
    Based on tag, set 'reasonable' bounds.
    '''
    bounds = []
    for tag in train_parameters:
        if len( tag )>=2 and tag[:2] == 'Kd':
            bounds.append( (np.log(100.0), np.log(10000.0) ) )
        elif len( tag )>=5 and tag[:5] == 'C_eff':
            bounds.append( (np.log(100.0), np.log(10000.0) ) )
        elif tag == 'C_init':
            bounds.append( (np.log(0.1), np.log(100.0) ) )
        elif len(tag)>0 and tag[:1] == 'l':
            bounds.append( (np.log(0.1), np.log(10.0) ) )
        elif tag == 'K_coax':
            bounds.append( (np.log(0.1), np.log(100.0) ) )
        else:
            print( 'Unrecognized tag!!', tag )
            exit()
    return bounds
