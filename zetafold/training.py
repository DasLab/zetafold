#!/usr/bin/python
from __future__ import print_function
import numpy as np
import sys
import os
from .parameters import get_params
from .partition import partition
from .score_structure import score_structure


def calc_dG_gap( training_example ):
    ( sequence, structure, force_base_pairs, params, train_parameters ) = ( training_example.sequence, training_example.structure, training_example.force_base_pairs, training_example.params, training_example.train_parameters )
    dG_structure = score_structure( sequence, structure, params = params )
    p = partition( sequence, params = params, suppress_all_output = True, mfe = True, force_base_pairs = force_base_pairs )
    dG = p.dG
    dG_gap = dG_structure - dG # will be a positive number, best case zero.
    print(p.struct_MFE, training_example.name, dG_gap)
    return dG_gap

def calc_dG_gap_deriv( training_example ):
    ( sequence, structure, force_base_pairs, params, train_parameters ) = ( training_example.sequence, training_example.structure, training_example.force_base_pairs, training_example.params, training_example.train_parameters )
    (dG_structure, log_derivs_structure ) = score_structure( sequence, structure, params = params, deriv_params = train_parameters )
    p = partition( sequence, params = params, suppress_all_output = True, mfe = True, force_base_pairs = force_base_pairs, deriv_params = train_parameters )
    log_derivs = p.log_derivs

    return np.array( log_derivs ) - np.array( log_derivs_structure )

def pack_variables( x, params, train_parameters, training_examples ):
    for n,param_tag in enumerate(train_parameters):
        assert( param_tag in params.parameter_tags )
        params.set_parameter( param_tag, np.exp(x[n]))
    for training_example in training_examples:
        training_example.params = params
        training_example.train_parameters = train_parameters

def free_energy_gap( x, params, train_parameters, training_examples, pool, outfile ):
    pack_variables( x, params, train_parameters, training_examples )
    print('\n',np.exp(x))
    all_dG_gap = pool.map( calc_dG_gap, training_examples )
    sum_dG_gap = sum( all_dG_gap )
    output_info( outfile, x, sum_dG_gap )
    return sum_dG_gap

def free_energy_gap_deriv( x, params, train_parameters, training_examples, pool ):
    pack_variables( x, params, train_parameters, training_examples )
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

