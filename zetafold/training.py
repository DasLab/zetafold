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
    '''
    For one training_example, calculate the delta-G gap between target structure and ensemble over all structure.
    Should be positive. Units will be kcal/mol
    '''
    ( sequence, structure, force_base_pairs, params, train_parameters, allow_extra_base_pairs ) = ( training_example.sequence, training_example.structure, training_example.force_base_pairs, training_example.params, training_example.train_parameters, training_example.allow_extra_base_pairs )
    dG_structure = score_structure( sequence, structure, params = params, allow_extra_base_pairs = allow_extra_base_pairs  )
    p = partition( sequence, params = params, suppress_all_output = True, mfe = True, structure = force_base_pairs, allow_extra_base_pairs = allow_extra_base_pairs )
    dG = p.dG
    dG_gap = dG_structure - dG # will be a positive number, best case zero.
    print(structure, training_example.name, '[target]')
    print(p.struct_MFE, training_example.name, '[mfe]', dG_gap)
    return dG_gap

def calc_dG_gap_deriv( training_example ):
    '''
    For one training_example, calculate the derivatives w.r.t. all training parameters of
     delta-G gap between target structure and ensemble over all structure.
    Units will be kcal/mol/parameter-unit
    '''
    ( sequence, structure, force_base_pairs, params, train_parameters, allow_extra_base_pairs ) = ( training_example.sequence, training_example.structure, training_example.force_base_pairs, training_example.params, training_example.train_parameters, training_example.allow_extra_base_pairs )
    (dG_structure, log_derivs_structure ) = score_structure( sequence, structure, params = params, deriv_params = train_parameters, allow_extra_base_pairs = allow_extra_base_pairs )
    p = partition( sequence, params = params, suppress_all_output = True, mfe = True, structure = force_base_pairs, allow_extra_base_pairs = allow_extra_base_pairs, deriv_params = train_parameters )
    log_derivs = p.log_derivs
    dG_gap = dG_structure - p.dG
    print(p.struct_MFE, training_example.name, dG_gap, ' in deriv' )
    return KT_IN_KCAL * ( np.array( log_derivs ) - np.array( log_derivs_structure ) )

def pack_variables( x, params, train_parameters, training_examples = None, allow_extra_base_pairs = False):
    '''
    Set parameters based on log-params values in x.
    Also stuff information on parameters and training settings into each training_example, so that it
     can be used when running partition() derivative calculations.
    '''
    for n,param_tag in enumerate(train_parameters):
        assert( param_tag in params.parameter_tags )
        params.set_parameter( param_tag, np.exp(x[n]))
    if not training_examples: return
    for training_example in training_examples:
        training_example.params = params
        training_example.train_parameters = train_parameters
        training_example.allow_extra_base_pairs = allow_extra_base_pairs

def free_energy_gap( x, params, train_parameters, training_examples, allow_extra_base_pairs, priors, pool, outfile ):
    '''
    Main Loss function: Sum delta-G gaps over all training-examples, wrapping around calc_dG_gap.
    Handles parallelization using multiprocessing 'pool'.
    '''
    pack_variables( x, params, train_parameters, training_examples, allow_extra_base_pairs )
    params.output_to_file( 'current.params' )
    print('\n',np.exp(x))
    all_dG_gap = pool.map( calc_dG_gap, training_examples )
    sum_dG_gap = sum( all_dG_gap )
    output_info( outfile, x, sum_dG_gap )
    loss = sum_dG_gap
    if priors: loss += priors(x)[0]
    return loss

def free_energy_gap_deriv( x, params, train_parameters, training_examples, allow_extra_base_pairs, priors, pool ):
    '''
    Main gradient function: Derivative of delta-G gaps w.r.t. all training parameters, summed over all training-examples, wrapping around calc_dG_gap.
    Handles parallelization using multiprocessing 'pool'.
    '''
    pack_variables( x, params, train_parameters, training_examples, allow_extra_base_pairs )
    all_dG_gap_deriv = pool.map( calc_dG_gap_deriv, training_examples )
    deriv = sum( all_dG_gap_deriv )
    if priors: deriv += priors(x)[1]
    return deriv

BOUND_DELTA = 1.0 # tighter deltas (e.g., 0.1) lead to overflow

def eval_priors( x_list, bounds_list ):
    '''
    A prior for the log parameters that is zero within two bounds, but then rises
    quadratically outside those bounds, with log-parameter scale delta.
    '''
    val = 0
    deriv = np.zeros( len( x_list ) )
    for i,(x,bounds) in enumerate(zip( x_list, bounds_list )):
        if x < bounds[0]:
            val += ( abs(x - bounds[0]) / BOUND_DELTA )**2
            deriv[ i ] +=  2 * ( x - bounds[0] )/BOUND_DELTA**2
        if x > bounds[1]:
            val += ( abs(x - bounds[1]) / BOUND_DELTA )**2
            deriv[ i ] +=  2 * ( x - bounds[1] )/BOUND_DELTA**2
    return (val,deriv)

def output_info( outfile, x, sum_dG_gap ):
    '''
    Write loss and all log-parameter values to outfile (typically named 'training_loss.txt').
    '''
    if outfile == None: return
    fid = open( outfile, 'a' )
    fid.write( '%12.6f' % sum_dG_gap )
    for val in x: fid.write( '%25.6f' % val )
    fid.write( '\n' )

def create_outfile( outfile, params, training_params ):
    '''
    Write header_lines to outfile (typically named 'training_loss.txt').
    '''
    if outfile == None: return
    fid = open( outfile, 'w' )
    fid.write( '%12s' % 'Loss' )
    for param_tag in training_params: fid.write( '%25s' % param_tag )
    fid.write( '\n' )

def train_deriv_check( x0, loss, grad, train_parameters ):
    '''
    Not enough output from scipy.check_grad, so I wrote my own deriv_check.
    Converts all derivs to dlogZ/dlog parameter for ease of comparison with --deriv_check in zetafold.py and score_structure.py
    '''
    #print( 'Deriv error: ', check_grad( loss, grad, x0 ) ) # pretty wimpy
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
    Based on tag, set 'reasonable' bounds -- pretty arbitrarily set for now.
    '''
    bounds = []
    for tag in train_parameters:
        if len( tag )>=2 and tag[:2] == 'Kd':
            bounds.append( (np.log(100.0), np.log(100000.0) ) )
        elif len( tag )>=5 and tag[:5] == 'C_eff':
            bounds.append( (np.log(0.1), np.log(100000.0) ) )
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

def get_priors( train_parameters ):
    '''
    setup up priors function based on get_bounds()
    '''
    bounds = get_bounds( train_parameters )
    return lambda x: eval_priors( x, bounds )

def initialize_training_examples( all_training_examples, training_sets, train_data ):
    '''
    Generate list of training examples (sequences and target structures) based on database of training examples and sets

    INPUTS:
    all_training_examples = list of TrainingExamples, available from zetafold.data.training_examples
    training_sets         = list of pre-curated sets of training examples, available from zetafold.data.training_examples
    train_data            = user-inputted name of desired training_set
    '''
    if train_data == None or not train_data in training_sets.keys():
        print( '\nMust specify --train_data. Options are:' )
        for set_name in training_set_names:  print( '%30s (%s)' % (set_name,len(training_sets[set_name]) ) )
        exit()
    training_examples = [ all_training_examples[ tag ] for tag in training_sets[ train_data ] ]
    return training_examples

def initialize_train_parameters( params, train_params = None, train_params_exclude = None, no_coax = False ):
    '''
    Generate list of parameter tags that will be optimized during training

    INPUTS:
    params       = ZetaFoldParameters object
    train_params = list of training parameter tags
    train_params_exclude = list of training parameter tags to exclude
    no_coax      = do not train K_coax, l_coax
    '''
    train_parameters = train_params
    if train_parameters == None:
        if no_coax:
            if train_params_exclude == None: train_params_exclude = []
            train_params_exclude += [ 'K_coax', 'l_coax' ]
        if train_params_exclude:
            train_parameters = []
            for param_exclude in train_params_exclude:  assert( param_exclude in params.parameter_tags )
            for param in params.parameter_tags:
                if not param in train_params_exclude: train_parameters.append( param )
        elif deriv_check:
            train_parameters = params.parameter_tags
        else:
            print( '\nMust specify which parameters to optimize with --train_params or --train_params_exclude' )
            params.show_parameters()
            exit()
    return train_parameters

def initialize_parameter_values( params, train_parameters, init_params = None, init_log_params = None, use_bounds = False ):
    '''
    Generate starting values of log-parameters x0. If not supplied from command-line, initialize
     from original values in parameters file, stored in params object

    INPUTS:
    params          = ZetaFoldParameters object
    train_parameters = list of training parameter tags
    init_params     = user-supplied list of float values for each parameter
    init_log_params = user-supplied list of log values for each parameter
    '''
    x0 = None
    bounds = get_bounds( train_parameters ) if use_bounds else None
    if init_log_params:
        if len( init_log_params ) != len( train_parameters ):
            print( 'Length of --init_log_params %d should exactly match length %d of train_parameters: ' % (len(init_log_params),len(train_parameters)), train_parameters )
            exit()
        x0 = [float(log_param) for log_param in init_log_params ]
    elif init_params:
        if len( init_params ) != len( train_parameters ):
            print( 'Length of --init_params %d should exactly match length %d of train_parameters: ' % (len(init_params),len(train_parameters)), train_parameters )
            exit()
        x0 = [ np.log(float(param)) for param in  init_params]
    else:
        x0 = np.zeros( len(train_parameters) )
        for n,param_tag in enumerate(train_parameters):

            val = params.get_parameter_value( param_tag )
            if val == None:
                print( 'Did not recognize parameter: ', param_tag )
                exit()
            if val == 0.0:
                x0[ n ] = -np.Inf # to achieve 0
                print( '\nlog %s is being set to -Infinity -- you may want to initialize it to be non-zero, or set --use_priors or --use_bounds.\n' % param_tag )
            else: x0[ n ] = np.log( val )

    if bounds:
        for n,param_tag in enumerate(train_parameters):
            minval = bounds[n][0] - BOUND_DELTA
            if x0[n] < minval:
                print( '\nlog %s is below bound %s by more than bound-delta %s; raising to %s\n' % (param_tag,bounds[n][0],BOUND_DELTA,minval) )
                x0[n] = minval
            maxval = bounds[n][1] + BOUND_DELTA
            if x0[n] > maxval:
                print( '\nlog %s is above bound %s by more than bound-delta %s; lowering to %s.\n' % (param_tag,bounds[n][1],BOUND_DELTA,maxval) )
                x0[n] = maxval
    return x0
