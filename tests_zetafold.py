#!/usr/bin/env python
from __future__ import print_function

import argparse

from zetafold.partition import *
from zetafold.util.output_util import *
from zetafold.parameters import get_params_from_file
from zetafold.score_structure import score_structure

def test_bad_params_none(verbose = False, use_simple_recursions = False):
    print( 'Check graceful response when requesting non-existent parameter file...' )
    assert( get_params_from_file( 'blah' ) == None ) # should not exist

def test_minimal_params_exist(verbose = False, use_simple_recursions = False):
    test_params = get_params_from_file( 'minimal' )
    assert( test_params )

def test_zetafold( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    assert( test_params )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    C_eff_stacked_pair = test_params.C_eff_stack[ test_params.base_pair_types[0] ][ test_params.base_pair_types[0] ]

def test_final_partition_1( verbose = False, use_simple_recursions = False ):
    # test of sequences where we know the final partition function.
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = 'CNNNGNN' # CIRCLE!
    p = partition( sequence, circle = True, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref   = C_init  * (l**7) * (1 + (C_init * l_BP**2) / Kd ) / C_std
    bpp_ref = (C_init * l_BP**2/ Kd) / ( 1 + C_init * l_BP**2/ Kd)
    deriv_parameters = ('Kd','Kd_matchlowercase','Kd_GC' ,'Kd_CG','l','l_BP','C_init','C_eff_stacked_pair')
    log_derivs_ref   = (-bpp_ref,0,-bpp_ref,-bpp_ref, 7, 2*bpp_ref, 1 + bpp_ref, 0 )
    output_test( p, Z_ref, [0,4], bpp_ref, deriv_parameters, log_derivs_ref )

def test_final_partition_2( verbose = False, use_simple_recursions = False ):
    # test of sequences where we know the final partition function.
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = 'CNNNGNN' # CIRCLE!
    structure= '(...)..'
    p = partition( sequence, circle = True, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions, structure = structure )
    Z_ref = C_init  * (l**7) * (C_init * l_BP**2) / Kd / C_std
    bpp_ref = 1.0
    deriv_parameters = ('Kd','Kd_matchlowercase','Kd_GC' ,'Kd_CG','l','l_BP','C_init','C_eff_stacked_pair')
    log_derivs_ref   = (-1, 0, -1, -1, 7, 2, 2, 0 )
    output_test( p, Z_ref, [0,4], bpp_ref, deriv_parameters, log_derivs_ref )

def test_final_partition_3( verbose = False, use_simple_recursions = False ):
    # test of sequences where we know the final partition function.
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = 'CNG'
    p = partition( sequence, params = test_params, mfe = True, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    assert( p.bps_MFE == [(0,2)] )
    Z_ref = 1 + C_init * l**2 * l_BP/ Kd
    bpp_ref = (C_init * l**2 * l_BP/Kd)/( 1 + C_init * l**2 * l_BP/Kd )
    output_test( p, Z_ref, [0,2], bpp_ref )

def test_final_partition_bimolecular_1( verbose = False, use_simple_recursions = False ):
    # test of sequences where we know the final partition function.
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequences = ['C','G']
    p = partition( sequences, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions ) # note that Z sums over only base pair (not dissociated strands!)
    output_test( p, C_std/ Kd, \
                 [0,1], 1.0 )

def test_final_partition_bimolecular_2( verbose = False, use_simple_recursions = False ):
    # test of sequences where we know the final partition function.
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    C_eff_stacked_pair = test_params.C_eff_stack[ test_params.base_pair_types[0] ][ test_params.base_pair_types[0] ]
    sequences = ['GC','GC']
    p = partition( sequences, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (C_std/Kd)*(2 + l**2 * l_BP**2 *C_init/Kd + C_eff_stacked_pair/Kd )
    bpp_ref = (1 + l**2 * l_BP**2 * C_init/Kd + C_eff_stacked_pair/Kd )/(2 + l**2 * l_BP**2 *C_init/Kd + C_eff_stacked_pair/Kd )
    log_deriv_C_init = (l**2 * l_BP**2 * C_init/Kd ) / (2 + (l**2 * l_BP**2 *C_init/Kd) + C_eff_stacked_pair/Kd )
    log_deriv_l = 2 *  log_deriv_C_init
    log_deriv_C_eff_stacked_pair = (C_eff_stacked_pair/Kd) / (2 + (l**2 * l_BP**2 *C_init/Kd) + C_eff_stacked_pair/Kd )
    deriv_parameters = ('C_init','l','l_BP','C_eff_stacked_pair','C_eff_stack_GC_GC','C_eff_stack_CG_CG','C_eff_stack_CG_GC','C_eff_stack_GC_CG')
    log_derivs_ref =  [ log_deriv_C_init, log_deriv_l, log_deriv_l, log_deriv_C_eff_stacked_pair, 0,0,0, log_deriv_C_eff_stacked_pair ]
    output_test( p, Z_ref, [0,3], bpp_ref, deriv_parameters, log_derivs_ref )

def test_without_all_extra_base_pairs( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    print( 'Testing structure calc, without and with all_extra_base_pairs' )
    sequences = ['GC','GC']
    structure = '(..)'
    p = partition( sequences, structure = structure, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (C_std/Kd)
    output_test( p, Z_ref, [0,3], 1.0 )
    output_test( p, Z_ref, [1,2], 0.0 )

def test_with_all_extra_base_pairs( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    C_eff_stacked_pair = test_params.C_eff_stack[ test_params.base_pair_types[0] ][ test_params.base_pair_types[0] ]
    sequences = ['GC','GC']
    structure = '(..)'
    p = partition( sequences, structure = structure, allow_extra_base_pairs = True, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (C_std/Kd)*(1 + l**2 * l_BP**2 *C_init/Kd + C_eff_stacked_pair/Kd )
    output_test( p, Z_ref, [0,3], 1.0 )
    output_test( p, Z_ref, [1,2], (l**2 * l_BP**2 *C_init/Kd + C_eff_stacked_pair/Kd)/( 1 + l**2 * l_BP**2 *C_init/Kd + C_eff_stacked_pair/Kd ) )

def test_nonuniform_C_eff_stacked_pair( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    # what if C_eff_stacked_pair is not uniform
    sequences = ['Ga','aC']
    test_params_C_eff_stack = get_params_from_file( 'minimal' )
    cross_C_eff_stacked_pair = 1.0  # default is 1.0e4. Now having the aa base pair next to the G-C base pair is worth 10,000-fold less.
    for base_pair_type_GC in test_params_C_eff_stack.base_pair_types[1:3]:
        test_params_C_eff_stack.C_eff_stack[ base_pair_type_GC ][  test_params_C_eff_stack.base_pair_types[0] ]= cross_C_eff_stacked_pair
        test_params_C_eff_stack.C_eff_stack[  test_params_C_eff_stack.base_pair_types[0] ][ base_pair_type_GC ] = cross_C_eff_stacked_pair
    p = partition( sequences, params = test_params_C_eff_stack, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (C_std/Kd)*(2 + l**2 * l_BP**2 *C_init/Kd + cross_C_eff_stacked_pair/Kd )
    bpp_ref = (1 + l**2 * l_BP**2 * C_init/Kd + cross_C_eff_stacked_pair/Kd )/(2 + l**2 * l_BP**2 *C_init/Kd + cross_C_eff_stacked_pair/Kd )
    log_deriv_l = 2 * (l**2 * l_BP**2 * C_init/Kd ) / (2 + (l**2 * l_BP**2 *C_init/Kd) + cross_C_eff_stacked_pair/Kd )
    log_deriv_C_eff_stacked_pair = (cross_C_eff_stacked_pair/Kd) / (2 + (l**2 * l_BP**2 *C_init/Kd) + cross_C_eff_stacked_pair/Kd )
    log_deriv_C_init = l**2 * l_BP**2 *C_init/Kd /  (2 + (l**2 * l_BP**2 *C_init/Kd) + cross_C_eff_stacked_pair/Kd )
    deriv_parameters = ('C_init','l','l_BP','C_eff_stack_GC_GC','C_eff_stack_CG_CG','C_eff_stack_CG_GC','C_eff_stack_GC_CG','C_eff_stack_GC_matchlowercase','C_eff_stack_matchlowercase_GC','Kd_CG','Kd_matchlowercase')
    log_derivs_ref =  [ log_deriv_C_init,log_deriv_l, log_deriv_l, 0,0,0,0, log_deriv_C_eff_stacked_pair, 0, -bpp_ref, -bpp_ref ]
    output_test( p, Z_ref, [0,3], bpp_ref, deriv_parameters, log_derivs_ref )

def test_CNGGC( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = 'CNGGC'
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose,  use_simple_recursions = use_simple_recursions )
    Z_ref = 1 + C_init * l**2 *l_BP/Kd * ( 2 + l )
    bpp_ref = C_init*l**2*l_BP/Kd /(  1+C_init*l**2*l_BP/Kd * ( 2 + l ))
    output_test( p, Z_ref, [0,2], bpp_ref )

def test_CNGGC_structure( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = 'CNGGC'
    structure= '(..).'
    p = partition( sequence, params = test_params, structure = structure, calc_bpp = True, suppress_bpp_output = True,
                   verbose = verbose,  use_simple_recursions = use_simple_recursions )
    output_test( p,  C_init * l**2 *l_BP/Kd * l, \
                 [0,2], 0.0 )

def test_CGNCG( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    C_eff_stacked_pair = test_params.C_eff_stack[ test_params.base_pair_types[0] ][ test_params.base_pair_types[0] ]
    sequence = 'CGNCG'
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions, mfe = True )
    Z_ref = 1 + C_init*l**2*l_BP/Kd + C_init*l**4*l_BP/Kd  + C_init**2 * (l_BP**3) * l**4 /Kd /Kd + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd /Kd
    bpp_ref = ( C_init*l**4*l_BP/Kd  + C_init**2 * (l_BP**3) * l**4 /Kd /Kd  + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd /Kd) / ( 1 + C_init*l**2*l_BP/Kd + C_init*l**4*l_BP/Kd  + C_init**2 * (l_BP**3) * l**4 /Kd /Kd + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd /Kd )
    output_test( p, Z_ref, [0,4], bpp_ref )

def test_MFE_tie( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    # an example with ties for MFE structure
    print( 'Example with ties for MFE structure...' )
    sequence = 'CNGNC'
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions, mfe = True )
    Z_ref = 1 + 2 * C_init*l**2*l_BP/Kd
    bpp_ref = C_init*l**2*l_BP/Kd/ Z_ref
    output_test( p, Z_ref, [0,2], bpp_ref )

def test_enumeration( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    print( 'Enumeration tests...' )
    sequence = 'CNGCNG'
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, do_enumeration = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (1 + C_init * l**2 *l_BP/Kd)**2  + C_init * l**5 * l_BP/Kd + (C_init * l**2 *l_BP/Kd)**2 * K_coax
    bpp_ref = (C_init * l**2 *l_BP/Kd*(1 + C_init * l**2 *l_BP/Kd) + (C_init * l**2 *l_BP/Kd)**2 * K_coax) / Z_ref
    deriv_parameters = ('C_eff_stacked_pair','Kd')
    log_deriv_Kd = -( 2*C_init * l**2 *l_BP/Kd + 2*(C_init * l**2 *l_BP/Kd)**2  + C_init * l**5 * l_BP/Kd + 2*(C_init * l**2 *l_BP/Kd)**2 * K_coax ) / Z_ref
    log_deriv_ref    = (0, log_deriv_Kd)
    output_test( p, Z_ref,[0,2], bpp_ref, deriv_parameters, log_deriv_ref )
    assert( set(p.struct_enumerate) == set(['......', '(.)...', '(....)', '...(.)', '(.)(.)']) )

def test_structure_constrained_score_deriv( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = 'CNGCNG'
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, do_enumeration = True, verbose = verbose, use_simple_recursions = use_simple_recursions )

    print( 'Stringent test of structure-constrained scores & derivs' )
    Z_tot_ref = p.Z
    Z_enumerate = []
    structures = ['......', '(.)...', '(....)', '...(.)', '(.)(.)']
    Z_refs     = [1, C_init * l**2 *l_BP/Kd, C_init * l**5 * l_BP/Kd,  C_init * l**2 *l_BP/Kd, (C_init * l**2 *l_BP/Kd)**2 * (1+K_coax)]
    bpp_refs_0_2=[0,1,0,0,1]
    deriv_params =   ['Kd_CG','C_init','l','l_BP','C_eff_stacked_pair','K_coax','l_coax']
    log_derivs_ref = [ [ 0,0,0,0,0,0,0],
                       [-1,1,2,1,0,0,0],
                       [-1,1,5,1,0,0,0],
                       [-1,1,2,1,0,0,0],
                       [-2,2,4,2,0,K_coax/(1+K_coax),0] ]
    for n,structure in enumerate( structures ):
        p = partition( sequence, structure = structure, params = test_params, calc_bpp = True, suppress_bpp_output = True, do_enumeration = False, verbose = verbose, use_simple_recursions = use_simple_recursions, deriv_params = deriv_params )
        output_test( p, Z_refs[n], [0,2], bpp_refs_0_2[n], deriv_params, log_derivs_ref[n] )
        # also throw in a test of score_structure here
        ( dG, log_derivs ) = score_structure( sequence, structure, params = test_params, deriv_params = deriv_params )
        assert_equal( dG, p.dG )
        for log_deriv,log_deriv_ref in zip( log_derivs, log_derivs_ref[n] ): assert_equal( log_deriv, log_deriv_ref )
        Z_enumerate.append( p.Z )
    assert_equal( sum(Z_enumerate), Z_tot_ref )

def test_extended_alphabet_coaxial_strained_3WJ( verbose = False, use_simple_recursions = False ):
    print( 'Testing extended alphabet & coaxial stacks...')
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    sequence = ['xy','yz','zx']
    params_allow_strained_3WJ = get_params_from_file( 'minimal' )
    params_allow_strained_3WJ.allow_strained_3WJ = True
    p = partition( sequence, params = params_allow_strained_3WJ, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = 3*(C_std/Kd)**2 * (1 + K_coax)  + \
            (C_std/Kd)**2 * (C_init/Kd) * l**3 * l_BP**3  + \
            3*(C_std/Kd)**2 * (C_init/Kd) * K_coax * l_coax*l**2 * l_BP
    bpp_ref = ( 2 * (C_std/Kd)**2 * (1 + K_coax) + \
                (C_std/Kd)**2 * (C_init/Kd) * l**3 * l_BP**3 + \
                3*(C_std/Kd)**2 * (C_init/Kd) * K_coax * l_coax*l**2 * l_BP ) / Z_ref
    deriv_parameters = [ 'l', 'K_coax', 'l_coax' ]
    log_derivs_ref = [ ( 3*(C_std/Kd)**2 * (C_init/Kd) * l**3 * l_BP**3  + \
                         2* 3*(C_std/Kd)**2 * (C_init/Kd) * K_coax * l_coax*l**2 * l_BP) / Z_ref,\
                       ( 3*(C_std/Kd)**2 * K_coax  +  \
                         3*(C_std/Kd)**2 * (C_init/Kd) * K_coax * l_coax*l**2 * l_BP ) / Z_ref,\
                       3*(C_std/Kd)**2 * (C_init/Kd) * K_coax * l_coax*l**2 * l_BP / Z_ref ]
    output_test( p, Z_ref, [1,2], bpp_ref, deriv_parameters, log_derivs_ref  )

def test_extended_alphabet_coaxial( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    # testing extended alphabet & coaxial stacks
    sequence = ['xy','yz','zx']
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = 3*(C_std/Kd)**2 * (1 + K_coax)  + \
            (C_std/Kd)**2 * (C_init/Kd) * l**3 * l_BP**3
    bpp_ref = ( 2 * (C_std/Kd)**2 * (1 + K_coax) + \
                (C_std/Kd)**2 * (C_init/Kd) * l**3 * l_BP**3 ) / Z_ref
    output_test( p, Z_ref, [1,2], bpp_ref  )

def test_Z_final_no_regression( verbose = False, use_simple_recursions = False ):
    test_params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd = test_params.base_pair_types[0].Kd
    # test that caught a bug in Z_final
    sequence = 'NyNyxNx'
    p = partition( sequence, params = test_params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (1 + C_init * l**2 *l_BP/Kd)**2  +(C_init * l**2 *l_BP/Kd)**2 * K_coax
    bpp_ref = ( C_init * l**2 *l_BP/Kd * (1 + C_init * l**2 *l_BP/Kd)  + (C_init * l**2 *l_BP/Kd)**2 * K_coax ) / Z_ref
    output_test( p, Z_ref, [1,3], bpp_ref  )

def test_sfxn_with_motif( verbose = False, use_simple_recursions = False ):
    # Try adding a 'motif' to our scorefunction
    print( 'Testing motif (single nt bulge)' )
    params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = params.get_variables()
    Kd = params.base_pair_types[0].Kd
    params.set_parameter( 'K_coax', 0.0 )
    C_eff_motif = 10.0
    params.set_parameter( 'C_eff_motif_CG_CAG', C_eff_motif )

    sequences = ['CG','CAG']
    p = partition( sequences, params = params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (C_std/Kd)*(2 + l**3 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    bpp_ref = (1 + l**3 * l_BP**2 * C_init/Kd + C_eff_motif/Kd )/(2 + l**3 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    log_deriv_C_init = (l**3 * l_BP**2 * C_init/Kd) / (2 + l**3 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    log_deriv_l = 3 *  log_deriv_C_init
    log_deriv_C_eff_stacked_pair = 0
    log_deriv_C_eff_motif = (C_eff_motif/Kd) / (2 + l**3 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    deriv_parameters = ('C_init','l','C_eff_stacked_pair','C_eff_motif_CG_CAG' )
    log_derivs_ref =  [ log_deriv_C_init, log_deriv_l, log_deriv_C_eff_stacked_pair, log_deriv_C_eff_motif ]
    output_test( p, Z_ref, [0,4], bpp_ref, deriv_parameters, log_derivs_ref )


def test_sfxn_with_symmetric_motif( verbose = False, use_simple_recursions = False ):
    # Try adding symmetric motif
    params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = params.get_variables()
    Kd = params.base_pair_types[0].Kd
    sequences = ['CAG','CAG']
    C_eff_motif = 10.0
    params.set_parameter( 'C_eff_motif_NNN_NNN', C_eff_motif )
    p = partition( sequences, params = params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions, do_enumeration = True )
    Z_ref = (C_std/Kd)*(2 + l**4 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    bpp_ref = (1 + l**4 * l_BP**2 * C_init/Kd + C_eff_motif/Kd )/(2 + l**4 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    log_deriv_C_init = (l**4 * l_BP**2 * C_init/Kd) /(2 + l**4 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    log_deriv_l = 4 * log_deriv_C_init
    log_deriv_C_eff_stacked_pair = 0
    log_deriv_C_eff_motif = (C_eff_motif/Kd) /(2 + l**4 * l_BP**2 *C_init/Kd + C_eff_motif/Kd )
    deriv_parameters = ('C_init','l','C_eff_stacked_pair','C_eff_motif_NNN_NNN' )
    log_derivs_ref =  [ log_deriv_C_init, log_deriv_l, log_deriv_C_eff_stacked_pair, log_deriv_C_eff_motif ]
    output_test( p, Z_ref, [0,5], bpp_ref, deriv_parameters, log_derivs_ref )

def test_sfxn_with_tetraloop_hairpin( verbose = False, use_simple_recursions = False ):
    print( 'Testing motif (tetraloop hairpin)' )
    params = get_params_from_file( 'minimal' )
    (C_init, l, l_BP, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = params.get_variables()
    Kd = params.base_pair_types[0].Kd

    C_eff_motif = 10.0
    params.set_parameter( 'C_eff_motif_NNN_NNN', C_eff_motif )
    params.set_parameter( 'C_eff_motif_CAAAAG', C_eff_motif )
    sequence = 'CAAAAG'
    p = partition( sequence, params = params, calc_bpp = True, suppress_bpp_output = True, verbose = verbose, use_simple_recursions = use_simple_recursions, do_enumeration = True )
    Z_ref = 1 + (C_init/Kd)*l**5 * l_BP  +  C_eff_motif/Kd
    bpp_ref = ((C_init/Kd)*l**5 * l_BP  +  C_eff_motif/Kd) / Z_ref
    log_deriv_C_init = (C_init/Kd)*l**5 * l_BP / Z_ref
    log_deriv_l = 5 * log_deriv_C_init
    log_deriv_C_eff_stacked_pair = 0
    log_deriv_C_eff_motif = (C_eff_motif/Kd) / Z_ref
    deriv_parameters = ('C_init','l','C_eff_stacked_pair','C_eff_motif_NNN_NNN','C_eff_motif_CAAAAG' )
    log_derivs_ref =  [ log_deriv_C_init, log_deriv_l, log_deriv_C_eff_stacked_pair, 0, log_deriv_C_eff_motif ]
    output_test( p, Z_ref, [0,5], bpp_ref, deriv_parameters, log_derivs_ref )

def test_deprecated_zetafold_v020( verbose = False, use_simple_recursions = False ):
    params = get_params( 'deprecated/zetafold_v0.20' )
    C_eff_motif = 10.0
    params.set_parameter( 'C_eff_motif_NN_NNN', C_eff_motif )
    p = partition( 'GCGACGUUG',deriv_check=True, params=params, structure = '.(.(...))',verbose=verbose, use_simple_recursions = use_simple_recursions )
    # rhiju -- you never do anything with this; is that intentional?

def test_secstruct(verbose = False, use_simple_recursions = False):
    # test secstruct
    assert( secstruct_from_bps( [(0,5),(1,4)],7 ) == '((..)).' )
    assert( bps_from_secstruct(  '((..)).' ) == [(0,5),(1,4)] )
    assert( parse_motifs( '.(((.)(.))).' ) == [[[1, 2], [9, 10]], [[2, 3], [5, 6], [8, 9]], [[3, 4, 5]], [[6, 7, 8]], [[10, 11, 0, 1]]] )
    assert( parse_motifs( '(((.)(.))).'  ) == [[[0, 1], [8, 9]], [[1, 2], [4, 5], [7, 8]], [[2, 3, 4]], [[5, 6, 7]], [[9, 10, 0]]] )
    assert( parse_motifs( '.(((.)(.)))'  ) == [[[1, 2], [9, 10]], [[2, 3], [5, 6], [8, 9]], [[3, 4, 5]], [[6, 7, 8]], [[10, 0, 1]]] )
    assert( parse_motifs( '(((.)(.)))'   ) == [[[0, 1], [8, 9]], [[1, 2], [4, 5], [7, 8]], [[2, 3, 4]], [[5, 6, 7]], [[9, 0]]] )

def test_score_structure( verbose = False, use_simple_recursions = False ):
    # score_structure
    sequence = 'GCUCAGUUGGGAGAGC'
    structure= '((((........))))'
    print("Testing score_structure on short sequence, full parameters: ", sequence, structure)
    dG = score_structure( sequence, structure, test_mode = True )

def test_score_structure_tRNA( verbose = False, use_simple_recursions = False ):
    print()
    print("Testing score_structure on tRNA folding, full parameters")
    sequence = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
    structure= '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
    dG = score_structure( sequence, structure )
    dG_new = score_structure( sequence+sequence, structure+structure )
    print("Check also double-sequence and double structure get 2*Z", dG_new, 2*dG)
    assert_equal( dG_new, 2*dG )

def test_deriv_check_deprecated_020( verbose = False, use_simple_recursions = False ):
    params = get_params( 'deprecated/zetafold_v0.20' )
    C_eff_motif = 10.0
    params.set_parameter( 'C_eff_motif_NN_NNN', C_eff_motif )
    print()
    print( 'Do deriv-check on small but complex case' )
    sequence = 'GCUCAGUGAGAGC'
    # This seems nonsensical, but back when this was all one function, this was the most recent
    # value of 'structure'
    structure= '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
    print("Testing partition on short sequence, full parameters with AA added in for good measure: ", sequence, structure)
    #params = get_params()
    params.set_parameter( 'Kd_AA', 1000 )
    params.set_parameter( 'C_eff_stack_CG_AA', 100000 )
    dG = partition( sequence, deriv_check=True, params = params  ) # deriv_check runs asserts

def test_deriv_check_deprecated_018( verbose = False, use_simple_recursions = False ):
    print()
    print( 'Do deriv-check on small but complex case' )
    sequence = 'GCUCAGUGAGAGC'
    # This seems nonsensical, but back when this was all one function, this was the most recent
    # value of 'structure'
    structure= '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
    print("Testing partition on short sequence, full parameters with AA added in for good measure, based of zetafold v0.18: ", sequence, structure)
    params = get_params( 'v0.18' )
    params.set_parameter( 'Kd_AA', 1000 )
    params.set_parameter( 'C_eff_stack_CG_AA', 100000 )
    dG = partition( sequence, deriv_check=True, params = params  ) # deriv_check runs asserts

def test_deriv_check_new_031( verbose = False, use_simple_recursions = False ):
    params = get_params( 'zetafold_v0.31' )
    sequence = 'GCUCAGUGAGAGC'
    structure= '(((((...)))))'
    print("Testing partition on short sequence, full zetafold_v031: ", sequence, structure)
    dG = partition( sequence, deriv_check=True, params = params  ) # deriv_check runs asserts

def all_tests_zetafold(verbose, use_simple_recursions):
    for key, value in globals().items():
        if callable(value) and key.startswith('test_'):
            value(verbose, use_simple_recursions)

if __name__=='__main__':
    parser = argparse.ArgumentParser( description = "Test nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument("-v","--verbose", action='store_true', default=False, help='output dynamic programming matrices')
    parser.add_argument("--simple", action='store_true', default=False, help='Use simple recursions (fast!)')
    args     = parser.parse_args()
    all_tests_zetafold( verbose = args.verbose, use_simple_recursions = args.simple )

