#!/usr/bin/python
import argparse
from zetafold.partition import *
from tests_zetafold import call_all_tests 

if __name__ =='__main__':

    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("-c","-circ","--circle", action='store_true', default=False, help='Sequence is a circle')
    parser.add_argument("-params","--parameters",type=str, default='', help='Parameter file to use [default: '', which triggers latest version]')
    parser.add_argument("-struct","--structure",type=str, default=None, help='force specific structure in dot-parens notation')
    parser.add_argument("--force_base_pairs",type=str, default=None, help='force base pairs (but allow any others) in dot-parens notation')
    parser.add_argument("--mfe", action='store_true', default=False, help='Get minimal free energy structure (approximately, backtracking through partition)')
    parser.add_argument("--bpp", action='store_true', default=False, help='Get base pairing probability')
    parser.add_argument("--stochastic", type=int, default=0, help='Number of Boltzman-weighted stochastic structures to retrieve')
    parser.add_argument("--enumerate",action='store_true', default=False, help='Backtrack to get all structures and their Boltzmann weights')
    parser.add_argument("--calc_deriv", action='store_true', default=False, help='Calculate derivative with respect to all parameters')
    parser.add_argument("--no_coax", action='store_true', default=False, help='Turn off coaxial stacking')
    parser.add_argument("-v","--verbose", action='store_true', default=False, help='output dynamic programming matrices')
    parser.add_argument("--simple", action='store_true', default=False, help='Use simple recursions (slow!)')
    parser.add_argument("--calc_Kd_deriv_DP", action='store_true', default=False, help='Calculate derivative with respect to Kd_BP inline with dynamic programming [rarely used]')
    parser.add_argument( "--deriv_params",help="Parameters for which to calculate derivatives. Default: None, or all params if --calc_deriv",nargs='*')
    args     = parser.parse_args()
    if args.calc_deriv and args.deriv_params == None: args.deriv_params = []

    if args.sequences != None: # run tests
        p = partition( args.sequences, circle = args.circle, params = args.parameters, verbose = args.verbose, mfe = args.mfe, calc_bpp = args.bpp, n_stochastic = int(args.stochastic), do_enumeration = args.enumerate, structure = args.structure, force_base_pairs = args.force_base_pairs, deriv_params = args.deriv_params, no_coax = args.no_coax, use_simple_recursions = args.simple )
    else:
        call_all_tests( verbose = args.verbose, use_simple_recursions = args.simple )
