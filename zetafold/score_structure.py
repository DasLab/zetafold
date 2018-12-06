#!/usr/bin/python
from __future__ import print_function
import argparse
from math import log
import sys
import os
if __package__ == None: sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from zetafold.partition import partition
from zetafold.parameters import get_params
from zetafold.util import sequence_util
from zetafold.util import secstruct_util
from zetafold.util.assert_equal import assert_equal
from zetafold.util.constants import KT_IN_KCAL
from zetafold.util.output_util import show_derivs

def score_structure( sequences, structure, circle = False, params = None, test_mode = False, deriv_params = None ):

    # What we get if we parse out motifs
    structure = secstruct_util.get_structure_string( structure )
    bps_list  = secstruct_util.bps_from_secstruct( structure )
    motifs = secstruct_util.parse_motifs( structure )
    sequence, ligated, sequences = sequence_util.initialize_sequence_and_ligated( sequences, circle )
    params = get_params( params, suppress_all_output = True )
    Kd_ref = params.base_pair_types[0].Kd # Kd[G-C], a la Turner rule convention
    C_std  = params.C_std

    # Now go through each motif parsed out of the target structure
    Z = 1.0
    log_derivs = None
    for motif in motifs:
        motif_res = []
        motif_sequences = []
        for strand in motif:
            strand_sequence = ''
            for i in strand:
                motif_res.append( i )
                strand_sequence += sequence[i]
                if not ligated[i]:
                    motif_sequences.append( strand_sequence )
                    strand_sequence = ''
            if len( strand_sequence ) > 0: motif_sequences.append( strand_sequence )
        motif_circle = ligated[ motif_res[-1] ] and ( (motif_res[0] - motif_res[-1]) % len(sequence) == 1 )
        motif_sequence = ''.join( motif_sequences )

        # each motif res better show up only once
        assert( len( set( motif_res ) ) == len( motif_res ) )

        motif_bps_list = []
        for i,j in bps_list:
            if motif_res.count( i ) == 0: continue
            if motif_res.count( j ) == 0: continue
            motif_bps_list.append( (motif_res.index(i), motif_res.index(j)) )
        motif_structure = secstruct_util.secstruct_from_bps( motif_bps_list, len( motif_res ) )

        p = partition( motif_sequences, circle = motif_circle, structure = motif_structure, params = params, suppress_all_output = True, deriv_params = deriv_params )
        Z_motif = p.Z
        log_derivs_motif = p.log_derivs

        # Need to 'correct' for half-terminal penalties (a la Turner rules) and also remove extra costs
        # for connecting these 'sub-strands' together.

        Z_motif *= ( Kd_ref / C_std ) ** sequence_util.get_num_strand_connections( motif_sequences, motif_circle )
        for i_motif, j_motif in motif_bps_list:
            # what kind of base pair is this?
            for base_pair_type in p.params.base_pair_types:
                if base_pair_type.is_match( motif_sequence[ i_motif ], motif_sequence[ j_motif ] ):
                    Z_motif *= ( base_pair_type.Kd / Kd_ref )**(0.5)
                    if deriv_params:
                        Kd_tags = ["Kd_"+base_pair_type.get_tag(), "Kd_"+base_pair_type.flipped.get_tag()]
                        for Kd_tag in Kd_tags:
                            if deriv_params.count( Kd_tag ): log_derivs_motif[ deriv_params.index( Kd_tag ) ] += 0.5
                    break

        if test_mode: print("Motif: ", Z_motif, motif, motif_sequences, motif_structure)

        Z *= Z_motif

        if deriv_params:
            if log_derivs == None: log_derivs = [0.0]*len( deriv_params )
            for n, log_deriv_motif in enumerate( log_derivs_motif ): log_derivs[n] += log_deriv_motif

    # Compute cost of connecting the strands into a complex
    Z_connect = ( C_std / Kd_ref ) ** sequence_util.get_num_strand_connections( sequences, circle )
    Z *= Z_connect

    if deriv_params and log_derivs == None: log_derivs = [0.0]*len( deriv_params )

    if test_mode:
        print("Connect strands: ", Z_connect)
        print('Product of motif Z      :', Z)

    if test_mode:
        # Reference value from 'hacked' dynamic programming, which takes a while.
        p = partition( sequences, circle = circle, structure = structure, params = params, suppress_all_output = True, deriv_params = deriv_params )
        print('From dynamic programming:', p.Z)
        assert_equal( Z, p.Z )
        print()
        if deriv_params:
            print( 'LOG-DERIVS FROM MOTIF-BY_MOTIF' )
            show_derivs( deriv_params, log_derivs )
            print( 'LOG-DERIVS FROM FULL PARTITION' )
            show_derivs( deriv_params, p.log_derivs )

    dG = -KT_IN_KCAL * log( Z )
    if deriv_params: return (dG,log_derivs)
    return dG

if __name__=='__main__':
    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("-struct","--structure",type=str, default=None, help='force structure in dot-parens notation')
    parser.add_argument("-c","-circ","--circle", action='store_true', default=False, help='Sequence is a circle')
    parser.add_argument("-params","--parameters",type=str, default='', help='Parameters to use [default: '']')
    parser.add_argument("-test","--test_mode",action="store_true", default=False, help='In test mode, also run (slow) dynamic programming calculation to get Z' )
    parser.add_argument("--calc_deriv", action='store_true', default=False, help='Calculate derivative with respect to all parameters')
    parser.add_argument( "--deriv_params",help="Parameters for which to calculate derivatives. Default: None, or all params if --calc_deriv",nargs='*')

    args     = parser.parse_args()
    if args.calc_deriv and args.deriv_params == None: args.deriv_params = []

    dG = score_structure( args.sequences, args.structure, circle = args.circle, params = args.parameters, test_mode = args.test_mode, deriv_params = args.deriv_params )
    if args.deriv_params: (dG,log_derivs) = dG
    print('dG = ',dG)
