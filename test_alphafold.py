#!/usr/bin/python
import argparse
from alphafold import partition, C_init, l, C_init_BP, Kd_BP, l_BP, C_std #output_DP, output_square, partition

def output_test( Z, Z_ref = 0, bpp = [], bpp_idx= [], bpp_expected = 0):
    print('Z =',Z_ref,' [expected]')
    assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
    print()
    print('bpp[0,4] = ',bpp[ bpp_idx[0] ][ bpp_idx[1] ])
    print('bpp[0,4] = ',bpp_expected,' [expected]')
    assert( abs( (bpp[ bpp_idx[0] ][ bpp_idx[1] ] - bpp_expected)/bpp[ bpp_idx[0] ][ bpp_idx[1] ] )  < 1e-5 )
    print()


def test_sequences():
    # test of sequences where we know the final partition function.
    sequence = 'CAAAGAA'
    (Z, bpp) = partition( sequence, circle = True )
    output_test( Z, C_init  * (l**7) * (1 + C_init_BP / Kd_BP ) / C_std, \
                 bpp, [0,4], (C_init**2 * (l**3) * l_BP/ Kd_BP) / ( C_init * (l**4) + C_init**2 * (l**3) * l_BP/ Kd_BP) )

    sequence = 'CAG'
    (Z, bpp) = partition( sequence )
    output_test( Z, 1 + C_init * l**2 / Kd_BP, \
                 bpp, [0,2], (C_init * l**2/Kd_BP)/( 1 + C_init * l**2/Kd_BP ) )

    sequences = ['C','G']
    (Z, bpp) = partition( sequences ) # note that Z sums over only base pair (not dissociated strands!)
    output_test( Z, C_std * l / Kd_BP/l_BP, \
                 bpp, [0,1], 1.0 )

    sequences = ['GC','GC']
    (Z, bpp) = partition( sequences )
    output_test( Z, (C_std/Kd_BP)*(l/l_BP)*(2 + l*l_BP*C_init/Kd_BP ), \
                 bpp, [0,3], (1 + l*l_BP*C_init/Kd_BP )/(2 + l*l_BP*C_init/Kd_BP ) )

    sequence = 'CAGGC'
    (Z, bpp) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
    output_test( Z, 1+C_init*l**2/Kd_BP * ( 2 + l ), \
                 bpp, [0,2], C_init*l**2/Kd_BP /(  1+C_init*l**2/Kd_BP * ( 2 + l )) )

