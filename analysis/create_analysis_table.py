#!/usr/bin/python
import argparse, os
from zetafold.data.training_examples import *
from zetafold.training import *
from zetafold.util.secstruct_util import *
from multiprocessing import Pool
import __builtin__

parser = argparse.ArgumentParser( description = "Test nearest neighbor model partitition function for RNA sequence" )
parser.add_argument("-p","--packages", type=str, help='Packages to use (must correspond to subdirectories)',nargs='*')
parser.add_argument("--data", type=str, help="Data to use. Give none to get list.")
args     = parser.parse_args()

examples = initialize_training_examples( all_training_examples, training_sets, training_set_names, args.data )


def read_bpp_file( bpp_file ):
    lines = open(bpp_file).readlines()
    bpp = []
    for line in lines:
        bpp_line = [ float(val) for val in line[:-1].split() ]
        bpp.append( bpp_line )
    return bpp

def read_posteriors_file( posteriors_file ):
    lines = open( posteriors_file ).readlines()
    N = len( lines )
    bpp = [None]*N
    for i in range( N ): bpp[ i ] = [0.0]*N
    for n,line in enumerate(lines):
        cols = line[:-1].split()
        i = int( cols[0] )
        assert( i == n + 1)
        for col in cols[2:]:
            tags = col.split(':')
            j = int( tags[0] )
            val = float( tags[1] )
            bpp[i-1][j-1] = val
            bpp[j-1][i-1] = val
    return bpp

ensemble_defects = {}
bpp_defects = {}
for package in args.packages:
    ensemble_defects[ package ] = {}
    bpp_defects[ package ] = {}
    for example in examples:
        subdirname = package+'/'+example.name
        if os.path.isfile( subdirname + '/bpp.txt' ):
            bpp = read_bpp_file( subdirname+'/bpp.txt' )
        elif package == 'contrafold' or package == 'contrafold-nc':
            assert( os.path.isfile( subdirname+'/posteriors.txt' ) )
            bpp = read_posteriors_file( subdirname+'/posteriors.txt'  )
        else:
            print 'COULD NOT FIND bpp.txt for: ', example.name, ' with package ', package
            exit()
        bps = bps_from_secstruct( example.structure )
        bpp_defect = 0.0
        for (i,j) in bps: bpp_defect += 1.0 - bpp[i][j]
        bpp_defects[ package ][ example ] = bpp_defect/len(bps)

        ensemble_defect = bpp_defect * 2.0
        N = 2* len(bps)
        bpp_sum = [sum(x) for x in bpp]
        for (i,char) in enumerate( example.structure ):
            if char == '.':
                N += 1
                ensemble_defect += bpp_sum[i]
        assert( N == len( example.structure ) )
        ensemble_defects[ package ][ example ] = ensemble_defect / N

# calculate overall
for package in args.packages:
    bpp_defects[ package ][ 'overall' ] = 0.0
    ensemble_defects[ package ][ 'overall' ] = 0.0
    tot_bps = 0
    tot_nts = 0
    for example in examples:
        num_bps = len( bps_from_secstruct( example.structure ) )
        tot_bps += num_bps
        tot_nts += len( example.structure )
        if example not in bpp_defects[ package ].keys():
            print 'HEY PROBLEM! %s not in package %s.' % (example.name,package)
            exit()
        bpp_defects[ package ][ 'overall' ]      += bpp_defects[ package ][ example ] * num_bps
        ensemble_defects[ package ][ 'overall' ] += ensemble_defects[ package ][ example ] * len( example.structure )
    bpp_defects[ package ][ 'overall' ] /= tot_bps
    ensemble_defects[ package ][ 'overall' ] /= tot_nts

def truncate( string, N ):
    if len( string ) > N: return string[:N]
    return string

def print_table_header():
    print '%30s' % 'RNA',
    for package in args.packages:
        print '%15s' % truncate(package,15),
    print

def print_table( data, tag ):
    print tag
    print_table_header()
    for example in examples:
        print '%30s' % example.name,
        for package in args.packages:
            print '%15.4f' % data[package][example],
        print
    print '%30s' % 'OVERALL',
    for package in args.packages:
        print '%15.4f' % data[ package ]['overall'],
    print
    print_table_header()
    print

print_table( ensemble_defects, 'Ensemble defect' )
print_table( bpp_defects, 'BPP defect' )
