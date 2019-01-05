#!/usr/bin/python
import argparse, os
from zetafold.data.training_examples import *
from zetafold.training import *
from zetafold.util.secstruct_util import *
from multiprocessing import Pool
import gzip
import __builtin__

parser = argparse.ArgumentParser( description = "Test nearest neighbor model partitition function for RNA sequence" )
parser.add_argument("-p","--packages", type=str, help='Packages to use (must correspond to subdirectories)',nargs='*')
parser.add_argument("--data", type=str, help="Data to use. Give none to get list.")
args     = parser.parse_args()

examples = initialize_training_examples( all_training_examples, training_sets, training_set_names, args.data )

def read_bpp_file( bpp_file ):
    if len( bpp_file ) > 3 and bpp_file[-3:] == '.gz':
        lines = gzip.open( bpp_file ).readlines()
    else:
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

def read_ppairs_file( ppairs_file ):
    lines = open( ppairs_file ).readlines()
    in_data = False
    for line in lines:
        if len( line ) < 2: continue
        if len( line ) > 2 and line[0] != '%' and not in_data:
            in_data = True
            N = int( line[:-1] )
            p_unpaired = [0.0]*N
            bpp = [None]*N
            for i in range( N ): bpp[ i ] = [0.0]*N
            continue
        if not in_data: continue
        cols = line[:-1].split()
        i = int( cols[0] )
        j = int( cols[1] )
        val = float( cols[2] )
        if ( j == N+1 ):
            p_unpaired[i-1] = val
        else:
            assert( j <= N )
            bpp[i-1][j-1] = val
            bpp[j-1][i-1] = val
    p_paired = [sum( x ) for x in bpp]
    p_tot = [ x[0]+x[1] for x in zip( p_paired, p_unpaired ) ]
    for p in p_tot: assert( abs(1.0-p)<1.0e-3 ) # Sanity Check
    return bpp

def read_dp_ps_file( dp_ps_file ):
    lines = open( dp_ps_file ).readlines()
    sequence_line = False
    sequence = ''
    for line in lines:
        if len( line ) < 2: continue
        if len( line ) > 2 and line[:9] == '/sequence':
            sequence_line = True
            continue
        if sequence_line:
            sequence += line[:-2]
            if line.count( '}' ):
                sequence_line = False
                N = len( sequence )
                bpp = [None]*N
                for i in range( N ): bpp[ i ] = [0.0]*N
        if line[0] == '\%': continue
        cols = line[:-1].split()
        if len( cols ) != 4: continue
        if cols[-1] != 'ubox': continue
        i = int( cols[0] )
        j = int( cols[1] )
        val = float( cols[2] )**2
        bpp[i-1][j-1] = val
        bpp[j-1][i-1] = val
    return bpp

def read_probability_plot_file( probability_plot_file ):
    lines = open( probability_plot_file ).readlines()
    N = int( lines[0][:-1] )
    bpp = [None]*N
    for i in range( N ): bpp[ i ] = [0.0]*N
    for line in lines[2:]:
        cols = line[:-1].split()
        i = int( cols[0] )
        j = int( cols[1] )
        val = 10.0 ** (-float( cols[2] ))
        bpp[i-1][j-1] = val
        bpp[j-1][i-1] = val
    return bpp

def read_simfold_file( simfold_file ):
    if len( simfold_file ) > 3 and simfold_file[-3:] == '.gz':
        lines = gzip.open( simfold_file ).readlines()
    else:
        lines = open( simfold_file ).readlines()
    N = len( lines[0][:-1].split()[-1] )
    bpp = [None]*N
    for i in range( N ): bpp[ i ] = [0.0]*N
    for line in lines[5:]:
        cols = line[:-1].split()
        i = int( cols[0] )
        j = int( cols[1] )
        val = float( cols[2] )
        bpp[i][j] = val
        bpp[j][i] = val
    return bpp

ensemble_defects = {}
bpp_defects = {}
for package in args.packages:
    ensemble_defects[ package ] = {}
    bpp_defects[ package ] = {}
    for example in examples:
        subdirname = package+'/'+example.name
        bpp = None
        readin_func = { 'bpp.txt':read_bpp_file,'bpp.txt.gz':read_bpp_file,\
                        'posteriors.txt':read_posteriors_file, '%s.ppairs' % example.name : read_ppairs_file,
                        '%s_dp.ps' % example.name : read_dp_ps_file, 'probability_plot.txt': read_probability_plot_file,
                        'simfold_pf.out' : read_simfold_file, 'simfold_pf.out.gz' : read_simfold_file }
        for bpp_file in readin_func.keys():
            if os.path.isfile( subdirname + '/' + bpp_file ):  bpp = readin_func[bpp_file]( subdirname + '/' + bpp_file )
        if bpp == None:
            print 'COULD NOT FIND ',readin_func.keys(),' for: ', example.name, ' with package ', package
            print ' Did you run run_packages.py --data ', args.data, ' --packages ', package, '?'
            exit()
        bps = bps_from_secstruct( example.structure )
        bpp_defect = 0.0
        for (i,j) in bps: bpp_defect += 1.0 - bpp[i][j]
        bpp_defects[ package ][ example ] = bpp_defect/len(bps)

        ensemble_defect = bpp_defect * 2.0
        N = 2* len(bps)
        bpp_sum = [sum(x) for x in bpp]
        for s in bpp_sum: assert( s >= 0.0 and s <= 1.001 )
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

# each RNA will be a row.
print_table( ensemble_defects, 'Ensemble defect' )
print_table( bpp_defects, 'BPP defect' )


# Summary table, just overall summary vals
print 'Overall summary'
print '%60s %20s %20s' % ('Package','Ensemble defect','BPP defect')
for package in args.packages:
    print '%60s %20.4f %20.4f' % (package, ensemble_defects[ package ]['overall'], bpp_defects[ package ]['overall'])
