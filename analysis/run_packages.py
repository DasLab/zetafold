#!/usr/bin/python
import argparse, os
from zetafold.data.training_examples import *
from zetafold.training import *
from multiprocessing import Pool
import __builtin__

parser = argparse.ArgumentParser( description = "Run nearest neighbor modeling package to get base pairing probability matrix for RNA sequence" )
parser.add_argument("-p","--packages", type=str, help='Packages to use -- give zetafold parameters file, or contrafold...',nargs='*')
parser.add_argument("--data", type=str, help="Data to use. Give none to get list.")
parser.add_argument("-f","--force", action='store_true', help="Overwrite prior output.")
parser.add_argument("--jobs","-j", type=int, default=1, help='Number of jobs to run in parallel')
args     = parser.parse_args()

examples = initialize_training_examples( all_training_examples, training_sets, training_set_names, args.data )

def make_fasta_file( subdirname, example ):
    fasta = '%s/%s.fasta' % (subdirname,example.name)
    fid = open( fasta, 'w' )
    fid.write( '>%s\n' % example.name )
    fid.write( '%s\n' % example.sequence )
    fid.close()
    return fasta

def run_package( example ):
    package = example.package
    dirname = os.path.basename(package).replace( '.params', '' )
    assert( os.path.exists( dirname ) )
    subdirname = dirname+'/'+example.name

    if not os.path.exists( subdirname ):
        print 'Creating directory: ', subdirname
        os.mkdir( subdirname )

    if dirname == 'contrafold' or dirname == 'contrafold-nc':
        fasta = make_fasta_file( subdirname, example )
        noncomplementaryflag = '--noncomplementary' if dirname == 'contrafold-nc' else ''
        cmdline = 'contrafold predict %s %s --parens %s/secstruct.txt --posteriors 0.00001 %s/posteriors.txt > %s/contrafold.out 2> %s/contrafold.err' % (fasta,noncomplementaryflag,subdirname,subdirname,subdirname,subdirname)
        outfile = '%s/contrafold.out' % subdirname
        if not args.force and os.path.exists( '%s/posteriors.txt' % subdirname ): return
    elif dirname == 'vienna':
        fasta = make_fasta_file( subdirname, example )
        cmdline = 'RNAfold -p %s >  %s/vienna.out 2> %s/vienna.err; mv %s_dp.ps %s' % (fasta,subdirname,subdirname,example.name,subdirname)
        outfile = '%s/vienna.out' % subdirname
        if not args.force and os.path.exists( '%s/%s_dp.ps' % (subdirname,example.name) ): return
    elif dirname == 'nupack' or dirname == 'nupack95':
        infile = '%s/%s.in' % (subdirname,example.name)
        fid = open( infile, 'w' )
        fid.write('%s\n' % example.sequence)
        fid.close()
        params_flag = '-material rna1999'
        if dirname == 'nupack95': params_flag = '-material rna1995'
        cmdline = 'pairs %s/%s %s --cutoff 0.00001 > %s/nupack.out 2> %s/nupack.err' % (subdirname,example.name,params_flag,subdirname,subdirname)
        outfile = '%s/nupack.out' % subdirname
        if not args.force and os.path.exists( '%s/%s.ppairs' % (subdirname,example.name) ): return
    elif dirname == 'rnastructure':
        infile = '%s/%s.seq' % (subdirname,example.name)
        fid = open( infile, 'w' )
        fid.write( ';\n%s\n\n%s1\n' % (example.name,example.sequence) )
        fid.close()
        cmdline = 'partition %s/%s.seq %s/%s.pfs > %s/rnastructure.out 2> %s/rnastructure.err; ProbabilityPlot %s/%s.pfs -T %s/probability_plot.txt' % (subdirname,example.name,subdirname,example.name,subdirname,subdirname,subdirname,example.name,subdirname)
        outfile = '%s/rnastructure.out' % subdirname
        if not args.force and os.path.exists( '%s/probability_plot.txt' % (subdirname) ): return
    else:
        cmdline = 'zetafold.py -s %s -params %s --bpp --stochastic 100 --mfe --calc_gap_structure "%s" --bpp_file  %s/bpp.txt.gz > %s/zetafold.out 2> %s/zetafold.err' % (example.sequence,package,example.structure,subdirname,subdirname,subdirname)
        outfile = '%s/zetafold.out' % subdirname
        if not args.force and os.path.exists( '%s/bpp.txt.gz' % subdirname ): return

    print cmdline
    os.system( cmdline )
    os.system( 'cat %s' % outfile  )

if args.packages == None:
    print 'Specify --packages. Available packages: '
    for package in ['nupack','nupack95','vienna','rnastructure','contrafold','contrafold-nc','zetafold_v0.18','Any zetafold parameter file']: print ' ',package
    exit()

pool = __builtin__
if args.jobs > 1: pool = Pool( args.jobs )

for package in args.packages:
    dirname = os.path.basename(package).replace( '.params', '' )
    if not os.path.exists(dirname):
        print 'Creating directory: ', dirname
        os.mkdir( dirname )

    for example in examples: example.package = package

    executable = 'zetafold.py'
    if package.count( 'contrafold' ): executable = 'contrafold'
    if package.count( 'nupack' ): executable = 'pairs'
    if package.count( 'vienna' ): executable = 'RNAfold'
    if package.count( 'rnastructure' ): executable = 'partition'
    if os.system( 'which '+executable+' > /dev/null' ) != 0:
        print '\nFor',package,', you need the executable ',executable,'in your path!\n'
        exit()

    pool.map( run_package, examples )

