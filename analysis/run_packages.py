#!/usr/bin/python
import argparse, os
from zetafold.data.training_examples import *
from zetafold.training import *
from multiprocessing import Pool
import __builtin__

parser = argparse.ArgumentParser( description = "Run nearest neighbor modeling package to get base pairing probability matrix for RNA sequence" )
parser.add_argument("-p","--packages", type=str, help='Packages to use -- could be zetafold .params files, or contrafold',nargs='*')
parser.add_argument("--data", type=str, help="Data to use. Give none to get list.")
parser.add_argument("-f","--force", action='store_true', help="Overwrite prior output.")
parser.add_argument("--jobs","-j", type=int, default=1, help='Number of jobs to run in parallel')
args     = parser.parse_args()

examples = initialize_training_examples( all_training_examples, training_sets, training_set_names, args.data )

def run_package( example ):
    package = example.package
    dirname = os.path.basename(package).replace( '.params', '' )
    assert( os.path.exists( dirname ) )
    subdirname = dirname+'/'+example.name

    if not os.path.exists( subdirname ):
        print 'Creating directory: ', subdirname
        os.mkdir( subdirname )

    if dirname == 'contrafold' or dirname == 'contrafold-nc':
        fasta = '%s/%s.fasta' % (subdirname,example.name)
        fid = open( fasta, 'w' )
        fid.write( '>%s\n' % example.name )
        fid.write( '%s\n' % example.sequence )
        fid.close()
        noncomplementaryflag = '--noncomplementary' if dirname == 'contrafold-nc' else ''
        cmdline = 'contrafold predict %s %s --parens %s/secstruct.txt --posteriors 0.00001 %s/posteriors.txt > %s/contrafold.out 2> %s/contrafold.err' % (fasta,noncomplementaryflag,subdirname,subdirname,subdirname,subdirname)
        outfile = '%s/contrafold.out' % subdirname
        if not args.force and os.path.exists( '%s/posteriors.txt' % subdirname ): return
    else:
        cmdline = 'zetafold.py -s %s -params %s --bpp --stochastic 100 --mfe --calc_gap_structure "%s" --bpp_file  %s/bpp.txt.gz > %s/zetafold.out 2> %s/zetafold.err' % (example.sequence,package,example.structure,subdirname,subdirname,subdirname)
        outfile = '%s/zetafold.out' % subdirname
        if not args.force and os.path.exists( '%s/bpp.txt.gz' % subdirname ): return

    print cmdline
    os.system( cmdline )
    os.system( 'cat %s' % outfile  )

pool = __builtin__
if args.jobs > 1: pool = Pool( args.jobs )

for package in args.packages:
    dirname = os.path.basename(package).replace( '.params', '' )
    if not os.path.exists(dirname):
        print 'Creating directory: ', dirname
        os.mkdir( dirname )

    for example in examples: example.package = package

    executable = 'contrafold' if package.count( 'contrafold' ) else 'zetafold.py'
    if os.system( 'which '+executable+' > /dev/null' ) != 0:
        print '\nHey you need the executable ',executable,'in your path!\n'
        exit()

    pool.map( run_package, examples )

