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
        cmdline = 'RNAfold -p %s >  %s/vienna.out 2> %s/vienna.err; mv %s_ss.ps %s_dp.ps %s' % (fasta,subdirname,subdirname,example.name,example.name,subdirname)
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
    elif dirname.count( 'rnasoft' ) > 0:
        exe = os.popen( 'which simfold_pf' ).readlines()[0]
        print exe
        if dirname == 'rnasoft99':          parameter_file = os.path.dirname( exe )+'/params/turner_parameters_fm363_constrdangles.txt'
        elif dirname == 'rnasoft07':        parameter_file = os.path.dirname( exe )+'/params/CG_best_parameters_ISMB2007.txt'
        elif dirname == 'rnasoftBLFRstar': # THIS IS GIVING WACKY RESULTS
            if 'ANDRONESCU_DATA' not in os.environ:
                print 'Need to defined ANDRONESCU DATA through a line like\nexport ANDRONESCU_DATA=$HOME/projects/RNA/zetafold_tests/other_packages/andronescu/Andronescu2007_results_data_script/data\n in your .bashrc'
            parameter_file = os.environ[ 'ANDRONESCU_DATA' ] +'/parameters_BLFRstar.txt'
        elif dirname == 'rnasoftBLstar':
            parameter_file = os.environ[ 'ANDRONESCU_DATA' ] +'/parameters_BLstar.txt'
        else:
            print 'Need to specify rnasoft99, rnasoft77, rnasoftBLFRstar'
        parameter_flag = '-p %s' % parameter_file
        cmdline = '`which simfold_pf` -s "%s" %s > %s/simfold_pf.out 2> %s/simfold_pf.err; gzip -f %s/simfold_pf.out' % (example.sequence,parameter_flag,subdirname,subdirname,subdirname)
        outfile = '%s/simfold_pf.err' % subdirname
        if not args.force and os.path.exists( '%s/simfold_pf.out.gz' ): return
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
    if package.count( 'rnasoft' ): executable = 'simfold_pf'
    if os.system( 'which '+executable+' > /dev/null' ) != 0:
        print '\nFor %s, you need the executable %s in your path!\n' % (package,executable)
        exit()

    pool.map( run_package, examples )

