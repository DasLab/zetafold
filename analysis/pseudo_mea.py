from mea import MEA
import numpy as np
from glob import glob
import argparse
import sys, os

def pseudoacc(matrix_list, true_struct_list, gamma_min=-7, gamma_max=7, verbose=True):
    '''Compute MEA over a range of gamma values to find optimum pseudoaccuracy.  

    Inputs:
    Important! Files in matrix_dir and true_structs need to have the same names corresponding
        to their same constructs, but suffixes don't matter.
    matrix_dir: list of NxN base pair probability matrices.
    true_structs: list of NxN true structure base pair matrices. Can be
        symmetric matrices or not; upper triangle is taken.
    gamma_min, gamma_max: minimum/maximum log_2(gamma) value used.
    verbose: print output or not (for command line use)

    Outputs:
    SEN: TP/(TP+FN), library keyed by gamma values used.
    PPV: TP/(TP+FP), "
    MCC: Mathews correlation coefficient
    Fscore: 2*TP/(2*TP + FP + FN)    
    '''

    matrices, structs = [], []

    if len(matrix_list) == 0:
        raise ValueError('No matrix files found!')

    for x in matrix_list:
        matrices.append(np.loadtxt(x))
        for s in true_struct_list:
            if os.path.basename(x).split('.')[0] in s:
                structs.append(np.loadtxt(s))

    gamma_vals = [x for x in range(gamma_min, gamma_max)]
    TP = {k:1e-6 for k in gamma_vals}
    FP = {k:1e-6 for k in gamma_vals}
    cFP = {k:1e-6 for k in gamma_vals}
    TN = {k:1e-6 for k in gamma_vals}
    FN = {k:1e-6 for k in gamma_vals}

    SEN = {}
    PPV = {}
    MCC = {}
    Fscore = {}

    for g in gamma_vals:
        for i in range(len(matrices)):
            cls = MEA(matrices[i], gamma=2**g)
            tp,fp,cfp, tn,fn = cls.score(structs[i][:cls.N])
            TP[g] += tp
            FP[g] += fp
            cFP[g] += cfp
            TN[g] += tn
            FN[g] += fn

        SEN[g] = TP[g]/(TP[g] + FN[g])
        PPV[g] = TP[g]/(TP[g] + FP[g]-cFP[g])
        MCC[g] = (TP[g]*TN[g] - (FP[g]-cFP[g])*FN[g])/np.sqrt((TP[g]+FP[g]-cFP[g])*(TP[g]+FN[g])*(TN[g]+FP[g]-cFP[g])*(TN[g]+FN[g]))
        Fscore[g] = 2*TP[g]/(2*TP[g] + FP[g] - cFP[g] + FN[g])

        if verbose: print("%d\t%.3f\t%.3f\t%.3f\t%.3f" % (g, SEN[g], PPV[g], MCC[g], Fscore[g]))
    return SEN, PPV, MCC, Fscore

if __name__ == '__main__':
    parser=argparse.ArgumentParser(
        description='''pseudoaccuracy structure prediction, using CentroidFold algorithm (Hamada 2009).  Important: base pair probability matrices (specified in --bp_matrices) need to have same names as structure matrices (specified in --true_structs), but the extensions don't matter.''')
    parser.add_argument('--bp_matrices','-p', nargs='+', help='path to NxN matrices of bp probabilities, i.e. `contrafold/*.bpps`.')
    parser.add_argument('--true_structs','-s', nargs='+', help='path to NxN matrices of true structures, i.e. `rnaview/*.struct`.  These can be symmetric or not -- upper triangle is taken.')
    parser.add_argument('--gamma_min',type=int, default=-7, help='minimum value for log_2(gamma)')
    parser.add_argument('--gamma_max',type=int, default=7, help='maximum value for log_2(gamma)')
    args = parser.parse_args()
    print('Pseudoaccuracy scan from MEA structure')
    print('Path to first base pair matrix: %s' % args.bp_matrices[0])
    print('Path to first true struct matrix: %s' % args.true_structs[0])
    print("log_2(gamma)\tSEN\tPPV\tMCC\tFscore")
    pseudoacc(args.bp_matrices, args.true_structs, gamma_min = args.gamma_min, gamma_max = args.gamma_max, verbose=True)

