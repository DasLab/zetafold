#!/usr/bin/python
from __future__ import print_function
from zetafold.parameters import get_params
from zetafold.data.training_examples import *
from zetafold.training import *
from scipy.optimize import minimize
import numpy as np
from multiprocessing import Pool
import math
import argparse
import __builtin__

parser = argparse.ArgumentParser( description = "Test nearest neighbor model partitition function for RNA sequence" )
parser.add_argument("-params","--parameters", type=str, help='Parameter file to use [default: use latest zetafold version]')
parser.add_argument("--train_data", type=str, help="Training data to use. Give none to get list.")
parser.add_argument("--train_params", help="Parameters to optimize. Give none to get list.", nargs='*')
parser.add_argument("--train_params_exclude", help="Parameters to optimize. Give none to get list.", nargs='*')
parser.add_argument("--jobs","-j", type=int, default=4, help='Number of jobs to run in parallel')
parser.add_argument("--outfile","-out","-o", type=str, help='Outfile to save loss/variables during training')
parser.add_argument("--final_params_file",type=str,default='final.params',help="Name of params file for outputting final",nargs='*')
parser.add_argument("--use_derivs","-d", action='store_true', help='Use analytical derivatives during training')
parser.add_argument("--init_params",help="Initial values for parameters (default are values in params file)",nargs='*')
parser.add_argument("--init_log_params",help="Initial values for log parameters (alternative to init_params)",nargs='*')
parser.add_argument("--no_coax", action='store_true', default=False, help='Turn off coaxial stacking')
parser.add_argument("--deriv_check", action='store_true', default=False, help='Run numerical vs. analytical deriv check')
parser.add_argument("-eval","--evaluate", action='store_true', default=False, help='Just evaluate at starting parameter values')
parser.add_argument("--allow_extra_base_pairs",action='store_true',default=False, help='allow extra base pairs compatible with --structure')
parser.add_argument("--use_priors",action='store_true', help='add priors to force log parameters to stay in reasonable bounds.')
parser.add_argument("--use_bounds",action='store_true', help='force log parameters to stay in reasonable bounds; not applied to BFGS')
parser.add_argument("--method",type=str,default='BFGS',help="Minimization routine")
args     = parser.parse_args()

# pool of CPU's to use. for testing on local machines, specify -j1 to get builtin CPU -- allows ctrl-c to cancel.
pool = __builtin__
if args.jobs > 1: pool = Pool( args.jobs )

# set up parameter file
params = get_params( args.parameters, suppress_all_output = True )
if args.no_coax: params.set_parameter( 'K_coax', 0.0 )

# set up training examples
training_examples = initialize_training_examples( all_training_examples, training_sets, training_set_names, args.train_data )
train_parameters  = initialize_train_parameters( params, args.train_params, args.train_params_exclude, args.no_coax )
x0                = initialize_parameter_values( params, train_parameters, args.init_params, args.init_log_params, (args.use_bounds or args.use_priors) )
bounds = get_bounds( train_parameters ) if args.use_bounds else None
priors = get_priors( train_parameters ) if args.use_priors else None

loss = lambda x:free_energy_gap(      x,params,train_parameters,training_examples,args.allow_extra_base_pairs,priors,pool,args.outfile)
grad = lambda x:free_energy_gap_deriv(x,params,train_parameters,training_examples,args.allow_extra_base_pairs,priors,pool)
jac = grad if args.use_derivs else None

if args.deriv_check: train_deriv_check( x0, loss, grad, train_parameters )
if args.evaluate:
    loss( x0 )
    exit(0)

create_outfile( args.outfile, params, train_parameters )
result = minimize( loss, x0, method = args.method, jac = jac, bounds = bounds )
final_loss = result.fun

print(result)
print('Final parameters:', result.x, 'Loss:',final_loss )

pack_variables( result.x, params, train_parameters )
params.output_to_file( args.final_params_file )
