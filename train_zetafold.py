#!/usr/bin/python
from zetafold.parameters import get_params
from zetafold.data.training_examples import *
from zetafold.training import *
from scipy.optimize import minimize
import numpy as np
from multiprocessing import Pool
import math
import argparse

parser = argparse.ArgumentParser( description = "Test nearest neighbor model partitition function for RNA sequence" )
parser.add_argument("-params","--parameters",type=str, default='', help='Parameter file to use [default: '', which triggers latest version]')
parser.add_argument( "--train_data",type=str,help="Training data to use. Give none to get list.")
parser.add_argument( "--train_params",help="Parameters to optimize. Give none to get list.",nargs='*')
parser.add_argument( "--train_params_exclude",help="Parameters to optimize. Give none to get list.",nargs='*')
parser.add_argument( "--init_params",help="Initial values for parameters",nargs='*')
parser.add_argument( "--init_log_params",help="Initial values for log parameters (alternative to init_params)",nargs='*')
parser.add_argument( "--method",type=str,default='BFGS',help="Minimization routine")
parser.add_argument("--no_coax", action='store_true', default=False, help='Turn off coaxial stacking')
parser.add_argument("--jobs","-j", type=int, default=4, help='Number of jobs to run in parallel')
parser.add_argument("--outfile","-out","-o", type=str, help='Outfile to save loss/variables during training')
parser.add_argument("--use_derivs","-d", action='store_true', help='Use analytical derivatives during training')
args     = parser.parse_args()

# set up parameter file
params = get_params( args.parameters, suppress_all_output = True )
if args.no_coax: params.set_parameter( 'K_coax', 0.0 )

# set up training examples
if args.train_data == None:
    print '\nMust specify training set. Options are:'
    for set_name in training_sets.keys():
        print '%30s:' % set_name,
        for training_example in training_sets[ set_name ]: print training_example.name,
        print
    exit()

training_examples = training_sets[ args.train_data ] #, P4P6_outerjunction ]
train_parameters = args.train_params
if train_parameters == None:
    if args.train_params_exclude:
        train_parameters = []
        for param_exclude in args.train_params_exclude:  assert( param_exclude in params.parameter_tags )
        for param in params.parameter_tags:
            if not param in args.train_params_exclude: train_parameters.append( param )
    else:
        print '\nMust specify which parameters to optimize'
        params.show_parameters()
        exit()

x0 = None
if args.init_log_params: x0 = [float(log_param) for log_param in args.init_log_params ]
elif args.init_params:   x0 = [ np.log(float(param)) for param in  args.init_params]
else:
    x0 = np.zeros( len(train_parameters) )
    for n,param_tag in enumerate(train_parameters):
        val = params.get_parameter_value( param_tag )
        if val == 0.0:  x0[ n ] = -5.0
        else: x0[ n ] = np.log( val )

pool = Pool( args.jobs )

loss = lambda x:free_energy_gap(      x,params,train_parameters,training_examples,pool,args.outfile)
grad = lambda x:free_energy_gap_deriv(x,params,train_parameters,training_examples,pool)
jac = grad if args.use_derivs else None

create_outfile( args.outfile, params, train_parameters )
result = minimize( loss, x0, method = args.method, jac = jac )
final_loss = loss( result.x )

print(result)
print('Final parameters:', result.x, 'Loss:',final_loss)
