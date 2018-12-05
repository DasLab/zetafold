#!/usr/bin/python
from zetafold.parameters import get_params
from zetafold.data.training_examples import *
from zetafold.data.config import *
from zetafold.training import *
from scipy.optimize import minimize
from multiprocessing import Pool

# set up starter parameters
params = get_params( suppress_all_output = True )
params.set_parameter( 'K_coax', 0.0 )

# set up training examples
training_examples = [ tRNA ] #, P4P6_outerjunction ]
train_parameters = ['C_init','C_eff_stack_WC_WC','C_eff_stack_WC_GU','C_eff_stack_WC_UG','C_eff_stack_GU_GU','C_eff_stack_UG_GU','C_eff_stack_GU_UG']
x0      = np.array( [     0.5,                5.3,                4.1,                5.6,                4.9,                4,                    4] ) * np.log(10.0)

train_parameters = ['C_init']
x0      = np.array( [     0.5 ] )

for training_example in training_examples:
    training_example.params = params
    training_example.train_parameters = train_parameters

pool = Pool(4)

loss = lambda x:free_energy_gap(      x,params,train_parameters,training_examples,pool)
grad = lambda x:free_energy_gap_deriv(x,params,train_parameters,training_examples,pool)

result = minimize( loss, x0, method = 'L-BFGS-B') #, jac = free_energy_gap_deriv )
final_loss = loss( result.x )
print( 'Deriv: ', grad( result.x ) )

print(result)
print('Final parameters:', result.x, 'Loss:',final_loss)
