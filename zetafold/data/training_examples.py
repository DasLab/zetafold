import os
from glob import glob

class TrainingExample:
    '''
    Info needed for train_zetafold
    '''
    def __init__( self, name, sequence, structure, force_base_pairs = None ):
        self.name = name
        self.sequence = sequence
        self.structure = structure
        self.force_base_pairs = force_base_pairs

all_training_examples = {}
training_sets = {}
training_set_names = []

def read_in_training_examples( training_example_file ):
    fid = open( training_example_file, 'r' )
    line = fid.readline()
    while line:
        name = line[:-1]
        line = fid.readline()
        sequence  = line[:-1]
        line = fid.readline()
        structure = line[:-1]
        line = fid.readline()
        if line and len(line[:-1].replace(' ','')) > 0:
            force_base_pairs = line[:-1]
            line = fid.readline()
        else:
            force_base_pairs = None
        all_training_examples[ name ] = TrainingExample( name, sequence, structure, force_base_pairs )
        line = fid.readline()
    fid.close()

def read_in_training_sets( training_set_file ):
    fid = open( training_set_file, 'r' )
    line = fid.readline()
    while line:
        training_set_name = line[:-1]
        line = fid.readline()
        training_examples  = line[:-1].split(' ')
        line = fid.readline()
        training_set_names.append( training_set_name )
        training_sets[ training_set_name ] = training_examples
        line = fid.readline()
    fid.close()

data_dir = os.path.dirname( os.path.abspath(__file__) )

training_example_files = glob( data_dir + '/*.examples.txt' )
training_example_files.sort()
for training_example_file in training_example_files: read_in_training_examples( training_example_file )

training_set_files = glob( data_dir + '/*.sets.txt' )
training_set_files.sort()
for training_set_file in training_set_files: read_in_training_sets( training_set_file )
