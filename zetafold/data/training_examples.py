class TrainingExample:
    '''
    Info needed for train_zetafold
    '''
    def __init__( self, sequence, structure, force_base_pairs = None ):
        self.sequence = sequence
        self.structure = structure
        self.force_base_pairs = force_base_pairs

tRNA = TrainingExample( 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA',
                        '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....' )

# P4-P6
P4P6_sequence =  'GGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA'
P4P6_structure = '.....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...))...'

# P5abc
P5abc_sequence =  'CCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGG'
P5abc_structure = '(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))'

# P4-P6 outer junction -- further minimized -- easier to read input
P4P6_outerjunction_sequence  = 'GGAAUUGCGGGAAAGG CUAACCACGCAGCCAAGUCCUAAGUC GAUAUGGAUGCAGUUCA'
P4P6_outerjunction_structure = '.....((((((...(( ))..)).))))((...((((...((( )))..))))...))...'
P4P6_outerjunction_force_bps = '...............( )........................( )................'

add_sequence  = 'CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG'
add_structure = '(((((((((...((((((.........))))))........((((((.......))))))..)))))))))'
