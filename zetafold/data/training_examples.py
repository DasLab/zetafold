class TrainingExample:
    '''
    Info needed for train_zetafold
    '''
    def __init__( self, name, sequence, structure, force_base_pairs = None ):
        self.name = name
        self.sequence = sequence
        self.structure = structure
        self.force_base_pairs = force_base_pairs

tRNA = TrainingExample( 'tRNA',
                        'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA',
                        '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....' )
P5abc = TrainingExample( 'P5abc',
                         'CCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGG',
                         '(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))' )
P4P6_outerjunction = TrainingExample( 'P4P6_outerjunction',
                                      'GGAAUUGCGGGAAAGG CUAACCACGCAGCCAAGUCCUAAGUC GAUAUGGAUGCAGUUCA',
                                      '.....((((((...(( ))..)).))))((...((((...((( )))..))))...))...',
                                      '...............( )........................( )................' )
add = TrainingExample( 'add',
                        'CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG',
                        '(((((((((...((((((.........))))))........((((((.......))))))..)))))))))' )
# P4-P6
P4P6 = TrainingExample( 'P4P6',
                        'GGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA',
                        '.....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...))...' )


training_sets = {}
training_sets[ 'tRNA' ]    = [tRNA]
training_sets[ 'minimal' ] = [tRNA,P5abc,P4P6_outerjunction,add]
training_sets[ 'full' ]    = [tRNA,P4P6,add]

