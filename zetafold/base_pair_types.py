class BasePairType:
    def __init__( self, nt1, nt2, Kd, match_lowercase = False ):
        '''
        Two sequence characters that get matched to base pair, e.g., 'C' and 'G';
           and the dissociation constant K_d (in M) associated with initial pair formation.
        match_lowercase means match 'x' to 'x', 'y' to 'y', etc.
        TODO: also store chemical modification info.
        '''
        self.nt1 = nt1
        self.nt2 = nt2
        self.Kd = Kd
        self.match_lowercase = ( nt1 == '*' and nt2 == '*' and match_lowercase )
        self.flipped = self # needs up be updated later.

    def is_match( self, s1, s2 ):
        if self.match_lowercase: return ( s1.islower() and s2.islower() and s1 == s2 )
        return ( s1 == self.nt1 and s2 == self.nt2 )

    def get_tag( self ):
        if self.match_lowercase: return 'matchlowercase'
        return self.nt1+self.nt2

def setup_base_pair_type( params, nt1, nt2, Kd, match_lowercase = False ):
    if not hasattr( params, 'base_pair_types' ): params.base_pair_types = []
    bpt1 = BasePairType( nt1, nt2, Kd, match_lowercase = match_lowercase )
    params.base_pair_types.append( bpt1 )
    if not match_lowercase:
        bpt2 = BasePairType( nt2, nt1, Kd, match_lowercase = match_lowercase )
        bpt1.flipped = bpt2
        bpt2.flipped = bpt1
        params.base_pair_types.append( bpt2 )

def get_base_pair_type_for_tag( params, tag ):
    if not hasattr( params, 'base_pair_types' ): return None
    for base_pair_type in params.base_pair_types:
        if (tag == 'matchlowercase' and base_pair_type.match_lowercase) or \
           (tag == base_pair_type.nt1 + base_pair_type.nt2 ):
            return base_pair_type
    #print( 'Could not figure out base_pair_type for ', tag )
    return None

def get_base_pair_types_for_tag( params, tag ):
    if not hasattr( params, 'base_pair_types' ): return None
    if tag == 'WC':
        WC_nts =  [('C','G'),('G','C'),('A','U'),('U','A'),('G','U'),('U','G')]
        base_pair_types = []
        for base_pair_type in params.base_pair_types:
            if (base_pair_type.nt1,base_pair_type.nt2) in WC_nts:
                base_pair_types.append( base_pair_type )
        return base_pair_types
    else:
        return [ get_base_pair_type_for_tag( params, tag ) ]
    return None

def initialize_base_pair_types( self ):
    self.base_pair_types = []

