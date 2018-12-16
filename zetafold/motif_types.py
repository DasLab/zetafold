class MotifType:
    def __init__( self, start_base_pair_type, strands, base_pair_types, C_eff ):
        '''
        MotifType stores information on specially stable motifs, including
         UA-handles, symmetric internal loops
        Currently only handling 2-way junctions
        '''
        self.start_base_pair_type = start_base_pair_type
        self.strands = strands
        self.base_pair_types = base_pair_types
        self.C_eff = C_eff
        self.permuted = None # for now.

    def is_match( self, sequence, ligated, i, j ):
        N = len( sequence )

        if not self.start_base_pair_type.is_match( sequence[i%N], sequence[j%N] ): return False

        for q in range( len( self.strands[0] ) ):
            if sequence[(i+q)%N] != self.strands[0][q]: return False
        for q in range( len( self.strands[0] )-1 ):
            if not ligated[(i+q)%N]: return False

        if not self.base_pair_types[0].is_match( sequence[(i+len(self.strands[ 0])-1)%N],
                                                 sequence[(j-len(self.strands[-1])+1)%N] ): return False

        for q in range( len( self.strands[-1] ) ):
            if sequence[(j-len(self.strands[-1])+1+q)%N] != self.strands[-1][q]: return False
        for q in range( len( self.strands[-1] )-1 ):
            if not ligated[ (j-len(self.strands[-1])+1+q)%N ]: return False

        # redundant with start_base_pair_type -- using as cross check
        if not self.base_pair_types[-1].is_match( sequence[j%N], sequence[i%N] ): return False

        return True

    def get_other_base_pair( self, i, j ):
        return ( self.base_pair_types[0], i+len(self.strands[0])-1, j-len(self.strands[-1])+1 )


    def get_tag( self ):
        return 'startbp'+self.start_base_pair_type.get_tag() + \
            '_strand'+self.strands[0]+'_bp'+self.base_pair_types[0].get_tag() + \
            '_strand'+self.strands[1]+'_bp'+self.base_pair_types[1].get_tag()

def get_motif_type_for_tag( params, tag ):
    for motif_type in params.motif_types:
        if motif_type.get_tag() == tag:
            return motif_type
    return None
