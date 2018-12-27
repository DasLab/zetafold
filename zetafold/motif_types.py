from base_pair_types import get_base_pair_types_for_tag

class MotifType:
    def __init__( self, strands, bp_tags, C_eff, params ):
        '''
        MotifType stores information on specially stable motifs, including eventually
         UA-handles, symmetric internal loops
        Currently only handling 2-way junctions.

        Example 1:
        ----------
             a
          5'-CAC-3'
          d  : : b
          3'-G-G-5'
             c

        has a, b, c, d going in a circle. These are recorded in order.

          strands = ['CAC', 'GG']   (a, c)
          bp_tags = ['CG', 'GC']   (b, d)

        and internally this will become:

          base_pair_type_sets = [ [BasePairType('CG')], [BasePairType('GC')]]   (b, d, go in a 'circle')

        Perhaps a little confusing because the last base pair type set has d (the 'start base pair', flipped)


        Example 2:
        ----------
             a
          5'-NG-3'
    WC--> d  :: b
          3'-NU-5'
             c

          strands = ['NG', 'UN']   (a, c)
          base_pair_type_sets = [ [BasePairType('GU')], [BasePairType('CG'),BasePairType('GC'),BasePairType('AU'),BasePairType('UA')] ]
                                       since WC is interpreted as (strict) Watson-Crick pair

        TODO: handle hairpins, and document example
        TODO: handle 3WJ, 4WJ, and document example
        '''
        self.strands = strands
        self.bp_tags = bp_tags
        N = len( self.strands )
        assert( len(bp_tags) == N)

        self.base_pair_type_sets = []
        for i,bp_tag in enumerate(bp_tags):
            if bp_tag == None: bp_tag = strands[i][-1] + strands[ (i+1)%N ][0]
            self.base_pair_type_sets.append( get_base_pair_types_for_tag( params, bp_tag ) )

        self.C_eff = C_eff
        self.permuted = None # for now.

    def get_match_base_pair_type_sets( self, sequence, ligated, i, j ):
        '''
        TODO: generalize to N-way junction
        '''
        N = len( sequence )
        match_base_pair_type_sets = []

        for q in range( len( self.strands[0] ) ):
            if self.strands[0][q] != 'N' and sequence[(i+q)%N] != self.strands[0][q]: return None
        for q in range( len( self.strands[0] )-1 ):
            if not ligated[(i+q)%N]: return None

        matches = []
        for base_pair_type in self.base_pair_type_sets[ 0 ]:
            i_next, j_next = (i+len(self.strands[0])-1)%N, (j-len(self.strands[-1])+1)%N
            if base_pair_type.is_match( sequence[i_next],sequence[j_next] ): matches.append( (base_pair_type,i_next,j_next) )
        if len( matches ) == 0: return None
        match_base_pair_type_sets.append( matches )

        for q in range( len( self.strands[-1] ) ):
            if self.strands[-1][q] != 'N' and sequence[(j-len(self.strands[-1])+1+q)%N] != self.strands[-1][q]: return None
        for q in range( len( self.strands[-1] )-1 ):
            if not ligated[ (j-len(self.strands[-1])+1+q)%N ]: return None

        # redundant with start_base_pair_type -- using as cross check
        matches = []
        for base_pair_type in self.base_pair_type_sets[ -1 ]:
            if base_pair_type.is_match( sequence[j%N], sequence[i%N] ): matches.append( (base_pair_type,j,i) )
        if len( matches ) == 0: return None
        match_base_pair_type_sets.append( matches )

        return match_base_pair_type_sets

    def get_tag( self ):
        tag = ''
        for i in range( len( self.strands ) ):
            if i > 0: tag += '_'
            tag += self.strands[i]
            if self.bp_tags[i]: tag += '_'+self.bp_tags[i]
        return tag

def get_motif_type_for_tag( params, tag ):
    for motif_type in params.motif_types:
        if motif_type.get_tag() == tag:
            return motif_type
    return None


