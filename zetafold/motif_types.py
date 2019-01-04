from zetafold.base_pair_types import get_base_pair_types_for_tag

class MotifType:
    def __init__( self, strands, bp_tags, C_eff, params ):
        '''
        MotifType stores information on specially stable motifs, including eventually
         UA-handles, symmetric internal loops

        Currently only handling 2-way junctions.

        Note that base pair stacks (C_eff_stack) could also be handled by this MotifType object, but
          it turns out that the get_match_base_pair_type_sets() function below is just too damn slow.

        For speed, its better to handle base-pair-step 'motifs' separately.

        Example 1:
        ----------
             a
          5'-CAC-3'
          d  : : b
          3'-G-G-5'
             c

        has a, b, c, d going in a circle. These are recorded in order.

          strands = ['CAC', 'GG']   (a, c)
          bp_tags = [None,  None]   (b, d)

        and internally this will become:

          base_pair_type_sets = [ [BasePairType('CG')], [BasePairType('GC')]]   (b, d, go in a 'circle')

        Perhaps a little confusing because the last base pair type set has d (the 'start base pair', flipped)


        Example 2:
        ----------
              a
          5'-NAG-3'
    WC--> d  : : b
          3'-N-U-5'
              c

          strands = ['NAG', 'UN']   (a, c)
          bp_tags = [None, 'WC']   (b, d)
          base_pair_type_sets = [ [BasePairType('GU')], [BasePairType('CG'),BasePairType('GC'),BasePairType('AU'),BasePairType('UA')] ]
                                       since WC is interpreted as (strict) Watson-Crick pair

        TODO: handle hairpins, and document example
        TODO: handle 3WJ, 4WJ, and document example
        '''
        self.strands = strands
        self.bp_tags = bp_tags
        M = len( self.strands )
        assert( len(bp_tags) == M)

        self.base_pair_type_sets = []
        for i,bp_tag in enumerate(bp_tags):
            if bp_tag == None: bp_tag = strands[i][-1] + strands[ (i+1)%M ][0]
            self.base_pair_type_sets.append( get_base_pair_types_for_tag( params, bp_tag ) )

        self.C_eff = C_eff
        self.permuted = None # for now.

    def get_match_base_pair_type_sets( self, sequence, ligated, i, j ):
        '''
        TODO: generalize to N-way junction
        '''
        N = len( sequence )
        M = len( self.strands )
        match_base_pair_type_sets = []

        # works for hairpins and internal loops both:
        for q in range( len( self.strands[0] ) ):
            if self.strands[0][q] != 'N' and sequence[(i+q)%N] != self.strands[0][q]: return None
        for q in range( len( self.strands[0] )-1 ):
            if not ligated[(i+q)%N]: return None

        if M == 1:
            # uh this is kind of explicit. must be a more general way to treat M = 1, 2, ...
            if not (( j - i) % N) == len( self.strands[0] )-1: return None

        if M == 2:
            # uh this is kind of explicit. must be a more general way to treat M = 1, 2, ...
            if (( j - i) % N) < len( self.strands[0] ) + len( self.strands[1] ) - 1: return None

            matches = []
            for base_pair_type in self.base_pair_type_sets[ 0 ]:
                i_next, j_next = (i+len(self.strands[0])-1)%N, (j-len(self.strands[-1])+1)%N
                if base_pair_type.is_match( sequence[i_next],sequence[j_next] ): matches.append( (base_pair_type,i_next,j_next) )
            if len( matches ) == 0: return None
            match_base_pair_type_sets.append( matches )

                # for the second strand of the internal loop
            for q in range( len( self.strands[-1] ) ):
                if self.strands[-1][q] != 'N' and sequence[(j-len(self.strands[-1])+1+q)%N] != self.strands[-1][q]: return None
            for q in range( len( self.strands[-1] )-1 ):
                if not ligated[ (j-len(self.strands[-1])+1+q)%N ]: return None

        if M > 2:
            print('%s\n' % 'Cannot handle 3WJ or 4WJ motifs yet!')
            exit()

        # redundant with start_base_pair_type -- using as cross check
        matches = []
        for base_pair_type in self.base_pair_type_sets[ -1 ]:
            if base_pair_type.is_match( sequence[j%N], sequence[i%N] ): matches.append( (base_pair_type,j,i) )
        if len( matches ) == 0: return None
        match_base_pair_type_sets.append( matches )

        return match_base_pair_type_sets

    def get_tag( self ): return make_motif_type_tag( self.strands, self.bp_tags )

def make_motif_type_tag( strands, bp_tags ):
    tag = ''
    for i in range( len( strands ) ):
        if i > 0: tag += '_'
        tag += strands[i]
        if bp_tags[i]: tag += '_'+bp_tags[i]
    return tag


def get_motif_type_for_tag( params, tag ):
    for motif_type in params.motif_types:
        if motif_type.get_tag() == tag:
            return motif_type
    return None

def parse_motif_type_tag( motif_type_tag ):
    '''

        Example
        ----------

        NAG_NU_WC

              a
          5'-NAG-3'
    WC--> d  : : b
          3'-N-U-5'
              c

       becomes:

          strands = ['NAG', 'UN']   (a, c)
          bp_tags = [None, 'WC']   (b, d)

    '''
    strands = []
    bp_tags = []
    tags = motif_type_tag.split( '_' )
    for tag in tags:
        if tag == 'WC':
            bp_tags[-1] = 'WC'
            continue
        strands.append( tag )
        bp_tags.append( None )
    return (strands, bp_tags)

def check_equivalent_C_eff_stack_for_motif_type( strands, bp_tags = None):
    '''
    Will convert CC_GG    to C_eff_stack_CG_CG  (bp1_bp2)
    Will convert CN_WC_NG to C_eff_stack_CG_WC  (bp1_bp2)
    '''
    if isinstance( strands, str ) and bp_tags == None:
        (strands,bp_tags) = parse_motif_type_tag( strands )

    if len(strands) == 2 and len( strands[0] ) == 2 and len( strands[1] ) == 2:
        tag1 = bp_tags[1] if bp_tags[1] != None else strands[0][0]+strands[1][1]
        tag2 = bp_tags[0] if bp_tags[0] != None else strands[0][1]+strands[1][0]
        return 'C_eff_stack_%s_%s' % (tag1, tag2)
    return None
