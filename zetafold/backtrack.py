from .recursions.explicit_recursions import *
import random
import sys


##################################################################################################
def backtrack( self, backtrack_info_input, mode = 'mfe' ):
    '''
    modes are:
      mfe = backtrack, following maximum boltzmann weight. note that this is not *quite* MFE
      stochastic  = choose track based on boltzmann weights
      enumerative = follow all tracks!
    '''
    #print_backtrack_info( backtrack_info_input )
    if len( backtrack_info_input ) == 0: return []
    contrib_sum = sum( contrib[0] for contrib in backtrack_info_input )
    if   mode == 'enumerative':
        backtrack_info = [ contrib for contrib in backtrack_info_input ] # like a deepcopy
    elif mode == 'mfe':         backtrack_info = [ max_contrib(backtrack_info_input) ]
    elif mode == 'stochastic' : backtrack_info = [ get_random_contrib( backtrack_info_input ) ]
    p_bps = [] # list of tuples of (p_structure, bps_structure) for each structure
    N = self.N

    for contrib in backtrack_info: # each option ('contribution' to this partition function of this sub-region)
        if ( contrib[0] == 0.0 ): continue
        p_contrib = contrib[0]/contrib_sum
        p_bps_contrib = [ [p_contrib,[]] ]

        # each 'branch'; e.g., C_eff(i,k) Z_BP(k+1, j) has a C_eff and a Z_BP branch
        for backtrack_info in contrib[1]:
            ( Z_backtrack, i, j )  = backtrack_info
            if ( i == j ): continue
            for base_pair_type in self.params.base_pair_types:
                if Z_backtrack == self.Z_BPq[ base_pair_type ]:
                    base_pair = [i%N,j%N]
                    base_pair.sort()
                    # TODO: could also add type of base pair here -- we have the info!
                    p_bps_contrib = [ [p_bp[0], p_bp[1]+[tuple( base_pair )] ] for p_bp in p_bps_contrib ]
            next_backtrack_info = Z_backtrack.get_backtrack_info(self,i%N,j%N)
            p_bps_component = backtrack( self, next_backtrack_info, mode ) # recursion!
            if len( p_bps_component ) == 0: continue
            # put together all branches
            p_bps_contrib_new = []
            for p_bps1 in p_bps_contrib:
                for p_bps2 in p_bps_component:
                    p_bps_contrib_new.append( [p_bps1[0]*p_bps2[0], p_bps1[1]+p_bps2[1]] )
            p_bps_contrib = p_bps_contrib_new

        p_bps += p_bps_contrib
    return p_bps

##################################################################################################
def mfe( self, Z_final_contrib ):
    p_bps = backtrack( self, Z_final_contrib, mode = 'mfe' )
    assert( len(p_bps) == 1 )
    p,bps = p_bps[0]
    bps.sort()
    return (bps,p)

##################################################################################################
def boltzmann_sample( self, Z_final_contrib ):
    p_bps = backtrack( self, Z_final_contrib, mode = 'stochastic' )
    assert( len(p_bps) == 1 )
    return (p_bps[0][1],p_bps[0][0])

##################################################################################################
def enumerative_backtrack( self ):
    return backtrack( self, self.Z_final.get_backtrack_info(self,0), 'enumerative' )


##################################################################################################
def get_random_contrib( backtrack_info ):
    # Random sample weighted by probability. Must be a simple function for this.
    contrib_cumsum = [ backtrack_info[0][0] ]
    for contrib in backtrack_info[1:]: contrib_cumsum.append( contrib_cumsum[-1] + contrib[0] )
    r = random.random() * contrib_cumsum[ -1 ]
    for (idx,psum) in enumerate( contrib_cumsum ):
        if r < psum: return backtrack_info[idx]

def max_contrib(backtrack_info):
    max_contrib_val = None
    best_contrib = None
    for c in backtrack_info:
        if max_contrib_val is None or \
           c[0] > max_contrib_val:
            best_contrib = c
            max_contrib_val = c[0]
    return best_contrib

def print_contrib( contrib ):
    sys.stdout.write('[')
    print('%s:' % contrib[0])
    for n,backtrack_info in enumerate(contrib[1]):
        sys.stdout.write( '%s(%d,%d)' % (backtrack_info[0].name,backtrack_info[1],backtrack_info[2]) )
        if n < len( contrib[1] )-1: sys.stdout.write(',')
    sys.stdout.write(']')

def print_backtrack_info( backtrack_info ):
    print('[ ')
    for contrib in backtrack_info[:-1]:
        print_contrib( contrib )
        print('; ')
    print_contrib( backtrack_info[-1] )
    print('%s\n' % ' ]')
    return
