from __future__ import print_function
from .backtrack  import mfe, boltzmann_sample, enumerative_backtrack
from .parameters import get_params
from .util.wrapped_array  import WrappedArray, initialize_matrix
from .util.secstruct_util import *
from .util.output_util    import _show_results, _show_matrices
from .util.sequence_util  import initialize_sequence_and_ligated, initialize_all_ligated, get_num_strand_connections
from .util.constants import KT_IN_KCAL
from .util.assert_equal import assert_equal
from .derivatives import _get_log_derivs

from math import log

##################################################################################################
def partition( sequences, circle = False, params = '', mfe = False, calc_bpp = False,
               n_stochastic = 0, do_enumeration = False, structure = None, force_base_pairs = None, no_coax = False,
               verbose = False,  suppress_all_output = False,
               deriv_params = None,
               calc_Kd_deriv_DP = False, use_simple_recursions = False  ):
    '''
    Wrapper function into Partition() class
    Returns Partition object p which holds results like:

      p.Z   = final partition function (where unfolded state has unity weight)
      p.bpp = matrix of base pair probabilities (if requested by user with calc_bpp = True)
      p.struct_MFE = minimum free energy secondary structure in dot-parens notation
      p.bps_MFE  = minimum free energy secondary structure as sorted list of base pairs
      p.dZ_dKd_DP = derivative of Z w.r.t. Kd computed in-line with dynamic programming (if requested by user with calc_Kd_deriv_DP = True)

    '''
    if isinstance(params,str): params = get_params( params, suppress_all_output )
    if no_coax:                params.K_coax = 0.0

    p = Partition( sequences, params )
    p.calc_all_elements = calc_bpp or (deriv_params != None)
    p.use_simple_recursions = use_simple_recursions
    p.circle    = circle
    p.options.calc_deriv_DP = calc_Kd_deriv_DP
    p.structure = get_structure_string( structure )
    p.force_base_pairs = get_structure_string( force_base_pairs )
    p.suppress_all_output = suppress_all_output
    p.deriv_params = deriv_params
    p.run()
    if calc_bpp:         p.get_bpp_matrix()
    if mfe:              p.calc_mfe()
    if n_stochastic > 0: p.stochastic_backtrack( n_stochastic )
    if do_enumeration:   p.enumerative_backtrack()
    if verbose:          p.show_matrices()
    if not suppress_all_output: p.show_results()
    p.run_cross_checks()

    return p

##################################################################################################
class Partition:
    '''
    Statistical mechanical model for RNA folding, testing a bunch of extensions and with lots of cross-checks.
    (C) R. Das, Stanford University, 2018
    '''
    def __init__( self, sequences, params ):
        '''
        Required user input.
        sequences = string with sequence, or array of strings (sequences of interacting strands)
        params    = AlphaFoldParams object
        '''
        self.sequences = sequences
        self.params = params
        self.circle = False  # user can update later --> circularize sequence
        self.use_simple_recursions = False
        self.calc_all_elements     = False
        self.calc_bpp = False
        self.base_pair_types = params.base_pair_types
        self.suppress_all_output = False
        self.structure = None
        self.force_base_pairs = None
        self.deriv_params = None
        self.options = PartitionOptions()

        # for output:
        self.Z       = 0
        self.dG      = None
        self.bpp     = []
        self.bps_MFE = []
        self.struct_MFE = ''
        self.struct_stochastic = []
        self.struct_enumerate  = []
        self.log_derivs = []
        self.derivs     = []
        return

    ##############################################################################################
    def run( self ):
        '''
        Do the dynamic programming to fill partition function matrices
        '''
        initialize_sequence_information( self ) # N, sequence, ligated, all_ligated
        initialize_dynamic_programming_matrices( self ) # ( Z_BP, C_eff, Z_linear, Z_cut, Z_coax, etc. )
        initialize_force_base_pair( self )

        # do the dynamic programming
        for offset in range( 1, self.N ): #length of subfragment
            for i in range( self.N ):     #index of subfragment
                if (not self.calc_all_elements) and ( i + offset ) >= self.N: continue
                j = (i + offset) % self.N;  # N cyclizes
                for Z in self.Z_all: Z.update( self, i, j )

        for i in range( self.N): self.Z_final.update( self, i )

        self.log_derivs = self.get_log_derivs( self.deriv_params )

        fill_in_outputs( self )

    # boring member functions -- defined later.
    def get_bpp_matrix( self ): _get_bpp_matrix( self ) # fill base pair probability matrix
    def calc_mfe( self ): _calc_mfe( self )
    def stochastic_backtrack( self, N ): _stochastic_backtrack( self, N )
    def enumerative_backtrack( self ): _enumerative_backtrack( self )
    def show_results( self ): _show_results( self )
    def show_matrices( self ): _show_matrices( self )
    def get_log_derivs( self, deriv_params ): return _get_log_derivs( self, deriv_params )
    def run_cross_checks( self ): _run_cross_checks( self )
    def num_strand_connections( self ):  return get_num_strand_connections( self.sequences, self.circle)

##################################################################################################
def fill_in_outputs( self ):
    self.Z  = self.Z_final.val(0)
    if self.Z > 0.0: self.dG = -KT_IN_KCAL * log( self.Z )
    self.dZ_dKd_DP = self.Z_final.deriv(0)
    self.derivs = []
    if self.deriv_params:
        for n,log_deriv in enumerate(self.log_derivs):
            self.derivs.append( log_deriv * self.Z / self.params.get_parameter_value( self.deriv_params[n] )  )

def initialize_sequence_information( self ):
    '''
    Create sequence information from sequences of strands:

    INPUT:
    sequences = sequences of interacting strands (array of strings)
    circle    = user asks for nucleotides N and 1 to be ligated ('circularized') (bool)

    OUTPUT:
    sequence     = concatenated sequence (string, length N)
    is_ligated   = is not a cut ('nick','chainbreak') (Array of bool, length N)
    all_ligated  = no cutpoint exists between i and j (N X N)
    '''
    # initialize sequence
    self.sequence, self.ligated, self.sequences = initialize_sequence_and_ligated( self.sequences, self.circle, use_wrapped_array = self.use_simple_recursions )
    self.N = len( self.sequence )
    self.all_ligated = initialize_all_ligated( self.ligated )

##################################################################################################
class PartitionOptions:
    def __init__( self ):
        self.calc_deriv_DP = False
        self.calc_contrib  = False

##################################################################################################
def initialize_dynamic_programming_matrices( self ):
    '''
    A bunch of zero matrices. Only non-trivial thing is
    initialization of (i,i) [diagonal]:
         Z_BP(i,i)     = 0
         C_eff(i,i)    = C_init (units of M)
         Z_linear(i,i) = 1
    And: the order of initialization in the list Z_all will
      determine the actual order of updates during dynamic programmming at each (i,j).
    '''

    from .recursions.explicit_recursions import update_Z_BPq, update_Z_BP, update_Z_cut, update_Z_coax, update_C_eff_basic, update_C_eff_no_BP_singlet, update_C_eff_no_coax_singlet, update_C_eff, update_Z_final, update_Z_linear
    from .recursions.explicit_dynamic_programming import DynamicProgrammingMatrix, DynamicProgrammingList
    if self.use_simple_recursions: # over-ride with simpler recursions that are easier for user to input.
        from .recursions.recursions import update_Z_BPq, update_Z_BP, update_Z_cut, update_Z_coax, update_C_eff_basic, update_C_eff_no_BP_singlet, update_C_eff_no_coax_singlet, update_C_eff, update_Z_final, update_Z_linear
        from .recursions.dynamic_programming import DynamicProgrammingMatrix, DynamicProgrammingList

    N = self.N

    # Collection of all N X N dynamic programming matrices -- order in this list will
    #  determine order of updates.
    self.Z_all = Z_all = []

    # some preliminary helpers
    self.Z_cut    = DynamicProgrammingMatrix( N, DPlist = Z_all, update_func = update_Z_cut, options = self.options );

    # base pairs and co-axial stacks
    self.Z_BPq = {}
    for base_pair_type in self.base_pair_types:
        # the bpt = base_pair_type holds the base_pair_type info in the lambda (Python FAQ)
        update_func = lambda partition,i,j,bpt=base_pair_type: update_Z_BPq(partition,i,j,bpt)
        self.Z_BPq[ base_pair_type ] = DynamicProgrammingMatrix( N, DPlist = Z_all,
                                                                 update_func = update_func, options = self.options )
    self.Z_BP     = DynamicProgrammingMatrix( N, DPlist = Z_all, update_func = update_Z_BP, options = self.options );
    self.Z_coax   = DynamicProgrammingMatrix( N, DPlist = Z_all, update_func = update_Z_coax, options = self.options );

    # C_eff makes use of information on Z_BP, so compute last
    C_init = self.params.C_init
    self.C_eff_basic           = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff_basic, options = self.options );
    self.C_eff_no_BP_singlet   = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff_no_BP_singlet, options = self.options );
    self.C_eff_no_coax_singlet = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff_no_coax_singlet, options = self.options );
    self.C_eff                 = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff, options = self.options );

    self.Z_linear = DynamicProgrammingMatrix( N, diag_val = 1.0, DPlist = Z_all, update_func = update_Z_linear, options = self.options );

    # Last DP 1-D list (not a 2-D N x N matrix)
    self.Z_final = DynamicProgrammingList( N, update_func = update_Z_final, options = self.options  )

    self.params.check_C_eff_stack()

##################################################################################################
def initialize_force_base_pair( self ):
    self.allow_base_pair     = None
    self.in_forced_base_pair = None
    if self.structure != None:
        assert( self.force_base_pairs == None )
        bp_list = bps( self.structure )
    elif self.force_base_pairs != None:
        assert( self.structure == None )
        bp_list = bps( self.force_base_pairs )
    else:
        return

    N = self.N
    self.in_forced_base_pair = [False] * N
    if self.use_simple_recursions: self.in_forced_base_pair = WrappedArray( N, False )

    for i,j in bp_list:
        self.in_forced_base_pair[ i ] = True
        self.in_forced_base_pair[ j ] = True
        self.Z_linear.set_val( i, i, 0.0 )
        self.Z_linear.set_val( j, j, 0.0 )
        self.C_eff.set_val( i, i, 0.0 )
        self.C_eff.set_val( j, j, 0.0 )

    if self.structure != None:
        self.allow_base_pair = initialize_matrix( N, False )
        for i,j in bp_list:
            self.allow_base_pair[ i ][ j ] = True
            self.allow_base_pair[ j ][ i ] = True
    else:
        assert( self.force_base_pairs != None )
        # allow the specified base pairs (indeed, force them), but
        #  now disallow any base pairs that might cross with them.
        self.allow_base_pair = initialize_matrix( N, True )
        for i,j in bp_list:
            # no crossing pairs
            for m in range( i+1, j ):
                for n in range( j+1, i+N ):
                    if (( m - n ) % N) == 0: continue
                    self.allow_base_pair[ m ][ n ] = False
                    self.allow_base_pair[ n ][ m ] = False
            # no other partners
            for m in range( N ):
                if m != j:
                    self.allow_base_pair[ i ][ m ] = False
                    self.allow_base_pair[ m ][ i ] = False
                if m != i:
                    self.allow_base_pair[ j ][ m ] = False
                    self.allow_base_pair[ m ][ j ] = False

##################################################################################################
def _get_bpp_matrix( self ):
    '''
    Getting base pair probability matrix.
    Gets carried out pretty fast since we've already computed the sum over structures in i..j encapsulated by a pair (i,j), as well
      as structures in j..i encapsulated by those pairs.
    So: it becomes easy to calculate partition function over all structures with base pair (i,j), and then divide by total Z.
    '''
    assert( self.calc_all_elements )
    self.bpp = initialize_matrix( self.N, 0.0 )
    for i in range( self.N ):
        for j in range( self.N ):
            self.bpp[i][j] = 0.0
            for base_pair_type in self.params.base_pair_types:
                self.bpp[i][j] += self.Z_BPq[base_pair_type].val(i,j) * self.Z_BPq[base_pair_type.flipped].val(j,i) * base_pair_type.Kd / self.Z_final.val(0)

##################################################################################################
def _calc_mfe( self ):
    #
    # Wrapper into mfe(), written out in backtrack.py
    #
    N = self.N
    p_MFE   = [0.0]*N
    bps_MFE = [[]]*N

    # there are actually numerous ways to calculate MFE if we did all N^2 elements -- let's check.
    n_test = N if self.calc_all_elements else 1
    if not self.suppress_all_output:
        print()
        print('Doing backtrack to get minimum free energy structure:')
        print(self.sequence)

    for i in range( n_test ):
        (bps_MFE[i], p_MFE[i] ) = mfe( self, self.Z_final.get_contribs(self,i) )
        assert_equal( p_MFE[i], p_MFE[0] )
        # actually this doesn't always hold -- in some parameter sets and sequences there are literally ties.
        # assert( bps_MFE[i] == bps_MFE[0] )

    if not self.suppress_all_output:
        print( secstruct(bps_MFE[0],N), "   ", p_MFE[0], "[MFE]")
        print()
    self.bps_MFE = bps_MFE[0]
    self.struct_MFE = secstruct( bps_MFE[0], N)

##################################################################################################
def _stochastic_backtrack( self, N_backtrack ):
    #
    # Get stochastic, Boltzmann-weighted structural samples from partition function
    #
    print()
    print('Doing',N_backtrack,'stochastic backtracks to get Boltzmann-weighted ensemble')
    print(self.sequence)
    for i in range( N_backtrack ):
        bps, p = boltzmann_sample( self, self.Z_final.get_contribs(self,0) )
        print(secstruct(bps,self.N), "   ", p, "[stochastic]")
        self.struct_stochastic.append( secstruct(bps,self.N) )
    print()

    return

##################################################################################################
def _enumerative_backtrack( self ):
    #
    # Enumerate all structures, and track their probabilities
    #
    print()
    print('Doing complete enumeration of Boltzmann-weighted ensemble...')
    print(self.sequence)
    p_bps = enumerative_backtrack( self )
    for (p,bps) in p_bps:
        print(secstruct(bps,self.N), "   ", p, "[enumerative]")
        self.struct_enumerate.append( secstruct(bps,self.N) )
    p_tot = sum( p_bp[0] for p_bp in p_bps )
    print('p_tot = ',p_tot)
    assert( abs(p_tot - 1.0) < 1.0e-5 )
    return

##################################################################################################
def _run_cross_checks( self ):
    # stringent test that partition function is correct -- all the Z(i,i) agree.
    if self.calc_all_elements:
        for i in range( self.N ): assert_equal( self.Z_final.val(0), self.Z_final.val(i) )

        if self.options.calc_deriv_DP and self.Z_final.deriv(0) > 0:
            for i in range( self.N ): assert_equal( self.Z_final.deriv(0), self.Z_final.deriv(i) )

    # calculate bpp_tot = -dlog Z_final /dlog Kd in up to three ways! wow cool test
    if len(self.bpp)>0:
        bpp_tot = 0.0
        for i in range( self.N ):
            for j in range( self.N ):
                bpp_tot += self.bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)

        # uh this is a hack -- only works for minimal model where all the Kd are the same:
        Kd = self.params.base_pair_types[0].Kd
        if self.options.calc_deriv_DP:
            bpp_tot_based_on_deriv = -self.Z_final.deriv(0) * Kd / self.Z_final.val(0)
            print('bpp_tot',bpp_tot,'bpp_tot_based_on_deriv',bpp_tot_based_on_deriv)
            if bpp_tot > 0: assert_equal( bpp_tot, bpp_tot_based_on_deriv )




