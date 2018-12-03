from __future__ import print_function
import math
from .base_pair_types import BasePairType, setup_base_pair_type, get_base_pair_types_for_tag, get_base_pair_type_for_tag
from .util.constants import KT_IN_KCAL
import os.path

class AlphaFoldParams:
    '''
    Parameters that define the statistical mechanical model for RNA folding
    '''
    def __init__( self ):
        self.C_std  = 1.0      # 1 M. drops out in end (up to overall scale factor).
        self.allow_strained_3WJ = False

    def get_variables( self ):
        if self.C_init == 0.0 and self.name == 'empty': print('WARNING! C_init not defined, and params appear empty. Look at get_minimal_params() or get_latest_params() for examples')
        return ( self.C_init, self.l, self.l_BP, self.K_coax, self.l_coax, self.C_std, self.min_loop_length, self.allow_strained_3WJ )

    def check_C_eff_stack( self ): _check_C_eff_stack( self )

def get_params( params = None, suppress_all_output = False ):
    params_object = None
    if isinstance(params,AlphaFoldParams): return params
    elif params == None or params =='': params_object = get_latest_params()
    elif params == 'minimal': params_object = get_params_from_file( 'minimal' )
    elif params == 'v0.1':    params_object = get_params_from_file( 'zetafold_v0.1' )
    elif params == 'v0.15':   params_object = get_params_from_file( 'zetafold_v0.15' )
    elif params == 'v0.16':   params_object = get_params_from_file( 'zetafold_v0.16' )
    elif params == 'v0.17':   params_object = get_params_from_file( 'zetafold_v0.17' )
    elif params == 'v0.171':  params_object = get_params_from_file( 'zetafold_v0.171' )
    else: print('unrecognized params requested: ', params)
    if not suppress_all_output: print('Parameters: ', params_object.name, ' version', params_object.version)
    return params_object

def get_latest_params():
    # TODO check all params files and pick latest version
    return get_params_from_file( 'zetafold_v0.171' )

def update_C_eff_stack( params, val = None ):
    if not hasattr( params, 'C_eff_stack' ): params.C_eff_stack = {}
    for bpt1 in params.base_pair_types:
        if not params.C_eff_stack.has_key( bpt1 ):  params.C_eff_stack[ bpt1 ] = {}
        for bpt2 in params.base_pair_types:
            params.C_eff_stack[ bpt1 ][ bpt2 ] = val

def _check_C_eff_stack( params ):
    for bpt1 in params.base_pair_types:
        for bpt2 in params.base_pair_types:
            if ( params.C_eff_stack[ bpt1 ][ bpt2 ] != params.C_eff_stack[ bpt2.flipped ][ bpt1.flipped ] ):
                print("PROBLEM with C_eff_stacked pair!!!", bpt1.nt1, bpt1.nt2, " to ", bpt2.nt1, bpt2.nt2, params.C_eff_stack[ bpt1 ][ bpt2 ],
                    ' does not match ' , \
                    bpt2.flipped.nt1, bpt2.flipped.nt2, " to ", bpt1.flipped.nt1, bpt1.flipped.nt2, params.C_eff_stack[ bpt2.flipped ][ bpt1.flipped ] )
            assert( params.C_eff_stack[ bpt1 ][ bpt2 ] == params.C_eff_stack[ bpt2.flipped ][ bpt1.flipped ] )

def setup_base_pair_type_by_tag( params, Kd_tag, val ):
    tag = Kd_tag[3:]
    if get_base_pair_type_for_tag( params, tag ) != None: return
    if tag == 'matchlowercase':
        setup_base_pair_type( params, '*', '*', val, match_lowercase = True )
    else:
        setup_base_pair_type( params, tag[0], tag[1], val, match_lowercase = False )

def set_parameter( params, tag, val ):
    if tag == 'name': params.name = val
    elif tag == 'version': params.version = val
    elif tag == 'min_loop_length':    params.min_loop_length = int( val )
    elif tag == 'allow_strained_3WJ': params.allow_strained_3WJ = (val == 'True')
    elif len( tag )>=2 and tag[:2] == 'Kd':
        setup_base_pair_type_by_tag( params, tag, float(val) )
        update_C_eff_stack( params )
    elif len( tag )>=11 and tag[:11] == 'C_eff_stack':
        if tag == 'C_eff_stacked_pair':
            for bpt1 in params.base_pair_types:
                for bpt2 in params.base_pair_types:
                    params.C_eff_stack[bpt1][bpt2] = float(val)
        else:
            assert( len(tag) > 11 )
            tags = tag[12:].split('_')
            assert( len( tags ) == 2 )
            bpts1 = get_base_pair_types_for_tag( params, tags[0] )
            bpts2 = get_base_pair_types_for_tag( params, tags[1] )
            for bpt1 in bpts1:
                for bpt2 in bpts2:
                    params.C_eff_stack[ bpt1 ][ bpt2 ] = float(val)
                    params.C_eff_stack[ bpt2.flipped ][ bpt1.flipped ] = float(val)
    else:
        setattr( params, tag, float( val ) )

def read_params_fields( params_file ):
    assert( os.path.isfile( params_file ) )
    lines = open( params_file, 'r' ).readlines()
    tags = []
    vals = []
    for line in lines:
        if len( line.replace( ' ','' ) ) == 1: continue # blank line
        elif line[0] == '#': continue # comment line
        else:
            cols = line.split()
            tags.append( cols[0] )
            vals.append( cols[1] )
            if len( cols ) > 2: assert( cols[2][0] == '#' ) # better be a comment
    return zip( tags, vals )

def get_params_from_file( params_file_tag ):
    params = AlphaFoldParams()
    params_file = os.path.dirname( os.path.abspath(__file__) ) + '/parameters/'+params_file_tag +'.params'
    params_fields = read_params_fields( params_file );
    for param_tag,param_val in params_fields:  set_parameter( params, param_tag, param_val )
    return params
