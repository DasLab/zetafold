from __future__ import print_function
import math
from .base_pair_types import BasePairType, setup_base_pair_type, get_base_pair_types_for_tag, get_base_pair_type_for_tag
from .util.constants import KT_IN_KCAL
import glob
import os.path

class AlphaFoldParams:
    '''
    Parameters that define the statistical mechanical model for RNA folding
    '''
    def __init__( self ):
        self.C_std  = 1.0      # 1 M. drops out in end (up to overall scale factor).
        self.parameter_tags   = [] # K_CG, etc.
        self.parameter_values = [] # floats
        self.string_tags   = [] # name, version, etc.
        self.string_values = [] # strings

    def get_variables( self ):
        if self.C_init == 0.0 and self.name == 'empty': print('WARNING! C_init not defined, and params appear empty. Look at get_params() for examples')
        return ( self.C_init, self.l, self.l_BP, self.K_coax, self.l_coax, self.C_std, self.min_loop_length, self.allow_strained_3WJ )

    def set_parameter( self, tag, val ):
        val = _set_parameter( self, tag, val )

    def get_parameter_value( self, param_tag ):
        if self.parameter_tags.count( param_tag ) == 0: return None
        return self.parameter_values[ self.parameter_tags.index( param_tag ) ]

    def check_C_eff_stack( self ): _check_C_eff_stack( self )

    def show_parameters( self ):
        print( '%25s %12s %12s' % ('Parameter','log val','val')  )
        for tag, val in zip( self.parameter_tags, self.parameter_values ): print( '%25s %12.7f %12.6f' % (tag,math.log(val),val) )

    def output_to_file( self, params_file ): _output_to_file( self, params_file )

def get_params( params = None, suppress_all_output = False ):
    '''
    master function to get parameters
    '''
    params_object = None
    if isinstance(params,AlphaFoldParams): return params
    elif params == None or params =='':    params_object = get_latest_params()
    else:
        assert( isinstance( params, str ) )
        params_object = get_params_from_file( params )
    if not suppress_all_output: print('Parameters: ', params_object.name, ' version', params_object.version)
    return params_object

def read_params_fields( params_file ):
    '''
    simply read lines like
      name minimal
    i.e.., tag then value from file.
    '''
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
    '''
    find the file (if it exists) and then load up into params variables.
    '''
    params = AlphaFoldParams()
    params_file = params_file_tag
    if not os.path.exists( params_file ): params_file = params_file_tag +'.params'
    if not os.path.exists( params_file ): params_file = os.path.dirname( os.path.abspath(__file__) ) + '/parameters/'+params_file_tag +'.params'
    if not os.path.exists( params_file ): params_file = os.path.dirname( os.path.abspath(__file__) ) + '/parameters/zetafold_'+params_file_tag +'.params'
    if not os.path.exists( params_file ):
        print()
        print( 'Could not find requested parameters:', params_file_tag )
        print( 'Options are: ' )
        for params_file in get_all_params_files(): print('  ',params_file)
        print()
        return None
    params_fields = read_params_fields( params_file );
    for param_tag,param_val in params_fields:
        val = _set_parameter( params, param_tag, param_val )

    return params

def get_latest_params():
    '''
    look for parameters/zetafold_v*.*.params and choose latest
    '''
    params_dir =  os.path.dirname( os.path.abspath(__file__) ) + '/parameters/'
    params_files = glob.glob( params_dir+'zetafold*.params' )
    params_files.sort()
    return get_params_from_file( params_files[-1] )

def get_all_params_files():
    '''
    list of all parameters/*.params
    '''
    params_dir =  os.path.dirname( os.path.abspath(__file__) ) + '/parameters/'
    params_files = glob.glob( params_dir+'*.params' )
    params_files.sort()
    return [ os.path.basename(x).replace('.params','') for x in params_files ]

def _output_to_file( self, params_file ):
    print( 'Outputting parameters to file: ', params_file )
    fid = open( params_file, 'w' )
    for tag,val in zip(self.string_tags,self.string_values): fid.write( '%s %s\n' % (tag, val) )
    for tag,val in zip(self.parameter_tags,self.parameter_values): fid.write( '%s %12.6f\n' % (tag, val) )
    fid.close()

#############################################################################################################
#  Following has basic logic for setting parameter values based on input parameter -- note that
#     in some cases, like C_eff_stacked_pair, multiple internal parameters need to get updated
#############################################################################################################
def _set_parameter( self, tag, val ):
    '''
    Based on tag, one or multiple parameters might be set.
    Returns float( val ) if an actual continuous parameter is set; None otherwise.
    '''
    float_parameter = False
    if tag == 'name': self.name = val
    elif tag == 'version': self.version = val
    elif tag == 'min_loop_length':    self.min_loop_length = int( val )
    elif tag == 'allow_strained_3WJ': self.allow_strained_3WJ = (val == 'True')
    elif len( tag )>=2 and tag[:2] == 'Kd':
        setup_base_pair_type_by_tag( self, tag, float(val) )
        update_C_eff_stack( self )
        float_parameter = True
    elif len( tag )>=11 and tag[:11] == 'C_eff_stack':
        if tag == 'C_eff_stacked_pair':
            for bpt1 in self.base_pair_types:
                for bpt2 in self.base_pair_types:
                    self.C_eff_stack[bpt1][bpt2] = float(val)
        else:
            assert( len(tag) > 11 )
            tags = tag[12:].split('_')
            assert( len( tags ) == 2 )
            bpts1 = get_base_pair_types_for_tag( self, tags[0] )
            bpts2 = get_base_pair_types_for_tag( self, tags[1] )
            for bpt1 in bpts1:
                for bpt2 in bpts2:
                    self.C_eff_stack[ bpt1 ][ bpt2 ] = float(val)
                    self.C_eff_stack[ bpt2.flipped ][ bpt1.flipped ] = float(val)
        float_parameter = True
    else:
        if not tag in ('name','version','C_init','l','l_BP','l_coax','K_coax'):
            print( 'Unrecognized tag!!', tag )
            exit()
        setattr( self, tag, float( val ) )
        float_parameter = True
    if float_parameter:
        if self.parameter_tags.count( tag ) == 0:
            self.parameter_tags.append( tag )
            self.parameter_values.append( None )
        self.parameter_values[ self.parameter_tags.index( tag ) ] = float(val)
    else:
        if self.string_tags.count( tag ) == 0:
            self.string_tags.append( tag )
            self.string_values.append( None )
        self.string_values[ self.string_tags.index( tag ) ] = val

def update_C_eff_stack( params, val = None ):
    if not hasattr( params, 'C_eff_stack' ): params.C_eff_stack = {}
    for bpt1 in params.base_pair_types:
        if not params.C_eff_stack.has_key( bpt1 ):  params.C_eff_stack[ bpt1 ] = {}
        for bpt2 in params.base_pair_types:
            if not params.C_eff_stack[ bpt1 ].has_key( bpt2 ): params.C_eff_stack[ bpt1 ][ bpt2 ] = None
            if val != None: params.C_eff_stack[ bpt1 ][ bpt2 ] = val

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
    base_pair_type = get_base_pair_type_for_tag( params, tag )
    if base_pair_type != None:
        base_pair_type.Kd = val
        base_pair_type.flipped.Kd = val
        return
    if tag == 'matchlowercase':
        setup_base_pair_type( params, '*', '*', val, match_lowercase = True )
    else:
        setup_base_pair_type( params, tag[0], tag[1], val, match_lowercase = False )
    assert( get_base_pair_type_for_tag( params, tag ) )

