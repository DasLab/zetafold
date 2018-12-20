from __future__ import print_function
import math
from .constants import KT_IN_KCAL
from .assert_equal import assert_equal
import sys

def _show_results( self ):
    fid = sys.stdout
    write_string_with_spaces( self.sequence, self.ligated, fid )
    fid.write(' sequence\n')
    if self.structure != None:
        write_string_with_spaces( self.structure, self.ligated, fid )
        fid.write(' input structure\n')
    if self.struct_MFE != None:
        write_string_with_spaces( self.struct_MFE, self.ligated, fid )
        fid.write(' (pseudo)MFE\n')
    if len( self.struct_stochastic ) > 0:
        for struct in self.struct_stochastic:
            write_string_with_spaces( struct, self.ligated, fid )
            fid.write(' stochastic\n')
    if len( self.struct_enumerate ) > 0:
        for struct in self.struct_enumerate:
            write_string_with_spaces( struct, self.ligated, fid )
            fid.write(' enumerate\n')
    print('Z =',self.Z)
    print('dG =',-KT_IN_KCAL * math.log( self.Z ))
    print()
    if self.deriv_params:
        show_derivs( self.deriv_params, self.log_derivs )
    if self.bpp and not self.suppress_bpp_output:
        output_bpp_matrix( self )
        output_bpp_plot( self )

def write_string_with_spaces( sequence, ligated, fid ):
    for i in range( len( sequence) ):
        fid.write( sequence[i] )
        if not ligated[i]: fid.write( ' ' )
    # to denote circularized
    if ligated[ len(sequence)-1 ]:
        fid.write( '*' )
        fid.write( ' [circularized]' )


def output_bpp_matrix( self ):
    bpp_file = 'bpp.txt'
    fid = open( bpp_file, 'w' )
    for i in range( self.N ):
        for j in range( self.N ):
            fid.write(' %25.12f' % self.bpp[i][j] )
        fid.write('\n')
    fid.close()
    print( 'Outputted base pair probability matrix  to: ', bpp_file )

def output_bpp_plot( self ):
    import matplotlib.pyplot as plt
    import seaborn as sns
    f, ax = plt.subplots(dpi=50)
    sns.heatmap( self.bpp, linewidths=0.1,square=True, vmin=0, vmax=1,ax=ax)
    bpp_fig_file = 'bpp.png'
    plt.savefig( bpp_fig_file )
    print( 'Outputted base pair probability heatmap to: ', bpp_fig_file )

def show_derivs( deriv_params, log_derivs ):
    print( '%20s %25s' % ('parameter','d(logZ)/d(log parameter)' ) )
    for i,parameter in enumerate(deriv_params):
           print( '%20s %25.12f' % (parameter, log_derivs[i] ) )
    print()
    # too much verbiage if we also include dZ/dparameter:
    #print( '%20s %25s %25s' % ('parameter','d(logZ)/d(log parameter)','dZ/dparameter' ) )
    #for i,parameter in enumerate(deriv_params):
    #       print( '%20s %25.12f %25.8f' % (parameter, log_derivs[i], derivs[i] ) )

def _show_matrices( self ):
    output_DP( "Z_BP", self.Z_BP )
    output_DP( "Z_cut", self.Z_cut )
    output_DP( "C_eff_basic", self.C_eff_basic )
    output_DP( "C_eff", self.C_eff, self.Z_final )
    #output_DP( "dC_eff", self.dC_eff, self.dZ_final )
    output_DP( "Z_coax", self.Z_coax )
    output_DP( "Z_linear", self.Z_linear )
    output_square( "BPP", self.bpp );

def output_DP( tag, X, X_final = []):
    N = len( X )
    print()
    print("-----", tag, "-----")
    for i in range( N ):
        for q in range( i ): print('          ', end='')# padding to line up)
        for j in range( N ):
            print(' %9.3f' % X.val(i,i+j), end='')
        if len( X_final ) > 0: print('==> %9.3f' % X_final.val(i), end='')
        print()

def output_square( tag, X ):
    N = len( X )
    print()
    print("-----", tag, "-----")
    for i in range( N ):
        for j in range( N ):
            print(' %9.3f' % X[i][j],end='')
        print()

def output_test( p, Z_ref = 0, bpp_idx= [], bpp_expected = 0,  deriv_parameters = None, log_derivs_ref = None, ):
    print('Z =',p.Z,' [calc]' )
    print('Z =',Z_ref,' [expected]')
    assert_equal( p.Z, Z_ref )
    print('bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),p.bpp[ bpp_idx[0] ][ bpp_idx[1] ], ' [calc]')
    print('bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),bpp_expected,' [expected]')
    assert_equal( p.bpp[ bpp_idx[0] ][ bpp_idx[1] ], bpp_expected )

    if deriv_parameters != None:
        print()
        print( 'd(logZ)/d(log parameter)' )
        log_derivs = p.get_log_derivs( deriv_parameters )
        for i,parameter in enumerate(deriv_parameters): print( parameter,':', log_derivs[i],'[calc] ',log_derivs_ref[i],'[expected]' )
        for log_deriv,log_deriv_ref in zip(log_derivs,log_derivs_ref):  assert_equal( log_deriv, log_deriv_ref )

    print()
