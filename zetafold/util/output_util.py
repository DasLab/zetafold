from __future__ import print_function
import math
from .constants import KT_IN_KCAL
from .assert_equal import assert_equal
import sys
import gzip

def _show_results( self ):
    fid = sys.stdout
    write_result( 'sequence',self.sequence, self.ligated, fid )
    write_result( 'input structure',self.structure, self.ligated, fid )
    write_result( 'calculate gap structure',self.calc_gap_structure, self.ligated, fid )
    write_result( '(pseudo)MFE',self.struct_MFE, self.ligated, fid )
    write_result( 'stochastic',self.struct_stochastic, self.ligated, fid )
    write_result( 'enumerate',self.struct_enumerate, self.ligated, fid )
    print('Z =',self.Z)
    print('dG (kcal/mol) =',self.dG, ' [full]')
    if self.dG_gap:
        print('dG (kcal/mol) =',self.dG_gap + self.dG, ' [input structure]' )
        print('dG (kcal/mol) =',self.dG_gap, ' [free energy gap]' )
    print()
    if self.deriv_params:
        show_derivs( self.deriv_params, self.log_derivs )
    if self.bpp and not self.suppress_bpp_output:
        output_bpp_matrix( self )
        if not bpp_file: output_bpp_plot( self )

def write_result( tag, variable, ligated, fid ):
    if variable == None: return
    if isinstance( variable, list ):
        for struct in variable: write_result( tag, struct, ligated, fid )
        return
    write_string_with_spaces( variable, ligated, fid )
    fid.write(' %s\n' % tag )

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
    if self.bpp_file: bpp_file = self.bpp_file
    if len(bpp_file)>3 and bpp_file[-3:] == '.gz':
        fid = gzip.open( bpp_file, 'w' )
    else:
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
    if self.bpp: output_square( "BPP", self.bpp );

def output_DP( tag, X, X_final = []):
    N = len( X )
    print()
    print("-----", tag, "-----")
    for i in range( N ):
        for q in range( i ): print('          ', end='')# padding to line up)
        for j in range( N ):
            print(' %9.3f' % X.val(i,i+j), end='')
        if len( X_final ) > 0: print('==> %12.8f' % X_final.val(i), end='')
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
