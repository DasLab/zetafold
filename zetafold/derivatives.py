from .base_pair_types import get_base_pair_type_for_tag, get_base_pair_types_for_tag

def _get_log_derivs( self, deriv_parameters = [] ):
    '''
    Output

       d( log Z )/ d( log parameter )

    using simple expressions that require O( N^2 ) time or less after
    the original O( N^3 ) dynamic programming calculations
    '''
    if deriv_parameters == None: return None
    if deriv_parameters == []:
        for tag in self.params.parameter_tags: deriv_parameters.append( tag )

    derivs = [None]*len(deriv_parameters)

    # Derivatives with respect to each Kd
    N = self.N
    Z = self.Z_final.val( 0 )
    for n,parameter in enumerate(deriv_parameters):
        if parameter == 'l':
            # Derivatives with respect to loop closure parameters
            num_internal_linkages = 0.0
            for i in range( N ):
                if not self.ligated[i]: continue
                num_internal_linkages += self.params.l * self.C_eff_no_coax_singlet.val( i+1, i ) / self.params.C_std / Z
            derivs[ n ] = num_internal_linkages
        elif parameter == 'l_BP':
            num_base_pairs_closed_by_loops = 0.0
            for i in range( N ):
                for j in range( N ):
                    if ( j - i ) % N < 2: continue
                    num_base_pairs_closed_by_loops += self.params.l**2 * self.params.l_BP * self.C_eff.val(i+1,j-1) * self.Z_BP.val(j,i) / Z
            derivs[ n ] = num_base_pairs_closed_by_loops
        elif parameter == 'C_init':
            num_closed_loops = get_bpp_tot( self ) - self.num_strand_connections()
            derivs[ n ] = num_closed_loops
        elif len(parameter)>=2 and  parameter[:2] == 'Kd':
            if parameter == 'Kd':
                # currently can only handle case where Kd controls *all* of the base pair types
                for base_pair_type in self.params.base_pair_types: assert( base_pair_type.Kd == self.params.base_pair_types[0].Kd )
                derivs[ n ] = - get_bpp_tot( self )
            else:
                Kd_tag = parameter[3:]
                derivs[ n ] = - get_bpp_tot_for_base_pair_type( self, get_base_pair_type_for_tag( self.params, Kd_tag ) )
        elif len(parameter)>=11 and parameter[:11] == 'C_eff_stack':
            # Derivatives with respect to motifs (stacked pairs first)
            if parameter == 'C_eff_stacked_pair':
                motif_prob = 0.0
                for base_pair_type in self.params.base_pair_types:
                    for base_pair_type2 in self.params.base_pair_types:
                        # currently can only handle case where C_eff_stacked_pair controls *all* of the stacked pairs
                        assert( self.params.C_eff_stack[base_pair_type][base_pair_type2] == self.params.C_eff_stack[self.params.base_pair_types[0]][self.params.base_pair_types[0]])
                        motif_prob += get_motif_prob( self, base_pair_type, base_pair_type2 )
                derivs[n] = motif_prob
            else:
                assert( len(parameter) > 11 )
                tags = parameter[12:].split('_')
                assert( len( tags ) == 2 )
                bpts1 = get_base_pair_types_for_tag( self.params, tags[0] )
                bpts2 = get_base_pair_types_for_tag( self.params, tags[1] )
                derivs[ n ] = 0.0
                for bpt1 in bpts1:
                    for bpt2 in bpts2:
                        derivs[ n ] += get_motif_prob( self, bpt1, bpt2 )
        elif parameter == 'K_coax':
            coax_prob = 0.0
            C_eff_for_coax = self.C_eff if self.params.allow_strained_3WJ else self.C_eff_no_BP_singlet
            for i in range( N ):
                for j in range( N ):
                    coax_prob += self.Z_coax.val(i,j) * self.params.l_coax * self.params.l**2 * C_eff_for_coax.val(j+1,i-1) / Z
                    coax_prob += self.Z_coax.val(i,j) * self.Z_cut.val(j,i) / Z
            derivs[ n ] = coax_prob
        elif parameter == 'l_coax':
            coax_prob = 0.0
            C_eff_for_coax = self.C_eff if self.params.allow_strained_3WJ else self.C_eff_no_BP_singlet
            for i in range( N ):
                for j in range( N ):
                    coax_prob += self.Z_coax.val(i,j) * self.params.l_coax * self.params.l**2 * C_eff_for_coax.val(j+1,i-1) / Z
            derivs[ n ] = coax_prob
        else:
            print "Did not recognize parameter ", parameter
            pass

    return derivs

def get_bpp_tot_for_base_pair_type( self, base_pair_type ):
    assert( self.calc_all_elements )
    bpp = 0.0
    N = self.N
    for i in range( N ):
        for j in range( N ):
            bpp += self.Z_BPq[base_pair_type].val(i,j) * self.Z_BPq[base_pair_type.flipped].val(j,i) * base_pair_type.Kd / self.Z_final.val(0)
    return bpp

def get_bpp_tot( self ):
    bpp_tot = []
    for base_pair_type in self.params.base_pair_types: bpp_tot.append( get_bpp_tot_for_base_pair_type( self, base_pair_type ) )
    return sum( bpp_tot ) / 2.0

def get_motif_prob( self, base_pair_type, base_pair_type2 ):
    # base pair forms a stacked pair with previous pair
    #
    #      bp2
    #  i+1 ... j-1
    #    |     |
    #    i ... j
    #      bp1
    #
    motif_prob = 0.0
    Z_BPq1 = self.Z_BPq[base_pair_type.flipped]
    Z_BPq2 = self.Z_BPq[base_pair_type2]
    N = self.N
    for i in range( N ):
        for j in range( N ):
            if ( j - i ) % N < 3: continue
            if not self.ligated[i]: continue
            if not self.ligated[(j-1)%N]: continue
            if Z_BPq1.val(j  ,  i) == 0: continue
            if Z_BPq2.val(i+1,j-1) == 0: continue
            motif_prob += self.params.C_eff_stack[base_pair_type][base_pair_type2] * Z_BPq1.val(j,i) * Z_BPq2.val(i+1,j-1) / self.Z_final.val(0) / 2.0
    return motif_prob
