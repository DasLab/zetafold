from zetafold.util.wrapped_array import WrappedArray

class DynamicProgrammingMatrix:
    '''
    Dynamic Programming 2-D Matrix that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    '''
    def __init__( self, N, val = 0.0, diag_val = 0.0, DPlist = None, update_func = None, options = None, name = None ):
        self.N = N
        self.data = WrappedArray(N)
        for i in range( N ):
            self.data[i] = WrappedArray( N )
            for j in range( N ):
                self.data[i][j] = DynamicProgrammingData( val, options = options )
                self.data[i][j].info.append( (self,i,j) )

        for i in range( N ): self.data[i][i].Q = diag_val
        if DPlist != None: DPlist.append( self )
        self.options = options
        self.update_func = update_func

        self.backtrack_info_updated = [None]*N
        for i in range( N ): self.backtrack_info_updated[i] = [False]*N
        self.name = name

    def __getitem__( self, idx ):
        return self.data[ idx ]

    def __len__( self ): return len( self.data )

    def val( self, i, j ): return self.data[i][j].Q
    def set_val( self, i, j, val ): self.data[i][j].Q = val

    def update( self, partition, i, j ):
        self.data[ i ][ j ].zero()
        self.update_func( partition, i, j )

    def get_backtrack_info( self, partition, i, j ):
        if not self.backtrack_info_updated[i][j]:
            partition.options.calc_backtrack_info = True
            self.update( partition, i, j )
            partition.options.calc_backtrack_info = False
            self.backtrack_info_updated[i][j] = True
        return self.data[i][j].backtrack_info

class DynamicProgrammingList:
    '''
    Dynamic Programming 1-D list that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    Used for Z_final
    '''
    def __init__( self, N, val = 0.0, update_func = None, options = None, name = None ):
        self.N = N
        self.data = WrappedArray( N )
        for i in range( N ):
            self.data[i] = DynamicProgrammingData( val, options = options )
        self.options = options
        self.update_func = update_func
        self.backtrack_info_updated = [False]*N
        self.name = name

    def __getitem__( self, idx ):
        return self.data[ idx ]

    def __setitem__( self, idx, val ):
        self.data[ idx ] = val

    def __len__( self ): return len( self.data )

    def val( self, i ): return self.data[i].Q

    def get_backtrack_info( self, partition, i ):
        if not self.backtrack_info_updated[i]:
            partition.options.calc_backtrack_info = True
            self.update( partition, i )
            partition.options.calc_backtrack_info = False
            self.backtrack_info_updated[i] = True
        return self.data[i].backtrack_info

    def update( self, partition, i ):
        self.data[ i ].zero()
        self.update_func( partition, i )

class DynamicProgrammingData:
    '''
    Dynamic programming object, with derivs and contribution accumulation.
     Q   = value
     contrib = contributions

    NB: dQ  = derivative (removed in branch rhiju/remove_inline_derivatives)
    '''
    def __init__( self, val = 0.0, options = None ):
        self.Q = val
        self.backtrack_info = []
        self.info = []
        self.options = options

    def zero( self ):
        self.Q = 0.0
        self.backtrack_info = []

    def __iadd__(self, other):
        if other.Q == 0.0: return self
        self.Q  += other.Q
        if self.options and self.options.calc_backtrack_info:
            if len( other.info ) > 0: self.backtrack_info.append( [other.Q, other.info] )
        return self

    def __mul__(self, other):
        prod = DynamicProgrammingData()
        prod.options = self.options
        if not prod.options: prod.options = other.options
        if isinstance( other, DynamicProgrammingData ):
            prod.Q  = self.Q * other.Q
            if self.options and self.options.calc_backtrack_info:
                info = self.info + other.info
                if len( info ) > 0:
                    prod.backtrack_info = [ [ prod.Q, info ] ]
                    prod.info = info
        else:
            prod.Q  = self.Q * other
            if self.options and self.options.calc_backtrack_info:
                for contrib in self.backtrack_info:
                    prod.backtrack_info.append( [contrib[0]*other, contrib[1] ] )
                prod.info = self.info
        return prod

    def __truediv__( self, other ):
        return self * ( 1.0/other )

    __rmul__ = __mul__
    __floordiv__ = __truediv__
    __div__ = __truediv__


class WrappedArray:
    '''
    For all the various cross-checks, like equality of partition function starting at any
     i and wrapping around to N and then back 1 and then i-1, need to keep applying modulo N.
    '''
    def __init__( self, N, val = None ):
        self.data = [val] * N
        self.N = N
    def __getitem__( self, idx ):
        return self.data[idx % self.N]
    def __setitem__( self, idx, item ):
        self.data[idx % self.N] = item
    def __len__( self ):
        return self.N

