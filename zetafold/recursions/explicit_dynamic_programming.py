#
# Much simpler (less intelligent) object for dynamic programming than in dynamic_programming.py --
#  forces code to explicitly figure out updates to values, derivatives, and contributions
#
class DynamicProgrammingMatrix:
    '''
    Dynamic Programming 2-D Matrix that automatically:
      knows how to update values at i,j
    '''
    def __init__( self, N, val = 0.0, diag_val = 0.0, DPlist = None, update_func = None, options = None, name = None ):
        self.N = N

        self.Q = [None]*N
        for i in range( N ): self.Q[i] = [val]*N
        for i in range( N ): self.Q[i][i] = diag_val

        self.backtrack_info = [None]*N
        for i in range( N ):
            self.backtrack_info[i] = []
            for j in range( N ): self.backtrack_info[i].append( [] )

        self.backtrack_info_updated = [None]*N
        for i in range( N ): self.backtrack_info_updated[i] = [False]*N

        if DPlist != None: DPlist.append( self )
        self.update_func = update_func

        self.name = name

    def val( self, i, j ): return self.Q[i%self.N][j%self.N]
    def set_val( self, i, j, val ): self.Q[i%self.N][j%self.N] = val

    def update( self, partition, i, j ):
        self.Q[ i ][ j ] = 0
        self.backtrack_info[ i ][ j ] = []
        self.update_func( partition, i, j )

    def get_backtrack_info( self, partition, i, j ):
        if not self.backtrack_info_updated[i][j]:
            partition.options.calc_backtrack_info = True
            self.update( partition, i, j )
            partition.options.calc_backtrack_info = False
            self.backtrack_info_updated[i][j] = True
        return self.backtrack_info[i][j]

    def __len__( self ):
        return len( self.Q )

class DynamicProgrammingList:
    '''
    Dynamic Programming 1-D list that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    Used for Z_final
    '''
    def __init__( self, N, val = 0.0, update_func = None, options = None, name = None ):
        self.N = N
        self.Q = [ val ]*N
        self.backtrack_info = [None] * N
        for i in range( N ): self.backtrack_info[i] = []
        self.backtrack_info_updated = [False]*N
        self.update_func = update_func
        self.name = name

    def __len__( self ): return self.N

    def val( self, i ): return self.Q[i]

    def update( self, partition, i ):
        self.Q[ i ] = 0.0
        self.backtrack_info[ i ] = []
        self.update_func( partition, i )

    def get_backtrack_info( self, partition, i ):
        if not self.backtrack_info_updated[i]:
            partition.options.calc_backtrack_info = True
            self.update( partition, i )
            partition.options.calc_backtrack_info = False
            self.backtrack_info_updated[i] = True
        return self.backtrack_info[i]
