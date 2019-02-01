import numpy as np
import argparse

class MEA:
    def __init__(self, bpps, gamma = 1.0, debug=False):
        self.debug = debug
        self.bpps = bpps
        self.N=self.bpps.shape[0]
        self.gamma = gamma
        self.W = np.zeros([self.N,self.N])
        self.MEA_bp_list = []
        self.structure = ['.']*self.N
        self.MEA_bp_matrix = np.zeros([self.N,self.N])
        self.tb = np.zeros([self.N,self.N])
        self.min_hp_length=3
        self.evaluated = False
        
    def fill_W(self, i, j):
        options = [self.W[i+1, j], self.W[i, j-1], (self.gamma+1)*self.bpps[i,j] + self.W[i+1, j-1] - 1,\
                np.max([self.W[i,k] + self.W[k+1, j] for k in range(i+1,j)])]
        self.W[i,j] = np.max(options) 
        self.tb[i,j] = np.argmax(options) #0: 5' pass, 1: 3' pass, 2: bp, 3: multiloop
        
    def run_MEA(self):
        # fill weight matrix
        for length in range(self.min_hp_length, self.N):
            for i in range(self.N-length):
                j = i + length
                self.fill_W(i,j)
                
        self.traceback(0,self.N-1)
        
        for x in self.MEA_bp_list:
            self.MEA_bp_matrix[x[0],x[1]]=1
            self.structure[x[0]]='('
            self.structure[x[1]]=')'
        
        self.structure = ''.join(self.structure)
        if not self.evaluated: self.evaluated = True
        
    def traceback(self, i, j):
        if j <= i:
            return
        elif self.tb[i,j] == 0: #5' neighbor
            if self.debug: print(i,j, "5'")
            self.traceback(i+1,j)
        elif self.tb[i,j] == 1: #3' neighbor
            if self.debug: print(i,j, "3'")
            self.traceback(i,j-1)
        elif self.tb[i,j] == 2: # base pair
            if self.debug: print(i,j,'bp')
            self.MEA_bp_list.append((i,j))
            self.traceback(i+1,j-1)
        else: #multiloop
            for k in range(i+1,j):
                if self.W[i,j] == self.W[i, k] + self.W[k+1,j]:
                    if self.debug: print(i,j,"multiloop, k=",k)
                    self.traceback(i,k)
                    self.traceback(k+1,j)
                    break

    def score(self, struct):

        if not self.evaluated: self.run_MEA()

        true_m = struct[np.triu_indices(self.N)]
        est_m = self.MEA_bp_matrix[np.triu_indices(self.N)]

        TP, FP, cFP, TN, FN = 0,0,0,0,0 
        for i in range(len(true_m)):
            if true_m[i] == 1:
                if est_m[i] == 1: 
                    TP += 1
                else:
                    FN += 1
            elif true_m[i] == 0:
                if est_m[i] == 0: 
                    TN += 1
                else: 
                    FP +=1
                    #check for compatible false positive
                    a,b = np.triu_indices(self.N)
                    if np.sum(struct,axis=0)[a[i]]+ np.sum(struct,axis=0)[b[i]]==0:
                       cFP +=1
            else:
                pass
 
        return TP,FP, cFP, TN,FN
 
def convert_dotbracket_to_matrix(s):
    triu_inds = np.triu_indices(len(s))
    m = np.zeros([len(s),len(s)])
    bp1=[]
    bp2=[]
    for i, char in enumerate(s):
        if char=='(':
            bp1.append(i)
        if char==')':
            bp2.append(i)
    for i in list(reversed(bp1)):
        for j in bp2:
            if j > i:
                m[i,j]=1.0
                bp2.remove(j)
                break
    return m[triu_inds]

def convert_matrix_to_dotbracket(m):
    pass

if __name__ == '__main__':
    parser=argparse.ArgumentParser(
        description='''MEA (Maximum Expected Accuracy) structure prediction, using CentroidFold algorithm (Hamada 2009).''')
    parser.add_argument('--bp_matrix','-p', help='NxN matrix of bp probabilities.')
    parser.add_argument('--true_struct','-s', help='text file containing true structure')

    args = parser.parse_args()
    bp_matrix=np.loadtxt("./%s" % args.bp_matrix)
    #struct=open("./%s" % args.true_struct,'r').read()
    struct = np.loadtxt(args.true_struct)
    print('Gamma\tstruct\tTP\tFP\tcFP\tTN\tFN')
    for g in range(-7,7):
        cls = MEA(bp_matrix, gamma=2**g)
        cls.run_MEA()
        TP, FP, cFP, TN, FN = cls.score(struct[:cls.N])
        print('%d\t%.3f\t%.3f\t%.3f\t%.3f' % (g, TP,FP,cFP, TN,FN))

