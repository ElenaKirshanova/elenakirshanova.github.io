#
# Recommend to run from conda environment with sage and fpylll
#


from fpylll import *
from fpylll.algorithms.bkz2 import BKZReduction
import time



def binomial_dist(eta):
  s = 0
  for i in range(eta):
    s += randrange(2)
  return s


#LWE parameters
q = 3329
n = 100
eta = 1 #Pr[-1]=Pr[1] = 1/4, Pr[0] = 1/2

#LWE samples
A = random_matrix(GF(q), n, n)
s = vector(ZZ,[binomial_dist(2*eta) - eta for _ in range(n)])
e = vector(ZZ,[binomial_dist(2*eta) - eta for _ in range(n)])
b = A*s + e

#Upper bound on the LWE error
target_norm = ceil(sqrt(2*n))

#Basis of L_q(A)
B = matrix.block(ZZ,[[q*matrix.identity(n), matrix.zero(n)],[-A.transpose(),-matrix.identity(n)]])


#Embed the target b
Kannan = B.stack(vector(ZZ,list(b) + n*[0]))
Kannan = Kannan.transpose().stack( vector(2*n*[0] + [1]) ).transpose()

#Generate the GSO object
M = IntegerMatrix.from_matrix(Kannan)
GSO_M=GSO.Mat(M,float_type="dd")
GSO_M.update_gso()


#Run BKZ
max_block_size = 43
flags = BKZ.AUTO_ABORT|BKZ.MAX_LOOPS|BKZ.GH_BND|BKZ.VERBOSE
bkz = BKZReduction(GSO_M)
for beta in range(2,max_block_size,1):    
    par = BKZ.Param(beta,
                    max_loops=5,
                    flags=flags,
                    strategies=BKZ.DEFAULT_STRATEGY)
    then_round=time.perf_counter()
    bkz(par)
    round_time = time.perf_counter()-then_round
    print(f"BKZ-{beta} finished in {round_time} sec. with shortest vector norm = ", GSO_M.B[0].norm())
    if GSO_M.B[0].norm()<=target_norm:
       break


short = vector(GSO_M.B[0])
e_, s_ = short[:n], short[n:2*n]

assert(e==e_ or e==-e_)
assert(s==s_ or s==-s_)

# print(e)
# print(e_)
# print(s)
# print(s_)

