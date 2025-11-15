from fpylll import *
from g6k import *
from g6k.utils.stats import SieveTreeTracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from fpylll.util import gaussian_heuristic
from math import sqrt


#FPLLL.set_random_seed(1337)



#-------------------------------------------------------
# SVP
#-------------------------------------------------------
"""
n = 60
A = IntegerMatrix.random(n, "qary", k=30, bits=30)
g6k = Siever(A, seed=0x1337) # g6k has its own seed

g6k.initialize_local(0, 0, n) # sieve the entire basis
g6k(alg="bgj1") #alg \in {"nv", "gauss", "hk3"}
v = g6k.best_lifts()[0][2]
print('best lifts', g6k.best_lifts())
w = A.multiply_left(v)  # w = v*A
print('norms:', sum(w_**2 for w_ in w), A[0].norm()**2)

"""
#-------------------------------------------------------
# Simulate BKZ
#-------------------------------------------------------
"""
n = 60
A = IntegerMatrix.random(n, "qary", k=30, bits=30)

g6k = Siever(A, seed=0x1337) # g6k has its own seed
g6k.initialize_local(0, 0, n)
sp = g6k.params.new()
sp.threads = 1
#g6k(alg="bgj1")

tracer = SieveTreeTracer(g6k, root_label=("bgj1"), start_clocks=True)
for b in (20, 30, 40, 50, 60): # progressive BKZ
    pump_n_jump_bkz_tour(g6k, tracer,
                         b, pump_params={"down_sieve": True})
    print(b, A[0].norm()**2)

"""

#-------------------------------------------------------
# many short vectors
#-------------------------------------------------------

n = 44
A = IntegerMatrix.random(n, "qary", k=n/2, bits=10)
g6k = Siever(A)
g6k.lll(0, n)

# Run a progressive sieve (preprocessing)
g6k.initialize_local(0, n/2, n)
while g6k.l > 0:
    # Extend the lift context to the left
    g6k.extend_left(1)
    # Sieve
    g6k()

with g6k.temp_params(saturation_ratio=.95, saturation_radius=1.4,
                     db_size_base=sqrt(1.7), db_size_factor=3):
    g6k()


gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])

db = list(g6k.itervalues())
found = 0

for x in db:
    v = A.multiply_left(x)
    l = sum(v_**2 for v_ in v)
    if l < 1.2 * gh:
        print(l/gh, v)
        found += 1
