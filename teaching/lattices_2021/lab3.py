from fpylll import *
from copy import copy

inp = open("hidden.txt", "r")
p = int(inp.readline())
#print(p)

t = []
msbs = []
for data in inp:
    lst = data.split(" ")
    assert(len(lst)==2)
    t.append(int(lst[0]))
    msbs.append(int(lst[1].rstrip()))
#print(t)
#print(msbs)
target = msbs+[0]


def build_basis(t, msbs, p):
    n = len(t)
    A = IntegerMatrix(n,n)
    for i in range(n-1):
        A[i,i] = p
        A[n-1, i] = t[i]
        #A[n, i] = msbs[n+1]
    A[n-1, n-1] = 1
    return A

A = build_basis(t, msbs, p)
Acopy = copy(A)

#print(A)
U = IntegerMatrix.identity(A.nrows)
UInv = IntegerMatrix.identity(A.nrows)
M = GSO.Mat(Acopy, float_type='mpfr', U=U, UinvT=UInv)
M.update_gso()

L = LLL.Reduction(M, delta = 0.98, eta=0.52)
L()

#print(UInv.transpose()*Acopy)

M.update_gso()
error = M.babai(target)
n = len(t)
coeff = IntegerMatrix.from_iterable(1, n, error)*Acopy
print(coeff)
