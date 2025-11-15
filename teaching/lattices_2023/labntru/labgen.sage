from fpylll import *
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
FPLLL.set_precision( 100 )

def str_to_ternary(s, max_len=128): #кодированиие буквенного сообщения
    l = []
    for c in s:
        n = ord(c)%3**5
        for i in range(5):
            l.append( n%3-1 )
            n = n//3
    assert len(l) < max_len, "Too long message!"
    if len(l)<max_len:
        l += [ randrange(-1,2) for i in range(max_len-len(l)) ]
    return l

def ternary_to_str( l ):  #декодированиие буквенного сообщения
    s = []
    for i in range(len(l)//5):
        n=0
        for j in range(5):
            n+=(l[5*i+j]+1)*3^j
        s.append(chr(n))
    return "".join( [ ss for ss in s ] )

def convolution(f,g):
    return (f * g).mod(modulo)

def balancedmod(f,q):
    g = list(((f[i] + q//2) % q) - q//2 for i in list(range(N)))
    return ZZ[x](g)

def gen_invertable_f(x,d,N):  #генерация обратимого f
    Rq=Integers(q)[x].quotient_ring(modulo)
    baseq=Rq.base_ring()

    #generated poly in T(d+1,d)
    fq=Rq(0)
    stop_iter=d^2
    nump, numm = 0, 0 #КОЛИЧЕСТВО +1 И -1 В КОЭФФ
    i=0
    while i<=d:
        coord=randrange(N)
        #если коэффициент свободен, то ставим его как 1
        if fq[coord]==baseq(0):
            nump+=1
            fq+=Rq(x^coord)
        else:
            continue
        stop_iter-=1
        i+=1
        if stop_iter==0:
            print('Error!')
            break
    i=0
    while i<d:
        coord=randrange(N)
        #если коэффициент свободен, то ставим его как -1
        if fq[coord]==baseq(0):
            numm+=1
            fq-=Rq(x^coord)
        else:
            continue
        i+=1
        stop_iter-=1
        if stop_iter==0:
            print('Error!')
            break
    print('+:', nump, '-:', numm)
    return balancedmod(fq,q)

def gen_g(x,d,N):
    Rq=Integers(q)[x].quotient_ring(modulo)
    baseq=Rq.base_ring()

    #generated poly in T(d+1,d)
    fq=Rq(0)
    stop_iter=d^2
    nump, numm = 0, 0 #КОЛИЧЕСТВО +1 И -1 В КОЭФФ
    i=0
    while i<d:
        coord=randrange(N)
        #если коэффициент свободен, то ставим его как 1
        if fq[coord]==baseq(0):
            nump+=1
            fq+=Rq(x^coord)
        else:
            continue
        stop_iter-=1
        i+=1
        if stop_iter==0:
            print('Error!')
            break
    i=0
    while i<d:
        coord=randrange(N)
        #если коэффициент свободен, то ставим его как -1
        if fq[coord]==baseq(0):
            numm+=1
            fq-=Rq(x^coord)
        else:
            continue
        i+=1
        stop_iter-=1
        if stop_iter==0:
            print('Error!')
            break
    print('+:', nump, '-:', numm)
    return balancedmod(fq,q)


def invertmodprime(f,p): #обратить f по модулю p
    Zx.<x> = ZZ[]
    T = Zx.change_ring(Integers(p)).quotient(modulo)
    return Zx(lift(1 / T(f)))

def encrypt( modulo,h,q,msg ):
    print("mod", modulo)
    N=int( modulo.degree(x) )

    num = str_to_ternary(msg, max_len=N)
    msg = ZZ[x](num)
    print("Message:",msg)

    r=gen_g(x,d,N)

    enc=balancedmod(convolution(h,r) + msg,q)
    print("Cyphertext:",enc)
    return enc

def decrypt( modulo, f, q, p, enc ):
    a = balancedmod(convolution(enc,f),q)
    m=balancedmod(convolution(a,fp),p)
    return m
