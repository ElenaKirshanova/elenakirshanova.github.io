Zx.<x> = ZZ['x']

def convolution(f,g,n):
      return (f * g) % (x^n-1)

def balancedmod(f,q,n):
    g = list(((f[i] + q//2) % q) - q//2 for i in range(n))
    return Zx(g)

def randomdpoly(n,d):
    assert d <= n
    result = n*[0]
    for j in range(d):
        while True:
            r = randrange(n)
            if not result[r]: break
        result[r] = 1-2*randrange(2)
    return Zx(result)

def invertmodprime(f,p,n):
    T = Zx.change_ring(Integers(p)).quotient(x^n-1)
    return Zx(lift(1 / T(f)))

def invertmodpowerof2(f,q,n):
    assert q.is_power_of(2)
    g = invertmodprime(f,2,n)
    while True:
        r = balancedmod(convolution(g,f,n),q,n)
        if r == 1: return g
        g = balancedmod(convolution(g,2 - r,n),q,n)

def keypair(n,d,q):
    while True:
        try:
            f = randomdpoly(n,d)
            f3 = invertmodprime(f,3,n)
            fq = invertmodpowerof2(f,q,n)
            break
        except:
            print('fail in generating invertible mod q')
            pass
    g = randomdpoly(n,d)
    publickey = balancedmod(3 * convolution(fq,g,n),q,n)
    secretkey = f,f3
    return publickey,secretkey

def randommessage(n):
    result = list(randrange(3) - 1 for j in range(n))
    return Zx(result)

def encrypt(n,d,q,message,publickey):
    r = randomdpoly(n,d)
    return balancedmod(convolution(publickey,r,n) + message,q,n)

def decrypt(n,d,q,ciphertext,secretkey):
    f,f3 = secretkey
    a = balancedmod(convolution(ciphertext,f,n),q,n)
    return balancedmod(convolution(a,f3,n),3,n)

n = 19
d = ceil(2./3 * n) # number of '1's and '-1's in f and g
q = 512


def attack(n,q,publickey):
    # your code here
