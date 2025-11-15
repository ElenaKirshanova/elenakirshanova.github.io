
# bit-size of p
n = 2048

# p determines the field
p = next_prime(2^(n-1))

# number of MSBs
k = ceil(sqrt(n)) + ceil(log(n, 2))

def msb(k, x):
    """
        Randomised MSBs that returns k bits of x
        p is global
    """
    while True:
        z = randint(1, p-1)
        answer = abs(query - z)
        if answer < p / 2^(k+1):
            break
    return z

def create_oracle(alpha, k):
    """
        The Hidden Number Problem oracle for alpha and k MSBs
    """
    alpha = alpha
    def oracle():
        random_t = randint(1, p-1)
        return random_t, msb(k, (alpha * random_t) % p)
    return oracle


# Example
alpha = randint(1, p-1)
oracle = create_oracle(alpha, k)
