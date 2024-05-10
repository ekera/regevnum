from dependencies.timer import Timer


def sample_smooth_prime_excluding_factors(l, B=DEFAULT_BOUND_B, factors=set()):

    """
    @brief  Samples an l-bit prime p such that p - 1 is B-smooth for B some
            small bound, and such that p - 1 has no factor that is in the
            factors set besides two.

    @param l  The bit length of each prime p.

    @param B  The bound B on the smoothness of p - 1.

    @param factors  The set of prime factors less than B already used. If two
                    is a part of this set, two will be ignored, since p - 1
                    must be even when p is a large prime.

    @return   [p, factors], for p the prime selected, and factors an updated
              set of used factors that also contains the factors of p - 1.
    """

    # Create a list F of all prime factors up to B not in factors.
    F = []
    for p in primes(3, B):
        if p not in factors:
            F.append(p)

    # Search exhaustively for a combination of factors that yields an l-bit
    # prime p such that p - 1 is B-smooth and return this prime.
    while True:
        shuffle(F)
        used = []

        x = 2

        for p in F:
            while x * p >= 2^l:
                x = x // used[0]
                used = used[1:]

            x *= p
            used.append(p)

            if (2^(l - 1) <= x) and (x < 2^l):
                q = x + 1

                if q.is_prime(proof=False):
                    return [q, factors.union(set(used)).union(set([2]))]


def sample_smooth_prime(l, B=DEFAULT_BOUND_B):

    """
    @brief  Samples an l-bit prime p such that p - 1 is B-smooth for B some
            small bound.

    @param l  The bit length of the prime p.

    @param B  The bound B on the smoothness of p - 1.

    @return   The prime p sampled.
    """

    return sample_smooth_prime_excluding_factors(l, B)[0]


def sample_integer(t, l, B=DEFAULT_BOUND_B, verbose=False):

    """
    @brief  Samples an integer N = p_1 * ... * p_t, with t distinct l-bit prime
            factors such that all p_i - 1 are B-smooth for B some small bound,
            and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    @param t  The number of distinct l-bit primes p_i.

    @param l  The bit length of each prime p_i.

    @param B  The bound B on the smoothness of each p_i - 1.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   [N, [p_1, .., p_t]] for N = p_1 * ... * p_t.
    """

    # Setup and start a timer.
    timer = Timer().start()

    if verbose:
        print("Sampling N on special form to enable efficient simulation...")

    # Primes.
    primes = set()
    factors = set()

    # Pick t distinct l-bit primes p_i, such that p_i - 1 are B-smooth.
    for _ in range(t):
        while True:
            [p, factors] = sample_smooth_prime_excluding_factors(l, B, factors)

            if p not in primes:
                primes.add(p)

                if verbose:
                    print(" Sampled factor:", p)

                break

    # Compute the product.
    N = prod(primes)

    # Stop the timer.
    timer.stop()

    if verbose:
        print("")
        print(" Sampled N =", N)
        print("")
        print(" Time required to sample:", timer)

    return [N, primes]


def find_smooth_order_mod_p(g, p):

    """
    @brief  Finds the multiplicative order of g mod p, for p a prime such that
            p - 1 is B-smooth for some small B.

    @param g  The generator g.

    @param p  The prime p such that p - 1 is B-smooth.

    @return   The multiplicative order of g mod p.
    """

    r = p - 1

    F = GF(p, proof=False)
    g = F(g)

    for [q, _] in factor(p - 1):
        if g^((p - 1) / q) == 1:
            r /= q

    return ZZ(r)


def sample_x(g, p):

    """
    @brief  Samples an element x = g^e mod p uniformly at random from <g> by
            sampling e uniformly at random from [0, r), for r the order of g,
            and for p a prime such that p - 1 is B-smooth for B some small
            bound.

    @param g  The generator g.

    @param p  The prime p.

    @return   [x, e] for x = g^e the element sampled and e the exponent.
    """

    R = IntegerModRing(p)

    r = find_smooth_order_mod_p(g, p)

    e = IntegerModRing(r).random_element().lift()
    x = (R(g)^e).lift()

    return [x, e]


def sample_safe_prime(l):

    """
    @brief  Samples an l-bit prime p such that (p - 1) / 2 is also a prime.

    @param l  The bit length of the prime p.

    @return   The prime p selected.
    """

    while True:
        p = random_prime(2^l-1, proof=False, lbound=2^(l-1))
        if ZZ((p - 1)/2).is_prime():
            return p


def sample_domain_parameters_safe_prime(l, B=DEFAULT_BOUND_B, verbose=False):

    """
    @brief  Samples domain parameters to simulate a safe-prime group.

    Specifically, this function first samples an l-bit prime p such that p - 1
    is B-smooth for B some small bound (to enable efficient simulation). It then
    picks the first generator g that is of order r = (p - 1) / 2.

    The above assuming that ACTUAL_SAFE_PRIMES is set to False in common.sage,
    otherwise this function samples an actual safe-prime; an l-bit prime p such
    that (p - 1) / 2 is also a prime.

    @param l  The bit length of the prime p.

    @param B  The bound B on the smoothness of p - 1.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   [g, p] for g the generator picked and p the prime sampled.
    """

    # Setup and start a timer.
    timer = Timer().start()

    if verbose:
        print("Sampling domain parameters...")

    if ACTUAL_SAFE_PRIMES:
        p = sample_safe_prime(l)
    else:
        p = sample_smooth_prime(l, B)

    if verbose:
        print("")
        print(" Sampled p =", p)

    g = 2

    while find_smooth_order_mod_p(g, p) != (p - 1) / 2:
        g += 1

        if g >= p:
            raise Exception("Error: Failed to sample g.")

    if verbose:
        print(" Sampled g =", g)
        print("")
        print(" Time required to sample:", timer)

    return [g, p]


def sample_domain_parameters_schnorr(l, k, B=DEFAULT_BOUND_B, verbose=False):

    """
    @brief  Samples domain parameters to simulate a Schnorr group.

    Specifically, this function first samples an l-bit prime p such that p - 1
    is B-smooth for B some small bound (to enable efficient simulation), and
    such that p - 1 = 2 * u * r where r is of length approximately k bits.

    It then picks a random generator g that is of order r.

    @param l  The bit length of the prime p.

    @param k  The approximate bit length of the order r of g.

    @param B  The bound B on the smoothness of p - 1.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   [g, p] for g the generator picked and p the prime sampled.
    """

    # Setup and start a timer.
    timer = Timer().start()

    if verbose:
        print("Sampling domain parameters...")

    p = sample_smooth_prime(l, B)

    if verbose:
        print(" Sampled p =", p)

    # Select r.
    factors = [q for [q, _] in factor((p - 1) / 2)]
    shuffle(factors)

    r = 1
    for q in factors:
        if r.nbits() < k:
            r *= q

    # Select g.
    while True:
        g = IntegerModRing(p).random_element()
        if g == 0:
            continue;

        g = g^((p - 1) / r)
        g = g.lift()

        if find_smooth_order_mod_p(g, p) == r:
            break

    if verbose:
        print(" Sampled g =", g)
        print("")
        print(" Time required to sample:", timer)

    return [g, p]
