# ------------------------------------------------------------------------------
# In [EG23p], an extension of Regev's factoring quantum algorithm [Regev23] to
# the discrete logarithm problem was introduced. This Sage script (and the
# associated supporting scripts) implements a simulator for the quantum
# algorithm in [EG23p], alongside the classical post-processing algorithm from
# [EG23p] that recovers the logarithm from simulated samples.
#
#   [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
#                           ArXiv 2308.06572v2 (2023).
#
#   [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev’s factoring algorithm
#                                         to compute discrete logarithms".
#                                         ArXiv 2311.05545v2 (2023–2024).

from datetime import datetime
import platform
from dependencies.timer import Timer

load("common.sage")
load("simulator.sage")
load("problem-instances.sage")
load("parameter-search.sage")
load("uids.sage")


def has_maximal_order(g, p):

    """
    @brief  Checks if g has maximal multiplicative order p - 1 mod p, for p a
            prime such that p - 1 is B-smooth for B some small bound.

    @param g  The generator g of order r.

    @param p  The prime p.

    @return   True if g has maximal multiplicative order mod p, False otherwise.
    """

    F = GF(p, proof=False)
    g = F(g)

    for [q, _] in factor(p - 1):
        if g^((p - 1) / q) == 1:
            return False
    return True



def find_smooth_logarithm_mod_p(g, x, p):

    """
    @brief  Finds the discrete logarithm e such that x = g^e mod p, for p a
            prime such that p - 1 is B-smooth for B some small bound.

    The logarithm returned is on [0, r) for r the order of g.

    @param g  The generator g of order r.

    @param x  The element x = g^e mod p.

    @param p  The prime p.

    @return   The discrete logarithm e on [0, r) for r the order of g.
    """

    F = GF(p, proof=False)

    g = F(g)
    x = F(x)
    r = find_smooth_order_mod_p(g, p)

    logs = []

    for [q, _] in factor(p - 1):
        gq = g^((p - 1) / q)
        xq = x^((p - 1) / q)

        eq = discrete_log(xq, gq)

        # Sanity check.
        if gq^eq != xq:
            raise Exception("Error: Internal error (1).")

        logs.append([eq, q])

    # Use the Chinese remainder theorem to compose the solutions.
    e = CRT([d for [d, _] in logs], [q for [_, q] in logs])

    # Sanity check.
    if g^e != x:
        raise Exception("Error: Internal error (2).")

    return e % r


def build_logarithm_finding_element_vector(d, u):

    """
    @brief  Builds the element vector (g_1, .., g_{d-k}, u_1, .., u_k), where
            g_1, .., g_{d-k} are the first d - k primes that are distinct from
            u_1, .., u_k over the integers.

    @param d  The dimension d of the elements vector to build.

    @param u  A list [u_1, .., u_k] of the element u_1, .., u_k.

    @return   The elements vector b = (g_1, .., g_{d-k}, u_1, .., u_k) built.
    """

    k = len(u)

    # Sanity check.
    if k > d:
        raise Exception("Error: The vector u has more than d elements.")

    gis = []

    gi = 2

    while len(gis) < d - k:
        if gi not in u:
            gis.append(gi)

        gi = gi.next_prime()

    return gis + u


def generate_basis_for_logarithm_finding(
    p,
    u,
    d=None,
    *,
    threads=DEFAULT_THREADS,
    verbose=False):

    """
    @brief  Generates a basis for the d-dimensional lattice L_{u_1, .., u_k}
            used for computing orders and discrete logarithms mod p, when
            given the dimension d, the prime p and the k elements u_1, .., u_k.

    The prime p must be such that p - 1 is B-smooth for B some small bound.

    The vectors (z_1, .., z_d) in L_{u_1, .., u_k} are such that

      g_1^{z_1} * ... * g_{d-k}^{z_{d-k}} * \
        u_1^{z_{d-k+1}} * ... * u_k^{z_d} = 1 (mod p),

    where g_1, .., g_{d-k} are the first d - k first primes distinct from the
    elements in u = [u_1, .., u_k].

    @param p  The prime p.

    @param u  A list [u_1, .., u_k] of the elements u_1, .., u_k.

    @param d  The dimension d of the lattice L_{u_1, .., u_k}. Set to
              ceil(sqrt(n)), for n the bit length of p, if omitted.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   A basis for the lattice L_{u_1, .., u_k} used for computing
              discrete logarithms and orders mod p.
    """

    if d == None:
        n = p.nbits()
        d = ceil(sqrt(n))

    b = build_logarithm_finding_element_vector(d, u)
    return generate_basis_for_L(p, [p], b, threads=threads, verbose=verbose)


def solve_samples_for_logarithm(
    samples,
    g,
    x,
    p,
    R,
    *,
    block_size=DEFAULT_BLOCK_SIZE,
    profile_file=None,
    verbose=False):

    """
    @brief  Solves a list of samples for the discrete logarithm e such that
            x = g^e mod p for g, x and p as passed to this function when not
            pre-computing the logarithms of the small generators and instead
            including both g and x in the lattice.

    This function also computes the order r of g from the sampled vectors.

    The prime p must be B-smooth for B some small bound.

    @param samples  The list of samples to solve for the discrete logarithm.

    @param R  The parameter R specifying the standard deviation of the noise.

    @param g  The generator g of order r.

    @param x  The element x = g^e mod p.

    @param p  The prime.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param profile_file   The path to a file in which to save the profile of the
                          Gram–Schmidt norms of the post-processing matrix. If
                          omitted no profile is saved.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The candidate for the discrete logarithm e on [0, r), or None if
              the classical post-processing failed to recover such a candidate.
    """

    # Setup and start a timer.
    timer = Timer().start()

    # Print status.
    if verbose:
        print("Post-processing the sampled vectors to find the logarithm...")
        print(" Building the post-processing matrix...")

    # Extract the dimension.
    d = samples[0].dimensions()[1]
    for i in range(1, len(samples)):
        if d != samples[i].dimensions()[1]:
            raise Exception("Error: The samples differ in their dimensions.")

    # Pick S as described in [Regev23]: In Corollary 4.5 of [Regev23], S = δ^-1,
    # and δ = sqrt(d) / (sqrt(2) R) on p. 7 of [Regev23].
    S = ceil(sqrt(2 / d) * R)
    M = build_M_matrix(samples, S)

    # Run a lattice basis reduction algorithm on the post-processing matrix, and
    # then extract the relevant submatrix.
    if block_size != 2:
        # Use BKZ with the prescribed block size to reduce the lattice basis.
        if verbose:
            print(" Running BKZ on the post-processing matrix...")

        denominator = M.denominator()
        X = (M * denominator).change_ring(ZZ).BKZ(
            block_size=block_size, algorithm="fpLLL", fp="rr", precision=128
        ) / denominator
    else:
        # Use LLL to reduce the lattice basis.
        if verbose:
            print(" Running LLL on the post-processing matrix...")

        X = M.LLL()

    if X.row_space() != M.row_space():
        raise Exception("Error: Failed to run BKZ/LLL.")

    if profile_file != None:
        if verbose:
            print(" Saving the profile after reduction...")

        XR = X.change_ring(RDF)
        (XRgs, _) = XR.gram_schmidt()

        with open(profile_file, "w") as f:
            for i in range(XR.dimensions()[0]):
                f.write(f"{log(abs(XR[i] * XRgs[i]), 2).n()}\n")

    LB = X[:, :d]

    # Build the elements vector b.
    b = build_logarithm_finding_element_vector(d, u = [x, g])

    Lxg = None

    for row in LB:
        if row == 0:
            break

        if not is_in_lattice(row, p, b):
            continue

        if Lxg == None:
            Lxg = matrix(ZZ, row)
            S = Lxg.row_space()
        else:
            if row not in S:
                Lxg = Lxg.stack(row)

    # Sanity check.
    if None == Lxg:
        raise Exception("Error: Failed to find d linearly independent vectors.")

    # Use LLL to remove any linear dependence in the generating set.
    Lxg = Lxg.LLL()[-d:, :]

    if verbose:
        print("  Found", Lxg.rank(), "/", d, "linearly independent vectors...")

    # Sanity checks.
    if Lxg.rank() != d:
        raise Exception("Error: Failed to find d linearly independent vectors.")

    # Solve for v_r = [0, .., 0, r] and v_e = [0, .., 0, a, b].
    S = matrix(ZZ, Lxg)
    Sp = S.hermite_form()

    v_r = Sp[d - 1, :]
    r = v_r[0, d - 1]

    v_d = Sp[d - 2, :]
    a = v_d[0, d - 2]
    b = v_d[0, d - 1]

    # N.B.: There is no point in adding Sp[d - 1, :] to v, since the last
    #       component of Sp[d - 1, :] must be an integer multiple of r.

    tau = gcd(a, r)

    if tau == 1:
        R = IntegerModRing(r)
        e = (R(-b) / R(a)).lift()
    else:
        return None

    if verbose:
        print("")
        print(" Found e =", e)
        print(" Found r =", r)
        print("")
        print(" Time required to post-process:", timer)

    # Return the logarithm.
    return e


def solve_samples_for_logarithm_with_precomputation(
    samples,
    g,
    x,
    p,
    R,
    *,
    block_size=DEFAULT_BLOCK_SIZE,
    profile_file=None,
    verbose=False):

    """
    @brief  Solves a list of samples for the discrete logarithm e such that
            x = g^e mod p for g, x and p as passed to this function when
            pre-computing the logarithms of the small generators.

    The prime p must be such that p - 1 is B-smooth for B some small bound.

    This function assumes that the order r of g mod p is known; it uses the
    fact that p - 1 is smooth to compute the order. Similarly, this function
    assumes that the logarithms e_i such that g = g_i^e_i mod p are known for
    all i in [0, d); it uses the fact p - 1 is smooth to compute the logarithms.

    @param samples  The list of samples to solve for the discrete logarithm.

    @param g  The generator g of order r.

    @param x  The element x = g^e mod p.

    @param p  The prime p.

    @param R  The parameter R specifying the standard deviation of the noise.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param profile_file   The path to a file in which to save the profile of the
                          Gram–Schmidt norms of the post-processing matrix. If
                          omitted no profile is saved.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The candidate for the discrete logarithm e on [0, r), or None if
              the classical post-processing failed to recover such a candidate.
    """

    # Setup and start a timer.
    timer = Timer().start()

    # Print status.
    if verbose:
        print("Post-processing the sampled vectors to find the logarithm...")
        print(" Building the post-processing matrix...")

    # Extract the dimension.
    d = samples[0].dimensions()[1]
    for i in range(1, len(samples)):
        if d != samples[i].dimensions()[1]:
            raise Exception("Error: The samples differ in their dimensions.")

    # Pick S as described in [Regev23]: In Corollary 4.5 of [Regev23], S = δ^-1,
    # and δ = sqrt(d) / (sqrt(2) R) on p. 7 of [Regev23].
    S = ceil(sqrt(2 / d) * R)
    M = build_M_matrix(samples, S)

    # Run a lattice basis reduction algorithm on the post-processing matrix, and
    # then extract the relevant submatrix.
    if block_size != 2:
        # Use BKZ with the prescribed block size to reduce the lattice basis.
        if verbose:
            print(" Running BKZ on the post-processing matrix...")

        denominator = M.denominator()
        X = (M * denominator).change_ring(ZZ).BKZ(
            block_size=block_size, algorithm="fpLLL", fp="rr", precision=128
        ) / denominator
    else:
        # Use LLL to reduce the lattice basis.
        if verbose:
            print(" Running LLL on the post-processing matrix...")

        X = M.LLL()

    if X.row_space() != M.row_space():
        raise Exception("Error: Failed to run BKZ/LLL.")

    if profile_file != None:
        if verbose:
            print(" Saving the profile after reduction...")

        XR = X.change_ring(RDF)
        (XRgs, _) = XR.gram_schmidt()

        with open(profile_file, "w") as f:
            for i in range(XR.dimensions()[0]):
                f.write(f"{log(abs(XR[i] * XRgs[i]), 2).n()}\n")

    LB = X[:d, :d]

    # Build the elements vector b.
    b = build_logarithm_finding_element_vector(d, [x])

    # Assume the order r of g mod p to be known.
    r = find_smooth_order_mod_p(g, p)

    # Pre-compute the logarithms e_i such that g^e_i = b_i mod p.
    e = [find_smooth_logarithm_mod_p(g, bi, p) for bi in b[ : d - 1]]

    # Setup the ring of integers mod p.
    R = IntegerModRing(p)

    if verbose:
        print("")
        print(" Processing the relations found...")

    e_found = None

    for v in LB:
        if not is_in_lattice(v, p, b):
            continue

        # We have:    g_1^v_1 * ... * g_{d-1}^v_{d-1} * x^v_d = 1 = g^0.
        #
        #   => v_1 e_1 + ... + v_{d-1} e_{d-1} + v_d e_d = 0 (mod r)
        #   => e_d = -(v_1 e_1 + ... + v_{d-1} e_{d-1}) / v_d (mod r)

        if gcd(r, v[d - 1]) != 1:
          continue

        Rr = IntegerModRing(r)
        sum_vi_ei = Rr(-sum([v[i] * e[i] for i in range(d - 1)]))
        e_found = (sum_vi_ei / Rr(v[d - 1])).lift()

        if verbose:
          print("  Found e =", e_found)

        break

    if verbose:
        print("")
        print(" Time required to post-process:", timer)

    # Return the logarithm.
    return e_found


def test_logarithm_finding_in_safe_prime_group(
    l=2048,
    *,
    B=DEFAULT_BOUND_B,
    C=DEFAULT_CONSTANT_C,
    dp=1,
    mp=None,
    failure_rate=0,
    block_size=DEFAULT_BLOCK_SIZE,
    threads=DEFAULT_THREADS,
    verbose=True):

    """
    @brief  A convenience function for testing the simulator for the quantum
            algorithm in [EG23p], and the associated classical post-processing,
            with respect to computing discrete logarithms in simulated
            safe-prime groups.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    This function first picks domain parameters for a simulated safe-prime
    group. Specifically, this function first samples an l-bit prime p such that
    p - 1 is B-smooth for some small bound B (to enable efficient simulation).
    It then picks the smallest generator g that is of order r = (p - 1) / 2.

    This function also sets up a problem instance x = g^e mod p by sampling e
    uniformly at random from [0, r) and computing x.

    It then sets up the simulator, and uses the simulator to sample m vectors
    representative of vectors that the quantum algorithm would output for the
    problem instance (g, x, p) according to the analysis in [EG23p].

    Finally, it solves the vectors sampled for the discrete logarithm e by
    using the lattice-based post-processing from [EG23p].

    @param l  The bit length of the prime p.

    @param B  The bound B on the smoothness of p - 1.

    @param C  The constant C that specifies the control register lengths.

    @param dp   A scaling factor d' such that d = ceil(d' * sqrt(n)).

    @param mp   A scaling factor m' such that m = ceil(m' * sqrt(n)). If
                omitted, m = ceil(sqrt(n)) + 4 as in Regev's original analysis.

    @param failure_rate   The failure rate on [0, 1). May be set to a value
                          greater than zero to simulate error-correction
                          failures resulting in bad vectors being output.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   Return True if the discrete logarithm was successfully recovered,
              returns False otherwise.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    # Sample the domain parameters [g, p].
    [g, p] = sample_domain_parameters_safe_prime(l, B, verbose=verbose)

    # Sample [x, e] such that x = g^e where e is uniformly selected from [0, r)
    # for r the order of g.
    [x, e] = sample_x(g, p)

    if verbose:
        print("")
        print(" Sampled e =", e)
        print(" Computed x = g^e mod p =", x)

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    n = p.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    B = generate_basis_for_logarithm_finding(
            p, u=[x, g], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    if verbose:
        print("")
        print("** Sampling vectors...")
        print("")

    R = get_regev_R(C, n)

    samples = simulator.sample_vectors_with_failure_rate(
                  R, m=m, failure_rate=failure_rate, verbose=verbose
              )

    if verbose:
        print("")
        print("** Solving for the logarithm...")
        print("")

    e_found = solve_samples_for_logarithm(
                  samples, g, x, p, R, block_size=block_size, verbose=verbose
              )

    return e_found == e


def test_logarithm_finding_in_schnorr_group(
    l=2048,
    k=224,
    *,
    B=DEFAULT_BOUND_B,
    C=DEFAULT_CONSTANT_C,
    dp=1,
    mp=None,
    failure_rate=0,
    block_size=DEFAULT_BLOCK_SIZE,
    threads=DEFAULT_THREADS,
    verbose=True):

    """
    @brief  A convenience function for testing the simulator for the quantum
            algorithm in [EG23p], and the associated classical post-processing,
            with respect to computing discrete logarithms in simulated Schnorr
            groups.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    This function first picks domain parameters for a simulated Schnorr group.
    Specifically, this function first samples an l-bit prime p such that p - 1
    is B-smooth for some small bound B (to enable efficient simulation), and
    such that p - 1 = 2 * u * r where r is of length approximately k bits. It
    then picks a generator g that is of order r.

    This function also sets up a problem instance x = g^e mod p by sampling e
    uniformly at random from [0, r) and computing x.

    It then sets up the simulator, and uses the simulator to sample m vectors
    representative of vectors that the quantum algorithm would output for the
    problem instance (g, x, p) according to the analysis in [EG23p].

    Finally, it solves the vectors sampled for the discrete logarithm e by
    using the lattice-based post-processing from [EG23p]. This post-processing
    also yields the order r of g, so r need not be assumed known.

    @param l  The bit length of the prime p.

    @param k  The approximate bit length of the order r of g.

    @param B  The bound B on the smoothness of p - 1.

    @param C  The constant C that specifies the control register lengths.

    @param dp   A scaling factor d' such that d = ceil(d' * sqrt(n)).

    @param mp   A scaling factor m' such that m = ceil(m' * sqrt(n)). If
                omitted, m = ceil(sqrt(n)) + 4 as in Regev's original analysis.

    @param failure_rate   The failure rate on [0, 1). May be set to a value
                          greater than zero to simulate error-correction
                          failures resulting in bad vectors being output.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   Return True if the discrete logarithm was successfully recovered,
              returns False otherwise.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    # Sample the domain parameters [g, p].
    [g, p] = sample_domain_parameters_schnorr(l, k, B, verbose=verbose)

    # Sample [x, e] such that x = g^e where e is uniformly selected from [0, r)
    # for r the order of g.
    [x, e] = sample_x(g, p)

    if verbose:
        print("")
        print(" Sampled e =", e)
        print(" Computed x = g^e mod p =", x)

    # 1. Solve for e = log_{g}(x).

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    n = p.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    B = generate_basis_for_logarithm_finding(
            p, u=[x, g], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    if verbose:
        print("")
        print("** Sampling vectors...")
        print("")

    R = get_regev_R(C, n)

    samples = simulator.sample_vectors_with_failure_rate(
                  R, m=m, failure_rate=failure_rate, verbose=verbose
              )

    if verbose:
        print("")
        print("** Solving for the logarithm...")
        print("")

    e_found = solve_samples_for_logarithm(
                  samples, g, x, p, R, block_size=block_size, verbose=verbose
              )

    return e_found == e


def test_logarithm_finding_in_safe_prime_group_with_precomputation(
    l=2048,
    *,
    B=DEFAULT_BOUND_B,
    C=DEFAULT_CONSTANT_C,
    dp=1,
    mp=None,
    failure_rate=0,
    block_size=DEFAULT_BLOCK_SIZE,
    threads=DEFAULT_THREADS,
    verbose=True):

    """
    @brief  A convenience function for testing the simulator for the quantum
            algorithm in [EG23p], and the associated classical post-processing,
            with respect to computing discrete logarithms in simulated
            safe-prime groups with pre-computation.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    This function first picks domain parameters for a simulated safe-prime
    group. Specifically, this function first samples an l-bit prime p such that
    p - 1 is B-smooth for some small bound B (to enable efficient simulation).
    It then picks a generator g that is of order r = (p - 1) / 2.

    This function also sets up a problem instance x = g^e mod p by sampling e
    uniformly at random from [0, r) and computing x.

    It then sets up the simulator, and uses the simulator to sample vectors
    representative of vectors that the quantum algorithm would output for the
    problem instance (g, x, p) according to the analysis in [EG23p].

    Specifically, this function sets up and uses the simulator two times:

    - It first picks a generator g_max of the full group Z_p^*.

    - Using the assumed knowledge of e_i such that g_i = g_max^e_i, it solves m
      vectors sampled for the problem instance (g_max, g, p) for the discrete
      logarithm e_g such that g = g_max^e_g using the lattice-based
      post-processing from [EG23p]. This step may be pre-computed.

    - Again, using the assumed knowledge of e_i such that g_i = g_max^e_i, it
      solves m vectors sampled for the problem instance (g_max, x, p) for the
      discrete logarithm e_x such that x = g_max^e_x using the lattice-based
      post-processing from [EG23p].

    Finally, using that r = (p - 1) / 2, it computes e = e_x / e_g mod r.

    @param l  The bit length of the prime p.

    @param B  The bound B on the smoothness of p - 1.

    @param C  The constant C that specifies the control register lengths.

    @param dp   A scaling factor d' such that d = ceil(d' * sqrt(n)).

    @param mp   A scaling factor m' such that m = ceil(m' * sqrt(n)). If
                omitted, m = ceil(sqrt(n)) + 4 as in Regev's original analysis.

    @param failure_rate   The failure rate on [0, 1). May be set to a value
                          greater than zero to simulate error-correction
                          failures resulting in bad vectors being output.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   Return True if the discrete logarithm was successfully recovered,
              returns False otherwise.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    # Sample the domain parameters [g, p].
    [g, p] = sample_domain_parameters_safe_prime(l, B=B, verbose=verbose)

    # Sample [x, e] such that x = g^e where e is uniformly selected from [0, r)
    # for r the order of g.
    [x, e] = sample_x(g, p)

    if verbose:
        print("")
        print(" Sampled e =", e)
        print(" Computed x = g^e mod p =", x)

    # Sample a generator g_max for the full group.
    g_max = 2
    while not has_maximal_order(g_max, p):
        g_max += 1

    if verbose:
        print("")
        print(" Selected g_max =", g_max)

    # 1. Solve for e_x = log_{g_max}(x).

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator to compute e = log_{g_max}(x) ...")
        print("")

    n = p.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    B = generate_basis_for_logarithm_finding(
            p, u=[x], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    # Use the simulator to sample vectors.
    if verbose:
        print("")
        print("** Using the simulator to sample vectors...")
        print("")

    R = get_regev_R(C, n)

    samples = simulator.sample_vectors_with_failure_rate(
                  R, m=m, failure_rate=failure_rate, verbose=verbose
              )

    # Solve the vectors sampled for the logarithm.
    if verbose:
        print("")
        print("** Post-processing the sampled vectors to find " + \
                "e = log_{g_max}(x)...")
        print("")

    e_x = solve_samples_for_logarithm_with_precomputation(
              samples,
              g_max,
              x,
              p,
              R,
              block_size=block_size,
              verbose=verbose)

    # 2. Solve for e_g = log_{g_max}(g).

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator for e = log_{g_max}(g)...")
        print("")

    B = generate_basis_for_logarithm_finding(
            p, u=[g], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    # Use the simulator to sample vectors.
    if verbose:
        print("")
        print("** Using the simulator to sample vectors...")
        print("")

    samples = simulator.sample_vectors_with_failure_rate(
                  R, m=m, failure_rate=failure_rate, verbose=verbose
              )

    # Solve the vectors sampled for the logarithm.
    if verbose:
        print("")
        print("** Post-processing the sampled vectors to find " + \
                "e = log_{g_max}(g)...")
        print("")

    e_g = solve_samples_for_logarithm_with_precomputation(
              samples,
              g_max,
              g,
              p,
              R,
              block_size=block_size,
              verbose=verbose)

    # 3. Solve for the logarithm e = e_x / e_g (mod r) where r is assumed known.

    if verbose:
        print("")
        print("** Combining the results to solve for e = log_{g}(x)...")
        print("")

    r = find_smooth_order_mod_p(g, p)
    R = IntegerModRing(r)

    e_found = R(e_x) / R(e_g)

    if verbose:
      print("Found e =", e_found)

    return e_found == e


def test_logarithm_finding_in_schnorr_group_with_precomputation(
    l=2048,
    k=224,
    *,
    B=DEFAULT_BOUND_B,
    C=DEFAULT_CONSTANT_C,
    dp=1,
    mp=None,
    failure_rate=0,
    block_size=DEFAULT_BLOCK_SIZE,
    threads=DEFAULT_THREADS,
    verbose=True):

    """
    @brief  A convenience function for testing the simulator for the quantum
            algorithm in [EG23p], and the associated classical post-processing,
            with respect to computing discrete logarithms in simulated
            Schnorr groups with pre-computation.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    This function first picks domain parameters for a simulated Schnorr group.
    Specifically, this function first samples an l-bit prime p such that p - 1
    is B-smooth for some small bound B (to enable efficient simulation), and
    such that p - 1 = 2 * u * r where r is of length approximately k bits. It
    then picks a generator g that is of order r.

    This function also sets up a problem instance x = g^e mod p by sampling e
    uniformly at random from [0, r) and computing x.

    It then sets up the simulator, and uses the simulator to sample vectors
    representative of vectors that the quantum algorithm would output for the
    problem instance (g, x, p) according to the analysis in [EG23p].

    Specifically, this function sets up and uses the simulator two times:

    - It first picks a generator g_max of the full group Z_p^*.

    - Using the assumed knowledge of e_i such that g_i = g_max^e_i, it solves m
      vectors sampled for the problem instance (g_max, g, p) for the discrete
      logarithm e_g such that g = g_max^e_g using the lattice-based
      post-processing from [EG23p]. This step may be pre-computed.

    - Again, using the assumed knowledge of e_i such that g_i = g_max^e_i, it
      solves m vectors sampled for the problem instance (g_max, x, p) for the
      discrete logarithm e_x such that x = g_max^e_x using the lattice-based
      post-processing from [EG23p].

    Finally, using the assumed knowlege of r, it computes e = e_x / e_g mod r.

    @param l  The bit length of the prime p.

    @param k  The approximate bit length of the order r of g.

    @param B  The bound B on the smoothness of p - 1.

    @param C  The constant C that specifies the control register lengths.

    @param dp   A scaling factor d' such that d = ceil(d' * sqrt(n)).

    @param mp   A scaling factor m' such that m = ceil(m' * sqrt(n)). If
                omitted, m = ceil(sqrt(n)) + 4 as in Regev's original analysis.

    @param failure_rate   The failure rate on [0, 1). May be set to a value
                          greater than zero to simulate error-correction
                          failures resulting in bad vectors being output.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   Return True if the discrete logarithm was successfully recovered,
              returns False otherwise.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    # Sample the domain parameters [g, p].
    [g, p] = sample_domain_parameters_schnorr(l, k, B=B, verbose=verbose)

    # Sample [x, e] such that x = g^e where e is uniformly selected from [0, r)
    # for r the order of g.
    [x, e] = sample_x(g, p)

    if verbose:
        print("")
        print(" Sampled e =", e)
        print(" Computed x = g^e mod p =", x)

    # Sample a generator g_max for the full group.
    g_max = 2
    while not has_maximal_order(g_max, p):
        g_max += 1

    if verbose:
        print("")
        print(" Selected g_max =", g_max)

    # 1. Solve for e_x = log_{g_max}(x).

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator to compute e = log_{g_max}(x) ...")
        print("")

    n = p.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    B = generate_basis_for_logarithm_finding(
            p, u=[x], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    # Use the simulator to sample vectors.
    if verbose:
        print("")
        print("** Using the simulator to sample vectors...")
        print("")

    R = get_regev_R(C, n)

    samples = simulator.sample_vectors_with_failure_rate(
                  R, m=m, failure_rate=failure_rate, verbose=verbose
              )

    # Solve the vectors sampled for the logarithm.
    if verbose:
        print("")
        print("** Post-processing the sampled vectors to find " + \
                "e = log_{g_max}(x)...")
        print("")

    e_x = solve_samples_for_logarithm_with_precomputation(
              samples,
              g_max,
              x,
              p,
              R,
              block_size=block_size,
              verbose=verbose)

    # 2. Solve for e_g = log_{g_max}(g).

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator for e = log_{g_max}(g)...")
        print("")

    B = generate_basis_for_logarithm_finding(
            p, u=[g], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    # Use the simulator to sample vectors.
    if verbose:
        print("")
        print("** Using the simulator to sample vectors...")
        print("")

    samples = simulator.sample_vectors_with_failure_rate(
                  R, m=m, failure_rate=failure_rate, verbose=verbose
              )

    # Solve the vectors sampled for the logarithm.
    if verbose:
        print("")
        print("** Post-processing the sampled vectors to find " + \
                "e = log_{g_max}(g)...")
        print("")

    e_g = solve_samples_for_logarithm_with_precomputation(
              samples,
              g_max,
              g,
              p,
              R,
              block_size=block_size,
              verbose=verbose)

    # 3. Solve for the logarithm e = e_x / e_g (mod r) where r is assumed known.

    if verbose:
        print("")
        print("** Combining the results to solve for e = log_{g}(x)...")
        print("")

    r = find_smooth_order_mod_p(g, p)
    R = IntegerModRing(r)

    e_found = R(e_x) / R(e_g)

    if verbose:
      print("Found e =", e_found)

    return e_found == e


def find_minimum_C_for_logarithm_finding(
    p,
    g,
    x,
    *,
    precision=DEFAULT_PRECISION,
    dp=1,
    mp=None,
    failure_rate=0,
    probabilistic_failures=False,
    log_file_prefix=None,
    block_size=DEFAULT_BLOCK_SIZE,
    threads=DEFAULT_THREADS,
    verbose=True):

    """
    @brief  Finds the minimum value of the constant C, with a given prescribed
            precision, such that the post-processing succeeds in solving a given
            discrete logarithm problem instance.

    Given a prime p, a generator g of unknown order r mod p, and x = g^e mod p,
    this function finds the minimum value of the constant C for which the
    discrete logarithm e can be recovered.

    Specifically, it performs a binary search to find the minimum value of C
    that allows the discrete logarithm e to be recovered. This when setting up a
    simulator, using it to sample m vectors, and post-processing the samples,
    for each value of C to test.

    Note that logarithm finding is attempted for the given problem instance
    only for each value of C tested. Additional tests may hence need to be
    performed for the minimum value of C returned by this function, and for
    slightly smaller and larger values of C.

    @param p  The prime p such that p - 1 is B-smooth for B some small bound.

    @param g  The generator g of order r.

    @param x  The element x = g^e mod p.

    @param precision  The level of precision with which to perform the search.

    @param dp   A scaling factor d' such that d = ceil(d' * sqrt(n)).

    @param mp   A scaling factor m' such that m = ceil(m' * sqrt(n)). If
                omitted, m = ceil(sqrt(n)) + 4 as in Regev's original analysis.

    @param failure_rate   The failure rate on [0, 1). May be set to a value
                          greater than zero to simulate error-correction
                          failures resulting in bad vectors being output.

    @param probabilistic_failures   A flag that shall be set to True if bad
                                    vectors shall be sampled probabilistically
                                    with probability, as indicated by the
                                    failure rate, and to False if a fixed
                                    portion of the runs shall be bad again as
                                    indicated by the failure rate.

    @param log_file_prefix  The prefix to use when constructing the path to log
                            files in which to store partial results during the
                            search. If omitted no log files are saved.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The minimum value of the constant C.
    """

    e_expected = find_smooth_logarithm_mod_p(g, x, p)

    n = p.nbits()
    d = ceil(sqrt(n) * dp)

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    B = generate_basis_for_logarithm_finding(
            p, u=[x, g], d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    log_file = None
    if None != log_file_prefix:
        # Append the run ID to the profile file prefix.
        log_file_prefix = log_file_prefix + uid_for_N(p) + "-" + uid() + "-"

        log_file = open(log_file_prefix + "log.txt", "w")

        print(
            "Start time:",
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            file=log_file,
        )
        print("Hostname:", platform.node(), "\n", file=log_file)

        print("p =", p, "\n", file=log_file)

        print("g =", g, file=log_file)
        print("x =", x, "\n", file=log_file)

        print("d =", simulator.d, file=log_file)
        print("m =", m, "\n", file=log_file)

        print("block_size =", block_size, "\n", file=log_file)

        print("failure_rate =", failure_rate, file=log_file)
        print("probabilistic_failures =", probabilistic_failures, "\n", \
                file=log_file)

        log_file.flush()

    # Define a function that tests if the post-processing is successful for a
    # given value of the constant C. Passed to min_val_for_C() below.
    def test_for_C(C):
        R = get_regev_R(C, n)

        if probabilistic_failures:
            samples = simulator.sample_vectors_with_failure_rate(
                R, m=m, failure_rate=failure_rate, verbose=verbose
            )
        else:
            m1 = ceil(m * (1 - failure_rate))
            m2 = floor(m * failure_rate)

            samples = simulator.sample_vectors_good_and_bad(
                R, m1, m2, verbose=verbose
            )

        if verbose:
            print("")

        profile_file = None
        if log_file_prefix != None:
            profile_file = \
              f"{log_file_prefix}profile-{simulator.d}-{m}-{C}.dat"

        try:
            e_found = solve_samples_for_logarithm(
                          samples,
                          g,
                          x,
                          p,
                          R,
                          block_size=block_size,
                          profile_file=profile_file,
                          verbose=verbose)

            # Test the logarithm.
            solved = (e_found == e_expected)
        except Exception as e:
            if verbose:
                print(f"Error: Finding e failed with exception: {e}")

            solved = False

        if log_file_prefix != None:
            # Write out the metadata file.
            meta_file = f"{log_file_prefix}meta-{simulator.d}-{m}-{C}.txt"

            with open(meta_file, "w") as f:
                print("C =", C, file=f)
                print("R =", R, "\n", file=f)

                print("success =", solved, file=f)

            # Sample the samples.
            samples_file = \
              f"{log_file_prefix}samples-{simulator.d}-{m}-{C}.sobj"
            save(samples, samples_file)

        return solved

    return find_minimum_C(test_for_C, precision, log_file, verbose=verbose)


def find_minimum_C_for_logarithm_finding_in_safe_prime_group(
    l=2048,
    *,
    precision=DEFAULT_PRECISION,
    B=DEFAULT_BOUND_B,
    dp=1,
    mp=None,
    failure_rate=0,
    probabilistic_failures=False,
    log_file_prefix=None,
    block_size=DEFAULT_BLOCK_SIZE,
    threads=DEFAULT_THREADS,
    verbose=True):

    """
    @brief  Finds the minimum value, with given precision, such that the
            post-processing successfully solves a random discrete logarithm
            problem instance in a randomly selected simulated safe-prime group.

    This function first picks domain parameters for a simulated safe-prime
    group. Specifically, this function first samples an l-bit prime p such that
    p - 1 is B-smooth for some small bound B (to enable efficient simulation).
    It then picks a generator g that is of order r = (p - 1) / 2.

    This function also sets up a problem instance x = g^e mod p by sampling e
    uniformly at random from [0, r) and computing x.

    It then performs a binary search to find the minimum value of C that allows
    the discrete logarithm e to be recovered. This when setting up a simulator,
    using it to sample m vectors, and post-processing the samples, for each
    value of C to test.

    Note that logarithm finding is attempted for the generated problem instance
    only for each value of C tested. Additional tests may hence need to be
    performed for the minimum value of C returned by this function, and for
    slightly smaller and larger values of C.

    Note furthermore that this is a convenience function: It generates domain
    parameters for a simulated safe-prime group, and a random problem instance.
    It then simply calls the find_minimum_C_for_logarithm_finding() function.

    @param l  The bit length l of the prime p.

    @param precision  The level of precision with which to perform the search.

    @param B  The bound B on the smoothness of p - 1.

    @param dp   A scaling factor d' such that d = ceil(d' * sqrt(n)).

    @param mp   A scaling factor m' such that m = ceil(m' * sqrt(n)). If
                omitted, m = ceil(sqrt(n)) + 4 as in Regev's original analysis.

    @param failure_rate   The failure rate on [0, 1). May be set to a value
                          greater than zero to simulate error-correction
                          failures resulting in bad vectors being output.

    @param probabilistic_failures   A flag that shall be set to True if bad
                                    vectors shall be sampled probabilistically
                                    with probability, as indicated by the
                                    failure rate, and to False if a fixed
                                    portion of the runs shall be bad again as
                                    indicated by the failure rate.

    @param log_file_prefix  The prefix to use when constructing the path to log
                            files in which to store partial results during the
                            search. If omitted no log files are saved.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The minimum value of the constant C.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    [g, p] = sample_domain_parameters_safe_prime(l, B=B, verbose=verbose)

    [x, e] = sample_x(g, p)

    if verbose:
        print("")
        print(" Sampled e =", e)
        print(" Computed x = g^e mod p =", x)

    result = find_minimum_C_for_logarithm_finding(
        p,
        g,
        x,
        precision=precision,
        dp=dp,
        mp=mp,
        failure_rate=failure_rate,
        probabilistic_failures=probabilistic_failures,
        log_file_prefix=log_file_prefix,
        block_size=block_size,
        threads=threads,
        verbose=verbose)

    return result
