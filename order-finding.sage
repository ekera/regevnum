# ------------------------------------------------------------------------------
# In [EG23p], an extension of Regev's factoring quantum algorithm [Regev23] to
# the order-finding problem was introduced. This Sage script (and the
# associated supporting scripts) implements a simulator for the quantum
# algorithm in [EG23p], alongside the classical post-processing algorithm from
# [EG23p] that recovers the order from simulated samples.
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

load("dependencies/factor.sage")


def find_smooth_order_mod_N(g, N, F):

    """
    @brief  Finds the multiplicative order of g mod N, for N = p_1 * ... * p_t a
            composite integer with t distinct prime factors such that p_i - 1
            is B-smooth for B some small bound.

    @param g  The generator g.

    @param N  The integer N.

    @param F  A list [p_1, .., p_t] of the factors of N.

    @return   The multiplicative order of g mod N.
    """

    R = IntegerModRing(N)

    g = R(g)

    r = lcm([p - 1 for p in F]) # = phi(N)

    if g^r != 1:
        raise Exception("Error: Unexpected order.")

    for [f, _] in factor(r):
        while True:
          if g^(r / f) == 1:
              r /= f
          else:
              break

          if 0 != (r % f):
              break

    return r


def build_order_finding_element_vector(d, u):

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


def generate_basis_for_order_finding(
    N,
    F,
    u=[],
    d=None,
    *,
    threads=DEFAULT_THREADS,
    verbose=False):

    """
    @brief  Generates a basis for the d-dimensional lattice L_{u_1, .., u_k}
            used for computing orders mod N, when given the dimension d, the
            integer N to factor, the factors of N and the elements u_1, .., u_k.

    The vectors (z_1, .., z_d) in L_{u_1, .., u_k} are such that

      g_1^{z_1} * ... * g_{d-k}^{z_{d-k}} * \
        u_1^{z_{d-k+1}} * ... * u_k^{z_d} = 1 (mod N),

    where g_1, .., g_{d-k} are the first d - k first primes distinct from the
    elements in u = [u_1, .., u_k], and where N = p_1 * ... * p_t is a composite
    integer with t distinct prime factors such that all p_i - 1 are B-smooth for
    B some small bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    @param N  The integer N to factor.

    @param F  A list [p_1, .., p_t] of the factors of N.

    @param u  A list [u_1, .., u_k] of the elements u_1, .., u_k.

    @param d  The dimension d of the lattice L_{u_1, .., u_k}. Set to
              ceil(sqrt(n)), for n the bit length of p, if omitted.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   A basis for the lattice L_{u_1, .., u_k} used for computing orders
              mod N.
    """

    if None == d:
        n = N.nbits()
        d = ceil(sqrt(n))

    b = build_order_finding_element_vector(d, u)
    return generate_basis_for_L(N, F, b, threads=threads, verbose=verbose)


def solve_samples_for_order(
    samples,
    g,
    N,
    R,
    *,
    block_size=DEFAULT_BLOCK_SIZE,
    profile_file=None,
    verbose=False):

    """
    @brief  Solves g for the order r, the least positive integer r such that
            g^r = 1 mod p, given a list of samples.

    @param samples  The list of samples to solve for the order.

    @param R  The parameter R specifying the standard deviation of the noise.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param profile_file   The path to a file in which to save the profile of the
                          Gram–Schmidt norms of the post-processing matrix. If
                          omitted no profile is saved.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   A candidate for the order r, or None if the classical
              post-processing failed to recover such a candidate.
    """

    # Setup and start a timer.
    timer = Timer().start()

    # Print status.
    if verbose:
        print("Post-processing the sampled vectors to find the order...")
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
        raise Exception("Error: Failed to run BKZ / LLL.")

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
    b = build_order_finding_element_vector(d, u=[g])

    Lg = None

    for row in LB:
        if row == 0:
            break

        if not is_in_lattice(row, N, b):
            continue

        if Lg == None:
            Lg = matrix(ZZ, row)
            S = Lg.row_space()
        else:
            if row not in S:
                Lg = Lg.stack(row)

    # Sanity check.
    if None == Lg:
        raise Exception("Error: Failed to find d linearly independent vectors.")

    # Ensure that Lg is of the right type.
    Lg = matrix(ZZ, Lg)

    # Use LLL to remove any linear dependence in the generating set.
    Lg = Lg.LLL()[-d:, :]

    if verbose:
        print("  Found", Lg.rank(), "/", d, "linearly independent vectors...")

    # Sanity check.
    if Lg.rank() != d:
        raise Exception("Error: Failed to find d linearly independent vectors.")

    # Solve for v = [0, .., 0, 0, r].
    v = Lg.hermite_form()[d - 1, :]

    r = v[0, d - 1]

    if verbose:
        print("")
        print(" Found r =", r)

    if verbose:
        print("")
        print(" Time required to post-process:", timer)

    return r


def solve_samples_for_phi(
    samples,
    N,
    R,
    *,
    block_size=DEFAULT_BLOCK_SIZE,
    profile_file=None,
    verbose=False):

    """
    @brief  Solve a list of samples for Euler's phi function.

    @param samples  The list of samples to solve for Euler's phi function.

    @param R  The parameter R specifying the standard deviation of the noise.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param profile_file   The path to a file in which to save the profile of the
                          Gram–Schmidt norms of the post-processing matrix. If
                          omitted no profile is saved.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   A candidate for Euler's phi function, or None if the classical
              post-processing failed to recover such a candidate.
    """

    # Setup and start a timer.
    timer = Timer().start()

    # Print status.
    if verbose:
        print("Post-processing the sampled vectors to find phi...")
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
        raise Exception("Error: Failed to run BKZ / LLL.")

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
    b = build_order_finding_element_vector(d, u=[])

    L = None

    for row in LB:
        if row == 0:
            break

        if not is_in_lattice(row, N, b):
            continue

        if L == None:
            L = matrix(ZZ, row)
            S = L.row_space()
        else:
            if row not in S:
                L = L.stack(row)

    # Sanity check.
    if None == L:
        raise Exception("Error: Failed to find d linearly independent vectors.")

    # Ensure that L is of the right type.
    L = matrix(ZZ, L)

    # Use LLL to remove any linear dependence in the generating set.
    L = L.LLL()[-d:, :]

    if verbose:
        print("  Found", L.rank(), "/", d, "linearly independent vectors...")

    # Sanity check.
    if L.rank() != d:
        raise Exception("Error: Failed to find d linearly independent vectors.")

    phi = abs(L.determinant())

    if verbose:
        print("")
        print(" Found phi =", phi)
        print("")
        print(" Time required to post-process:", timer)

    return phi


def test_order_finding(
    t=2,
    l=1024,
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
            with respect to computing the order of elements in the
            multiplicative group of the ring of integers mod N.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    This function first picks an integer N = p_1 * ... * p_t, with t distinct
    l-bit prime factors such that all p_i - 1 are B-smooth for B some small
    bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    It furthermore picks an element g uniformly at random from the
    multiplicative group of the ring of integers mod N.

    It then sets up the simulator, and uses the simulator to sample m vectors
    representative of vectors that the quantum computer would output according
    to Regev's analysis of the quantum algorithm in [Regev23], and the extended
    analysis of Ekerå and Gärtner in [EG23p], when performing order finding.

    Finally, it solves the vectors sampled for the order r of g mod N, the least
    positive integer such that g^r = 1 mod N.

    @param t  The number of distinct prime factors of N.

    @param l  The bit length l of each distinct prime factor of m.

    @param B  The bound B on the smoothness of p_i - 1.

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

    @return   Return True if the order was successfully recovered, returns False
              otherwise.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    [N, F] = sample_integer(t, l, B=B, verbose=verbose)

    if verbose:
        print("")
        print("Sampling a generator...")

    R = IntegerModRing(N)
    while True:
        g = R.random_element()
        if gcd(g.lift(), N) == 1:
            break;

    if verbose:
        print(" Sampled g =", g)

    n = N.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    B = generate_basis_for_order_finding(
            N, F, u=[g], d=d, threads=threads, verbose=verbose
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
        print("** Solving for the order...")
        print("")

    order_found = solve_samples_for_order(
                      samples, g, N, R, block_size=block_size, verbose=verbose
                  )

    order_expected = find_smooth_order_mod_N(g, N, F)

    return order_found == order_expected


def test_order_finding_phi(
    t=2,
    l=1024,
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
            with respect to computing the order phi of the multiplicative group
            of the ring of integers mod N.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                            ArXiv 2308.06572v2 (2023).

    This function first picks an integer N = p_1 * ... * p_t, with t distinct
    l-bit prime factors such that all p_i - 1 are B-smooth for B some small
    bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    It then sets up the simulator, and uses the simulator to sample m vectors
    representative of vectors that the quantum computer would output according
    to Regev's analysis of the quantum algorithm in [Regev23], and the extended
    analysis of Ekerå and Gärtner in [EG23p], when performing order finding.

    Finally, it solves the vectors sampled for the order phi of the
    multiplicative group of the ring of integers mod N.

    @param t  The number of distinct prime factors of N.

    @param l  The bit length l of each distinct prime factor of m.

    @param B  The bound B on the smoothness of p_i - 1.

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

    @return   Return True if the order was successfully recovered, returns False
              otherwise.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    [N, F] = sample_integer(t, l, B=B, verbose=verbose)

    n = N.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    # Setup the simulator.
    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    B = generate_basis_for_order_finding(
            N, F, u=[], d=d, threads=threads, verbose=verbose
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
        print("** Solving for phi...")
        print("")

    phi_found = solve_samples_for_phi(
                    samples, N, R, block_size=block_size, verbose=verbose
                )

    phi_expected = prod([p - 1 for p in F])

    return phi_found == phi_expected


def find_minimum_C_for_order_finding(
    t=2,
    l=1024,
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
    @brief  Finds the minimum value of the constant C, with a given prescribed
            precision, such that the post-processing succeeds in recovering the
            order of an element selected uniformly at random from the
            multiplicative group of the ring of integers mod N.

    This function first picks an integer N = p_1 * ... * p_t, with t distinct
    l-bit prime factors such that all p_i - 1 are B-smooth for B some small
    bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    It furthermore picks an element g uniformly at random from the
    multiplicative group of the ring of integers mod N.

    It then performs a binary search to find the minimum value of C that allows
    the order of g mod N to be recovered via order finding. This when setting up
    a simulator, using it to sample m vectors, and post-processing the samples,
    for each value of C to test.

    Note that order finding is attempted for a single random problem instance
    only for each value of C tested. Additional tests may hence need to be
    performed for the minimum value of C returned by this function, and for
    slightly smaller and larger values of C.

    @param t  The number of distinct prime factors of N.

    @param l  The bit length l of each distinct prime factor of m.

    @param precision  The level of precision with which to perform the search.

    @param B  The bound B on the smoothness of p_i - 1, for p_1, .., p_t the t
              distinct prime factors of N.

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

    [N, F] = sample_integer(t, l, B=B, verbose=verbose)

    R = IntegerModRing(N)
    while True:
        g = R.random_element()
        if gcd(g.lift(), N) == 1:
            break;

    order_expected = find_smooth_order_mod_N(g, N, F)

    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    n = N.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    B = generate_basis_for_order_finding(
            N, F, u=[g], d=d, threads=threads, verbose=verbose
        )

    simulator = Simulator(B, verbose=verbose, threads=threads)

    log_file = None
    if None != log_file_prefix:
        # Append the run ID to the profile file prefix.
        log_file_prefix = log_file_prefix + uid_for_N(N) + "-" + uid() + "-"

        log_file = open(log_file_prefix + "order.txt", "w")

        print(
            "Start time:",
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            file=log_file,
        )
        print("Hostname:", platform.node(), "\n", file=log_file)

        print("N =", N, file=log_file)
        print("F =", F, "\n", file=log_file)

        print("g =", g, "\n", file=log_file)

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
            order_found = solve_samples_for_order(
                              samples,
                              g,
                              N,
                              R,
                              block_size=block_size,
                              profile_file=profile_file,
                              verbose=verbose)

            # Test the order.
            solved = (order_found == order_expected)
        except Exception as e:
            if verbose:
                print(f"Error: Finding r failed with exception: {e}")

            solved = False

        if log_file_prefix != None:
            # Write out the metadata file.
            meta_file = f"{log_file_prefix}meta-{simulator.d}-{m}-{C}.txt"

            with open(meta_file, "w") as f:
                print("C =", C, file=f)
                print("R =", R, "\n", file=f)

                print("success =", solved, file=f)

            samples_file = \
              f"{log_file_prefix}samples-{simulator.d}-{m}-{C}.sobj"
            save(samples, samples_file)

        return solved

    return find_minimum_C(test_for_C, precision, log_file, verbose=verbose)


def find_minimum_C_for_order_finding_phi(
    t=2,
    l=1024,
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
    @brief  Finds the minimum value of the constant C, with a given prescribed
            precision, such that the post-processing succeeds in recovering the
            order phi of the multiplicative group of the ring of integers
            mod N.

    This function first picks an integer N = p_1 * ... * p_t, with t distinct
    l-bit prime factors such that all p_i - 1 are B-smooth for B some small
    bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    It then performs a binary search to find the minimum value of C that allows
    the order of phi of the multiplicative group of the ring of integers mod N
    to be recovered via order finding. This when setting up a simulator, using
    it to sample m vectors, and post-processing the samples, for each value of
    C to test.

    Note that order finding is attempted for a single random problem instance
    only for each value of C tested. Additional tests may hence need to be
    performed for the minimum value of C returned by this function, and for
    slightly smaller and larger values of C.

    @param t  The number of distinct prime factors of N.

    @param l  The bit length l of each distinct prime factor of N.

    @param precision  The level of precision with which to perform the search.

    @param B  The bound B on the smoothness of p_i - 1, for p_1, .., p_t the t
              distinct prime factors of N.

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

    [N, F] = sample_integer(t, l, B=B, verbose=verbose)

    phi_expected = prod([p - 1 for p in F])

    if verbose:
        print("")
        print("** Setting up the simulator...")
        print("")

    n = N.nbits()
    d = ceil(dp * sqrt(n))

    if mp == None:
        m = ceil(sqrt(n)) + 4
    else:
        m = ceil(mp * sqrt(n))

    B = generate_basis_for_order_finding(
            N, F, u=[], d=d, threads=threads, verbose=verbose
        )

    simulator = Simulator(B, verbose=verbose, threads=threads)

    log_file = None
    if None != log_file_prefix:
        # Append the run ID to the profile file prefix.
        log_file_prefix = log_file_prefix + uid_for_N(N) + "-" + uid() + "-"

        log_file = open(log_file_prefix + "order-phi.txt", "w")

        print(
            "Start time:",
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            file=log_file,
        )
        print("Hostname:", platform.node(), "\n", file=log_file)

        print("N =", N, file=log_file)
        print("F =", F, "\n", file=log_file)

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
            phi_found = solve_samples_for_phi(
                            samples,
                            N,
                            R,
                            block_size=block_size,
                            profile_file=profile_file,
                            verbose=verbose)

            # Test the order.
            solved = (phi_found == phi_expected)
        except Exception as e:
            if verbose:
                print(f"Error: Finding phi failed with exception: {e}")

            solved = False

        if log_file_prefix != None:
            # Write out the metadata file.
            meta_file = f"{log_file_prefix}meta-{simulator.d}-{m}-{C}.txt"

            with open(meta_file, "w") as f:
                print("C =", C, file=f)
                print("R =", R, "\n", file=f)

                print("success =", solved, file=f)

            samples_file = \
              f"{log_file_prefix}samples-{simulator.d}-{m}-{C}.sobj"
            save(samples, samples_file)

        return solved

    return find_minimum_C(test_for_C, precision, log_file, verbose=verbose)
