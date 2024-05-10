# ------------------------------------------------------------------------------
# This Sage script (and the associated supporting scripts) implements a
# simulator for the quantum algorithm in [Regev23], alongside the classical
# post-processing algorithm from [Regev23] that recovers the factors from
# simulated samples. It furthermore implements some of the improvements to
# and extensions of [Regev23] that are described in [EG23p].
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


def build_factoring_element_vectors(d):

    """
    @brief  Builds the a = (a_1, .., a_d) and b = (b_1, .., b_d) element
            vectors, where b_1, .., b_d are the first d primes, and a_i = b_i^2.

    @param d  The dimension d of the element vectors to build.

    @return   [a, b] for a and b the element vectors built.
    """

    # Set b = [b_1, .., b_d] to the first d primes.
    b = []

    p = 2
    while len(b) < d:
        b.append(p)
        p = next_prime(p)

    # Set a = [a_1, .., a_d] to a_i = b_i^2.
    a = [b_i^2 for b_i in b]

    return [a, b]


def generate_basis_for_factoring(
    N,
    F,
    d=None,
    *,
    threads=DEFAULT_THREADS,
    verbose=False):

    """
    @brief  Generates a basis for the d-dimensional lattice L used for
            factoring, when given the dimension d, the integer N to factor and
            the factors of N.

    The integer N = p_1 * ... * p_t is a composite such that p_i - 1 is B-smooth
    for B some small bound, and such that p_i shares only a factor of two with
    p_j for i ≠ j. Factors may not occur with multiplicity in p_i - 1.

    @param N  The integer N to factor.

    @param F  A list [p_1, .., p_t] of the factors of N.

    @param d  The dimension d of the lattice L. Set to ceil(sqrt(n)), for n the
              bit length of N, if omitted.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   A basis for the lattice L used for factoring.
    """

    if None == d:
        n = N.nbits()
        d = ceil(sqrt(n))

    [a, _] = build_factoring_element_vectors(d)
    return generate_basis_for_L(N, F, a, threads=threads, verbose=verbose)


def solve_samples_for_factors(
    samples,
    N,
    R,
    *,
    block_size=DEFAULT_BLOCK_SIZE,
    profile_file=None,
    verbose=False):

    """
    @brief  Solves a list of samples for the factors of N by using the
            post-processing from [Regev23].

    [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                            ArXiv 2308.06572v2 (2023).

    As a minor improvement over the procedure in [Regev23], this function
    post-processes all relations found, so as to attempt to achieve more than a
    split. In many cases, this yields the complete factorization, depending on
    how parameter such as C, d and the number of runs m are selected. The
    exact procedure is described in [EG23p]. Compare to [E21b] for Shor's
    factoring algorithm.

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev’s factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023–2024).

    [E21b]  Ekerå, M.: "On completely factoring any integer efficiently in
                        a single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

    @param samples  The list of samples to solve for factors.

    @param N  The integer N to factor.

    @param R  The parameter R specifying the standard deviation of the noise.

    @param block_size   The blocksize for lattice reduction used during
                        post-processing. For the default value of 2, LLL is
                        used. For larger block-sizes, BKZ is used instead.

    @param profile_file   The path to a file in which to save the profile of the
                          Gram–Schmidt norms of the post-processing matrix. If
                          omitted no profile is saved.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The set of factors of N found.
    """

    # Setup and start a timer.
    timer = Timer().start()

    # Print status.
    if verbose:
        print("Post-processing the sampled vectors to find factors...")
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

    if block_size != 2:
        if verbose:
            print(" Running BKZ on the post-processing matrix...")

        denominator = M.denominator()
        X = (M * denominator).change_ring(ZZ).BKZ(
            block_size=block_size, algorithm="fpLLL", fp="rr", precision=128
        ) / denominator
    else:
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

    # Extract relevant submatrix of the post-processing matrix
    LB = X[:d, :d]

    [a, b] = build_factoring_element_vectors(d)

    R = IntegerModRing(N)

    if verbose:
        print("")
        print(" Processing the factoring relations found...")

    factors = FactorCollection(N)

    for v in LB:
        if not is_in_lattice(v, N, a):
            continue

        B = prod([R(b[i])^v[i] for i in range(d)])

        if B not in [R(-1), R(1)]:

            # We have:    B^2 = 1 (mod N) but B \not \in {-1, 1}
            #          => B^2 - 1 = 0 (mod N)
            #          => (B - 1) (B + 1) = 0 (mod N),
            #
            # where B ± 1 ≠ 0 (mod N), and so we have a split.

            factor = gcd(B.lift() - 1, N)

            if verbose:
                print("  Found factor:", factor)

            factors.add(factor)

    if verbose:
        print("")
        print(" Time required to post-process:", timer)

    return factors


def test_factoring(
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
            algorithm in [Regev23] and the associated classical post-processing.

    [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                            ArXiv 2308.06572v2 (2023).

    This function first picks an integer N = p_1 * ... * p_t, with t distinct
    l-bit prime factors such that all p_i - 1 are B-smooth for B some small
    bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    It then sets up the simulator, and uses the simulator to sample m vectors
    representative of vectors that the quantum computer would output according
    to Regev's analysis of the quantum algorithm in [Regev23].

    Finally, it solves the vectors sampled for factors of N by using Regev's
    lattice-based post-processing from [Regev23] with improvements from [EG23p].

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    @param t  The number of distinct prime factors of N.

    @param l  The bit length l of each distinct prime factor of N.

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

    @return   The factors of N found.
    """

    # Setup the problem instance.
    if verbose:
        print("** Setting up the problem instance...")
        print("")

    [N, F] = sample_integer(t, l, B=B, verbose=verbose)

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

    B = generate_basis_for_factoring(
            N, F, d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    if verbose:
        print("")
        print("** Using the simulator to sample vectors...")
        print("")

    R = get_regev_R(C, n)

    samples = simulator.sample_vectors_with_failure_rate(
        R, m=m, failure_rate=failure_rate, verbose=verbose
    )

    if verbose:
        print("")
        print("** Post-processing the sampled vectors to find factors...")
        print("")

    factors = solve_samples_for_factors(
                  samples, N, R, block_size=block_size, verbose=verbose
              )

    if verbose:
        print("")

    return factors


def find_minimum_C_for_factoring(
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
            precision, such that the post-processing succeeds in completely
            factoring a random problem instance with the given parameters.

    This function first picks an integer N = p_1 * ... * p_t, with t distinct
    l-bit prime factors such that all p_i - 1 are B-smooth for B some small
    bound, and such that gcd(p_i, p_j) = 2 for all i ≠ j.

    It then performs a binary search to find the minimum value of C that allows
    the integer N to be factored completely. This when setting up a simulator,
    using it to sample m vectors, and post-processing the samples, for each
    value of C to test.

    Note that factorization is attempted for a single random problem instance
    only for each value of C tested. Additional tests may hence need to be
    performed for the minimum value of C returned by this function, and for
    slightly smaller and larger values of C.

    @param t  The number of distinct prime factors of N.

    @param l  The bit length l of each distinct prime factor of m.

    @param precision  The level of precision with which to perform the search.

    @param B  The bound B on the smoothness of p_i - 1, for p_1, .., p_t the
              t distinct prime factors of N.

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

    B = generate_basis_for_factoring(
            N, F, d=d, threads=threads, verbose=verbose
        )

    if verbose:
        print("")

    simulator = Simulator(B, threads=threads, verbose=verbose)

    log_file = None
    if None != log_file_prefix:
        # Append the run ID to the profile file prefix.
        log_file_prefix = log_file_prefix + uid_for_N(N) + "-" + uid() + "-"

        log_file = open(log_file_prefix + "log.txt", "w")

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

        factors = solve_samples_for_factors(
                      samples,
                      N,
                      R,
                      block_size=block_size,
                      profile_file=profile_file,
                      verbose=verbose)

        if log_file_prefix != None:
            # Write out the metadata file.
            meta_file = f"{log_file_prefix}meta-{simulator.d}-{m}-{C}.txt"

            with open(meta_file, "w") as f:
                print("C =", C, file=f)
                print("R =", R, "\n", file=f)

                print("success =", factors.is_complete(), file=f)

            samples_file = \
              f"{log_file_prefix}samples-{simulator.d}-{m}-{C}.sobj"
            save(samples, samples_file)

        return factors.is_complete()

    return find_minimum_C(test_for_C, precision, log_file, verbose=verbose)
