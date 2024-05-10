# ------------------------------------------------------------------------------
# In [EG23p], an extension of Regev's factoring quantum algorithm [Regev23] to
# the discrete logarithm problem was introduced. This Sage script (and the
# associated supporting scripts) implements a simulator for the quantum
# algorithm in [Regev23] and for the extensions in [EG23p].
#
# This script is used by the factoring.sage, logarithm-finding.sage and
# order-finding.sage scripts. It is not meant to be used directly; only via the
# aforementioned scripts.
#
#   [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
#                           ArXiv 2308.06572v2 (2023).
#
#   [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev’s factoring algorithm
#                                         to compute discrete logarithms".
#                                         ArXiv 2311.05545v2 (2023–2024).

from scipy.stats import binom
from dependencies.timer import Timer
from sage.stats.distributions.discrete_gaussian_lattice import \
    DiscreteGaussianDistributionLatticeSampler

load("common.sage")
load("uids.sage")


def is_in_lattice(z, N, element_vector):

    """
    @brief  A help function for testing whether a vector z is in the lattice L
            defined by the elements vector element_vector and the modulus N.

    @param  z   The vector z to test for membership in the lattice.

    @param  N   The modulus N.

    @param  element_vector  The elements vector [g_1, .., g_d]. The lattice
                            consists of all vectors (v_1, .., v_d) such that
                            g_1^v_1 * ... * g_d^v_d = 1 (mod N).

    @return   True if the vector is in the lattice, False otherwise.
    """

    R = IntegerModRing(N)

    return prod([R(ei)^zi for (zi, ei) in zip(z, element_vector)]) == R(1)


def get_regev_R(C, n):

    """
    @brief  For a given constant C and bit length n, this function returns the
            value for the parameter R in Regev's analysis.

    The parameter R controls the standard deviation of the noise that is added
    to sampled vectors by the simulator.

    @param C  The constant C that specifies the control register lengths.

    @param n  The bit length n.

    @return   Corresponding value for the R parameter.
    """

    return ceil(2^(C * sqrt(n)))


# Note: The below function is copied from Ekerå's simulators.
def _inner_product(a, b):

    """
    @brief  Computes the inner product between two vectors a and b.

    @param a  The vector a.
    @param b  The vector b.

    @return   The inner product of a and b.
    """

    if a.dimensions() != b.dimensions():
        raise Exception("Error: Incompatible dimensions.")

    n = a.dimensions()[1]

    result = 0
    for j in range(n):
        result += a[0, j] * b[0, j]

    return result


# Note: The below function is copied from Ekerå's simulators.
def _babai_cvp(B, t, Bgs=None):

    """
    @brief  Uses Babai's nearest plane algorithm to find a vector in the lattice
            L that is close to a target vector t given a reduced basis B for L.

    @param B  The reduced basis B for L.
    @param t  The target vector t.

    @param Bgs  The Gram–Schmidt orthogonal basis for B, or None in which case
                this function will call B.gram_schmidt() to compute the
                Gram–Schmidt orthogonal basis for B.

                When calling this function repeatedly, it is advantageous to
                pre-compute the Gram–Schmidt orthogonal basis for B.

    @return   A vector in L close to the target vector t.
    """

    # Compute the Gram–Schmidt orthogonal basis of B.
    if Bgs == None:
        (Bgs, _) = B.gram_schmidt()

    # Let n be the number of rows in B.
    n = B.dimensions()[0]

    b = copy(t)
    for j in range(n - 1, -1, -1):
        bj = B[j, :]
        bjs = Bgs[j, :]

        cj = round(_inner_product(b, bjs) / _inner_product(bjs, bjs))
        b = b - cj * bj

    return t - b


def build_M_matrix(samples, S):

    """
    @brief  Builds the post-processing matrix M given noisy samples from the
            dual lattice and a scaling parameter S.

    @param samples  The list of samples from which to build the matrix.

    @param S  The scaling parameter S.

    @return   The post-processing matrix M built.
    """

    samples_matrix = samples[0]
    for sample in samples[1:]:
        samples_matrix = samples_matrix.stack(sample)

    d = samples[0].dimensions()[1]

    m = len(samples)

    basis_matrix = (
        identity_matrix(d)
        .augment(zero_matrix(d, m))
        .stack(S * (samples_matrix.augment(identity_matrix(m))))
    ).transpose()

    return basis_matrix


def sample_from_L(N, F, b, *, max_tries = 10^4):

    """
    @brief  Samples a vector from the d-dimensional lattice L that is such that
            z = (z_1, ..., z_d) is in L if b_1^z_1 * ... * b_d^z_d = 1 (mod N).

    @param N  The integer N = p_1 * ... * p_t.

    @param F  A list [p_1, .., p_t] of the factors of N. For the sampling to be
              efficient, p_i - 1 must be smooth for all p_i in F.

    @param b  The d elements [b_1, .., b_d] in the ring of integer mod N.

    @param max_tries  The number of tries to find a vector in L. Dependening on
                      the vector b and integer N, the sampling may require many
                      tries. This parameter limits this number of tries.

    @return   The d-dimensional vector sampled.
    """

    d = len(b)

    for i in range(max_tries):
        # Select all z_i in z = [z_1, .., z_d] uniformly from [0, 10^6).
        z = [IntegerModRing(10^6).random_element().lift() for _ in range(d)]

        # Attempt to select the entry z_j of z so as to force
        #
        #   prod_{i = 1}^{d} b_i^z_i = 1 (mod N).
        j = randint(0, d - 1)

        z[j] = 0

        R = IntegerModRing(N)

        # We start by computing the residual factor
        #
        #   residual = prod_{i = 1}^{d} b_i^z_i (mod N),
        #
        # and the target = residual^-1 (mod N) that we want b_j^z_j to match.

        residual = prod([R(b[i])^z[i] for i in range(d)])
        target = residual^-1

        # Check if we can find z_j that yields b_j^z_j = target (mod N).
        b_j = R(b[j])

        # Solve independently for each prime power that divides the order of the
        # element b_j. In practice, we can solve independently for each prime
        # power that divides a positive multiple of the order.
        #
        # Below, we take lcm(p_1 - 1, .., p_t - 1) as such a multiple.
        order_multiple = lcm([p - 1 for p in F])
        factors = factor(order_multiple)

        logs = []

        try:
            for [f, e] in factors:
                b_j_constrained = b_j^(order_multiple / f^e)
                target_constrained = target^(order_multiple / f^e)

                logs.append(
                    [
                        discrete_log(
                            target_constrained, b_j_constrained, ord=f^e
                        ),
                        f^e
                    ]
                )
        except ValueError:
            # Failed to solve one of the discrete logarithms: Try again with
            # new vector z and index j.
            continue


        # Use the Chinese remainder theorem to compose the solutions.
        z_j = CRT([d for [d, _] in logs], [f_pow for [_, f_pow] in logs])

        # Sanity checks.
        if b_j^z_j != target:
            raise Exception("Error: Incorrect solution.")

        # Store the solution.
        z[j] = z_j

        # Sanity check.
        if prod([R(b[i])^z[i] for i in range(d)]) != 1:
            raise Exception("Error: Failed to meet the target.")

        return z

    raise Exception("Error: Falied to find solution.")


def generate_basis_for_L(N, F, b, *, threads=DEFAULT_THREADS, verbose=False):

    """
    @brief  Generates a basis for a d-dimensional lattice L with (z_1, ..., z_d)
            in L if b_1^z_1 * ... * b_d^z_d = 1 (mod N).

    @param N  The integer N = p_1 * ... * p_t.

    @param F  A list [p_1, .., p_t] of the factors of N. For the lattice to be
              generated efficiently, p_i - 1 must be smooth for all p_i in F.

    @param b  The d elements [b_1, .., b_d] in the ring of integer mod N.

    @param threads  The number of threads to use when sampling.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The basis for the lattice L.
    """

    # Setup and start a timer.
    timer = Timer().start()
    d = len(b)

    if verbose:
        print("Generating a basis for the lattice L...")
        print(" Sampling vectors from L...")

    c = 8
    B = []
    for _ in range(d + c):
        if verbose:
            print(
                "  Sampling vectors",
                len(B) + 1,
                "to",
                min(len(B) + threads, d + c),
                "of",
                d + c,
                "using",
                threads,
                "threads..."
            )

        if threads > 1:
            # Use a parallel implementation to speed up the sampling.

            @parallel(ncpus=threads)
            def sample_from_L_in_parallel(seed):
                set_random_seed(seed)
                return sample_from_L(N, F, b)

            seeds = [generate_new_random_seed() for _ in range(threads)]
            vs = [v for v in sample_from_L_in_parallel(seeds)]
            vs = [vs[i][1] for i in range(threads)]

            for v in vs:
                B.append(v)

                if len(B) >= d + c:
                    break

            if len(B) >= d + c:
                break
        else:
            # Use a basic single-threaded implementation.
            v = sample_from_L(N, F, b)
            B.append(v)

    # Reduce the basis for L.
    if verbose:
        print("")
        print(" Reducing the basis for L, this may take a moment...")

    B = matrix(B)
    A = B.LLL()[-d:]
    # Sanity checks.
    if A.row_space() != B.row_space():
        raise Exception("Error: Failed to run LLL (1).")

    if verbose:
        print("")
        print(" Time required to generate a basis for L:", timer)

    return A


class Simulator:

    """
    @brief  A class that implements a simulator for the quantum part of Regev's
            factoring algorithm [Regev23], and Ekerå–Gärtner's extensions
            [EG23p] to discrete logarithm finding, order finding and factoring
            completely via order finding.

    [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                            ArXiv 2308.06572v2 (2023).

    [EG23p]   Ekerå, M. and Gärtner, J.: "Extending Regev's factoring algorithm
                                          to compute discrete logarithms".
                                          ArXiv 2311.05545v2 (2023).

    The simulator requires the modulus to be of special form to enable a lattice
    basis to be efficiently constructed. It may therefore only be used to
    simulate problem instances that are classically tractable.
    """

    def __init__(self, B, *, threads=DEFAULT_THREADS, verbose=False):

        """
        @brief  Initializes the simulator.

        @param B  A full rank basis for the d-dimensional lattice underlying the
                  quantum algorithm.

        @param threads  The number of threads to use when sampling.

        @param verbose  A flag that may be set to True to print verbose status
                        messages, or to False not to print such messages.
        """

        # Setup and start a timer.
        timer = Timer().start()

        self.threads = threads

        dims = B.dimensions()

        if dims[0] != dims[1] or dims[0] != B.rank():
            raise Exception("Error: Simulator require full rank matrix B")

        if verbose:
            print("Setting up the simulator...")

        self.B = B
        self.d = dims[0]

        if verbose:
            print(" Computing the basis for the dual L^* of L...")

        B_dual_t = B * (B.transpose() * B)^-1
        B_dual = B_dual_t.LLL()
        self.B_dual = B_dual

        # Sanity checks.
        if B_dual_t.row_space() != B_dual.row_space():
            raise Exception("Error: Failed to run LLL (2).")

        if verbose:
            print(
                " Computing the Gram–Schmidt orthogonalization of the basis..."
            )

        (self.Bgs, _) = self.B_dual.gram_schmidt()

        if verbose:
            print("")
            print(" Time required to setup the simulator:", timer)


    def _sample_noisy_vector(self, D, R):

        """
        @brief  Samples a vector close to a random vector in L^* / Z^d.

        This function first samples a random vector in L^* / Z^d, and then
        returns a sample from a discrete Gaussian distribution centered around
        this vector.

        @param D  The parameter D that specifies the discretization.

        @param R  The parameter R that specifies the standard deviation.

        @return   A d x 1-dimensional matrix close to L^* / Z^d.
        """

        t = matrix(
            [
                IntegerModRing(D).random_element().lift() / D
                for _ in range(self.d)
            ]
        )
        v = _babai_cvp(self.B_dual, t, self.Bgs)

        # Discrete Gaussian distribution centered around the sampled vector.
        sample_distribution = DiscreteGaussianDistributionLatticeSampler(
            ZZ^self.d, D / (2 * sqrt(pi) * R), (D * v[0]).apply_map(round)
        )

        sample = matrix(sample_distribution() % D) / D

        return sample


    def _sample_bad_vector(self, D):

        """
        @brief  Samples a bad vector, by sampling a d-dimensional vector such
                that each component is selected uniformly at random from the set
                {0, 1/D, .., (D-1)/D}.

        @param D  The parameter D that specifies the discretization.

        @return   The vector sampled.
        """

        return matrix(
            [
                IntegerModRing(D).random_element().lift() / D
                for _ in range(self.d)
            ]
        )


    def sample_vectors(self, R, m=None, verbose=False):

        """
        @brief  Samples vectors close to L^* / Z^d.

        More specifically, the sampled vectors are from a discrete Gaussian
        distribution with parameter 1 / (2 sqrt(pi) R) centered around a random
        vector in L^* / Z^d, see [Regev23]. The discrete Gaussian distribution
        is discretized to {0, 1 / D, .., (D - 1) / D}^d.

        [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                                ArXiv 2308.06572v2 (2023).

        The actual quantum algorithm would output a vector with elements in
        {0, 1, .., D - 1}. This simulator divides all components by D.

        @param R  The parameter R specifying the standard deviation.

        @param m  The number of vectors to sample. May be set to None, in which
                  case the number of samples is set to d + 4.

        @param verbose  A flag that may be set to True to print verbose status
                        messages, or to False not to print such messages.

        @return   The m vectors that are sampled from the distribution that
                  the quantum algorithm is expected to induce.
        """

        if None == m:
            m = self.d + 4

        return self.sample_vectors_good_and_bad(R, m, 0, verbose)


    def sample_vectors_with_failure_rate(
        self,
        R,
        failure_rate,
        m=None,
        verbose=False):

        """
        @brief  Samples vectors close to L^* / Z^d.

        More specifically, the sampled vectors are from a discrete Gaussian
        distribution with parameter 1 / (2 sqrt(pi) R) centered around a random
        vector in L^* / Z^d, see [Regev23]. The discrete Gaussian distribution
        is discretized to {0, 1 / D, .., (D - 1) / D}^d. Furthermore, to
        simulate error failures of the quantum computer, each vector has a
        probability to be bad, in which case it instead is sampled from the
        uniform distribution over {0, 1 / D, .., (D - 1) / D}^d.

        [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                                ArXiv 2308.06572v2 (2023).

        The actual quantum algorithm would output a vector with elements in
        {0, 1, .., D - 1}. This simulator divides all components by D.

        @param R  The parameter R specifying the standard deviation.

        @param m  The number of vectors to sample. May be set to None, in which
                  case the number of samples is set to

                        (d + 4) / (1 - failure_rate)

                  rounded up.

        @param failure_rate   The failure rate on [0, 1). May be set to a value
                              greater than zero to simulate error-correction
                              failures resulting in bad vectors being output.

        @param verbose  A flag that may be set to True to print verbose status
                        messages, or to False not to print such messages.

        @return   The m vectors that are sampled from the distribution that
                  the quantum algorithm is expected to induce.
        """

        if m == None:
            m = ceil((self.d + 4) / (1 - failure_rate))

        bad_runs = binom(m, failure_rate).rvs()
        good_runs = m - bad_runs

        return self.sample_vectors_good_and_bad(
                  R, good_runs, bad_runs, verbose=verbose
               )


    def sample_vectors_good_and_bad(self, R, m1, m2=0, verbose=False):

        """
        @brief  Samples vectors close to L^* / Z^d.

        More specifically, the sampled vectors are from a discrete Gaussian
        distribution with parameter 1 / (2 sqrt(pi) R) centered around a random
        vector in L^* / Z^d, see [Regev23]. The discrete Gaussian distribution
        is discretized to {0, 1 / D, .., (D - 1) / D}^d.

        Additionally, to simulate error-correct failures, the sampled vectors
        also include m2 bad vectors sampled from the uniform distribution over
        {0, 1 / D, .., (D - 1) / D}^d.

        [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                                ArXiv 2308.06572v2 (2023).

        The actual quantum algorithm would output a vector with elements in
        {0, 1, .., D - 1}. This simulator divides all components by D.

        @param R  The parameter R specifying the standard deviation.

        @param m1   The number of good vectors to sample.

        @param m2   The number of bad vectors to sample.

        @param verbose  A flag that may be set to True to print verbose status
                        messages, or to False not to print such messages.

        @return   The m = m1 + m2 good and bad vectors sampled, shuffled so as
                  to randomize the order of the good and bad vectors.
        """

        # Setup and start a timer.
        timer = Timer().start()

        d = self.d
        D = 2^(ceil(log(2 * sqrt(d) * R) / log(2)))

        # Sanity checks.
        if not (R > sqrt(2 * d)):
            raise Exception("Error: Sanity check error (1).")

        if not D.is_power_of(2):
            raise Exception("Error: Sanity check error (2).")

        if not (2 * sqrt(d) * R <= D <= 4 * sqrt(d) * R):
            raise Exception("Error: Sanity check error (3).")

        if verbose:
            print(f"Sampling {m1} good vectors and {m2} bad vectors...")

        samples = []

        while len(samples) < m1:
            vs = []
            threads = min(self.threads, m1 - len(samples))
            if threads > 1:
                if verbose:
                    print(
                        " Sampling vectors",
                        len(samples) + 1,
                        "to",
                        len(samples) + threads,
                        "of",
                        m1,
                        "using",
                        threads,
                        "threads..."
                    )

                # Use a parallel implementation to speed up the sampling.
                @parallel(ncpus=threads)
                def sample_vector_in_parallel(seed):
                    set_random_seed(seed)
                    return self._sample_noisy_vector(D, R)

                seeds = [generate_new_random_seed() for _ in range(threads)]
                vs = [v[1] for v in sample_vector_in_parallel(seeds)]
            else:
                # Use a basic single-threaded implementation.
                vs = [self._sample_noisy_vector(D, R)]

            for sample in vs:
                samples.append(sample)

        while len(samples) < m1 + m2:
            samples.append(self._sample_bad_vector(D))

        # Shuffle samples so bad and good vectors are mixed
        shuffle(samples)

        if verbose:
            print("")
            print(" Time required to sample:", timer)

        return samples


    def __repr__(self):

        """
        @brief  Represents the simulator as a string.

        @return   A string representation of the simulator.
        """

        return f"Simulator for a {self.d}-dimensional lattice."
