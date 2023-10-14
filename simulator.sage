# ------------------------------------------------------------------------------
# This Sage script implements a simulator for the quantum algorithm in:
#
#   [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
#                           ArXiv 2308.06572v2 (2023).
#
# It furthermore implements the post-processing from [Regev23] with some minor
# tweaks to attempt to use all available relations to not only split N, but to
# factor N a completely as possible. In many cases, the available information is
# enough to completely factor N, depending on how the parameters are selected.

from timer import Timer;

def sample_smooth_prime(l, B = 10^5, factors = set()):

  """ @brief  Samples an l-bit prime p such that p - 1 is B-smooth, and such
              that p - 1 has no factor that is in the factors set beside two.

      @param l  The bit length of each prime p.
      @param B  The bound B on the B-smoothness of p - 1.

      @param factors  The set of prime factors less than B already used. If two
                      is a part of this set, two will be ignored, since p - 1
                      must be even when p is a large prime.

      @return   [p, factors], for p the prime selected, and factors an updated
                set of used factors that also contains the factors of p - 1. """

  # Create a list F of all prime factors up to B not in factors.
  F = [];
  for p in primes(3, B):
    if p not in factors:
      F.append(p);

  # Search exhaustively for a combination of factors that yields B-smooth p - 1
  # such that p is an l-bit prime and return this prime.
  while True:
    shuffle(F);
    used = set();

    x = 2;
    used.add(2);

    for p in F:
      x *= p;
      used.add(p);

      if x >= 2^(l-1) <= x < 2^l:
        p = x + 1;

        if p.is_prime(proof = False):
          return [p, factors.union(used)];

      if x >= 2^l:
        break;

def sample_integer(m, l, B = 10^5, verbose = True):

  """ @brief  Samples m distinct l-bit primes p_1, .., p_m such that p_i - 1 is
              B-smooth for all i in [1, m], and gcd(p_i - 1, p_j - 1) = 2 for
              any choice of distinct i, j in [1, m].

      @param m  The number of distinct l-bit primes p_1, .., p_m.
      @param l  The bit length of each prime p_1, .., p_m.
      @param B  The bound B on the B-smoothness of p_1 - 1, .., p_m - 1.

      @param verbose  A flag that may be set to True to print verbose status
                      messages, or to False not to print such messages.

      @return   [N, [p_1, .., p_m]] for N = p_1 * .. * p_m. """

  # Setup and start a timer.
  timer = Timer().start();

  if verbose:
    print("Sampling N on special form to enable efficient simulation...");

  # Primes.
  primes = set();
  factors = set();

  # Pick m primes p_i of length l bits each, all such that p_i is B-smooth.
  for _ in range(m):
    while True:
      [p, factors] = sample_smooth_prime(l, B, factors);

      if p not in primes:
        primes.add(p);

        if verbose:
          print(" Sampled factor:", p);

        break;

  # Compute the product.
  N = prod(primes);

  # Stop the timer.
  timer.stop();

  if verbose:
    print("");
    print(" Sampled N =", N);
    print("");
    print(" Time required to sample:", timer);

  return [N, primes];

def sample_from_L(N, F, d):

  """ @brief  Samples a vector from the lattice L, given N = p_1 * .. * p_m and
              F = [p_1, .., p_m] where p_1 - 1, .., p_m - 1 are B-smooth.

      @param N  The integer N = p_1 * .. * p_m.
      @param F  The factors F = [p_1, .., p_m] of N.
      @param d  The dimension d of the vector to be sampled.

      @return   The d-dimensional vector sampled. """

  # Compute the Carmichael function of N.
  lambda_N = lcm(list(F));

  # Set b = [b_1, .., b_d] and a = [a_1, .., a_d], where a_i = b_i^2, and where
  # the sequence b_1, .., b_d are the first d primes.
  [a, _] = build_a_b_vectors(d);

  # Select the first d - 1 entries of z = [z_1, .., z_d] uniformly at random
  # from the interval [0, \lambda(N)).
  z = [IntegerModRing(lambda_N).random_element().lift() for _ in range(d - 1)];

  # Attempt to select the last entry z_d of z so as to force
  #
  #   prod_{i = 1}^{d} a_i^z_i = 1  (mod N).

  R = IntegerModRing(N);

  # We start by computing the residual factor
  #
  #   residual = prod_{i = 1}^{d - 1} a_i^z_i (mod N),
  #
  # and the target = residual^-1 (mod N) that we would like a_d^z_d to match.

  residual = prod([R(a[i])^z[i] for i in range(d - 1)]);
  target = residual^-1;

  # Check if we can find z_d that yields a_d^z_d = target (mod N).
  a_d = R(a[d - 1]);

  # Solve independently for each factor p of N.
  z_d_parts = [];

  for p in F:
    # Reduce everything mod p.
    R_mod_p = IntegerModRing(p);
    a_d_mod_p = R_mod_p(a_d);
    target_mod_p = R_mod_p(target);

    # Solve independently for each factor f of p - 1.
    logs = [];

    for [f, _] in factor(p - 1):
      a_d_mod_p_constrained = a_d_mod_p^((p - 1) / f);
      target_mod_p_constrained = target_mod_p^((p - 1) / f);
      logs.append([discrete_log(target_mod_p_constrained,
                                a_d_mod_p_constrained), f]);

    # Use the Chinese remainder theorem to compose the solutions.
    z_d_mod_p_1 = CRT([d for [d, _] in logs], [f for [_, f] in logs]);

    # Sanity checks.
    if a_d_mod_p^z_d_mod_p_1 != target_mod_p:
      raise Exception("Error: Incorrect solution.");

    residual_mod_p = R_mod_p(residual);
    if residual_mod_p * a_d_mod_p^z_d_mod_p_1 != 1:
      raise Exception("Error: Incorrect solution.");

    # Store the solution.
    z_d_parts.append([z_d_mod_p_1, p - 1]);

  # Use the Chinese remainder theorem to compose the solutions.
  z_d = CRT([d for [d, _] in z_d_parts], [f / 2 for [_, f] in z_d_parts]);
  z.append(z_d);

  # Sanity check.
  if a_d^z_d != target:
    raise Exception("Error: Failed to meet the target (1).");

  if prod([R(a[i])^z[i] for i in range(d)]) != 1:
    raise Exception("Error: Failed to meet the target (2).");

  # Return the solution.
  return z;

def build_a_b_vectors(d):

  """ @brief  Builds the a = (a_1, .., a_d) and b = (b_1, .., b_d) vectors,
              where b_1, .., b_d are the first d primes, and a_i = b_i^2.

      @param d  The dimension d of the vectors to build.

      @return   [a, b] for a and b the vector built. """

  # Set b = [b_1, .., b_d] to the first d primes.
  b = [];

  p = 2;
  while len(b) < d:
    b.append(p);
    p = next_prime(p);

  # Set a = [a_1, .., a_d] to a_i = b_i^2.
  a = [b_i^2 for b_i in b];

  return [a, b];

def build_M_matrix(samples, S):

  """ @brief  Builds the post-processing matrix M given noisy samples from the
              dual lattice and a scaling parameter S.

      @param sample   The list of samples.
      @param S        The scaling parameter S.

      @return   The post-processing matrix M. """

  samples_matrix = samples[0];
  for sample in samples[1:]:
    samples_matrix = samples_matrix.stack(sample);

  d = samples[0].dimensions()[1];

  m = len(samples);

  basis_matrix = identity_matrix(d).\
                    augment(zero_matrix(d, m)).\
                    stack(S*(samples_matrix.augment(identity_matrix(m))));

  return basis_matrix;

class Simulator:

  """ @brief  A class that implements a simulator for Regev's factoring
              algorithm for integers N on special form, allowing the quantum
              part of the algorithm to be simulated for large cryptographically
              relevant (but insecure) integers N. """

  def __init__(self, N, F, C = 2, c = 8, threads = 8, verbose = True):

    """ @brief  Initializes the simulator for Regev's factoring algorithm.

        @param N  The integer N = p_1 * .. * p_m of length n bits.
        @param F  The factors F = [p_1, .., p_m] of N.

        @param C  The constant C such that R = ceil(2^(C * sqrt(n))).

        @param c  The constant c that controls how many vectors are sampled to
                  construct the basis for L. In total d + c vectors are sampled,
                  for d the dimension of each vector.

        @param threads  The number of threads to use to sample vectors from L.

        @param verbose  A flag that may be set to True to print verbose status
                        messages, or to False not to print such messages. """

    # Setup and start a timer.
    timer = Timer().start();

    # Print status.
    if verbose:
      print("Setting up the simulator...");

    self.N = N;
    self.F = F;

    self.C = C;

    # Get the bit length n of N.
    self.n = n = N.nbits();

    self.d = d = ceil(sqrt(n));
    self.R = R = ceil(2^(C * sqrt(n)));

    # Sanity check.
    if not (R > sqrt(2 * d)):
      raise Exception("Error: Internal error (1).");

    self.D = D = 2^(ceil(log(2 * sqrt(d) * R) / log(2)));

    # Sanity checks.
    if not (2 * sqrt(d) * R <= D <= 4 * sqrt(d) * R):
      raise Exception("Error: Internal error (2).");

    if not D.is_power_of(2):
      raise Exception("Error: Internal error (3).");

    # Print status.
    if verbose:
      print(" Picking parameters...");
      print("  n =", n);
      print("  d =", d);
      print("");

      print("  C =", C);
      print("  R =", R);
      print("  D =", D);
      print("");

    # Sample vectors from L.
    B = [];

    # Print status.
    if verbose:
      print(" Sampling vectors from L...");

    for _ in range(d + c):
      if verbose:
        print("  Sampling vectors", len(B) + 1, "to", \
                min(len(B) + threads, d + c), "of", d + c, \
                  "using", threads, "threads...");

      if threads > 1:
        # Use a parallel implementation to speed up the sampling.

        @parallel(ncpus = threads)
        def sample_from_L_in_parallel(seed):
          set_random_seed(seed);
          return sample_from_L(N, F, d);

        seeds = [(len(B) + i) for i in range(threads)];
        vs = [v for v in sample_from_L_in_parallel(seeds)];
        vs = [vs[i][1] for i in range(threads)];

        for v in vs:
          B.append(v);

          if len(B) >= d + c:
            break;

        if len(B) >= d + c:
          break;
      else:
        # Use a basic single-threaded implementation.
        v = sample_from_L(N, F, d);
        B.append(v);

    # Reduce the basis for L.
    if verbose:
      print("")
      print(" Reducing the basis for L, this may take a moment...");

    B = matrix(B);
    A = B.LLL()[-d:];
    self.B = A;

    # Sanity checks.
    if A.row_space() != B.row_space():
      raise Exception("Error: Failed to run LLL (1).");

    if not (abs(A.det()) < 2^n):
      raise Exception("Error: The determinant of the lattice is too large.")

    # Compute and reduce the basis for the dual L^* of L.
    if verbose:
      print(" Computing the basis for the dual L^* of L...");

    B_dual_t = A * (A.transpose() * A)^-1;
    B_dual = B_dual_t.LLL();
    self.B_dual = B_dual;

    # Sanity checks.
    if B_dual_t.row_space() != B_dual.row_space():
      raise Exception("Error: Failed to run LLL (2).");

    # Pre-compute the Gram–Schmidt orthogonalization of the basis.
    if verbose:
      print(" Computing the Gram–Schmidt orthogonalization of the basis...");

    (self.Bgs, _) = self.B_dual.gram_schmidt();

    if verbose:
      print("");
      print(" Time required to setup the simulator:", timer);

  def sample_vectors(
        self,
        eta = None,
        failure_rate = 0,
        verbose = True):

    """ @brief  Samples vectors from L^* / Z^d perturbed by Gaussian noise of
                standard deviation 1 / (2 sqrt(pi) R), see [Regev23]. The
                vectors are discretized to {0, 1 / D, .., (D - 1) / D}^d.

        [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                                ArXiv 2308.06572v2 (2023).

        The actual quantum algorithm would output a vector with elements in
        {0, 1, .., D - 1}. This simulator divides all components by D.

        @param eta  The number of vectors to sample. May be set to None, in
                    which case the number of samples is set to sqrt(n) + 4
                    rounded up, for n the bit length of the integer N to factor.

        @param failure_rate   The failure rate on [0, 1). May be set to a value
                              greater than zero to simulate error correction
                              failures resulting in bad vectors being output.

        @param verbose  A flag that may be set to True to print verbose status
                        messages, or to False not to print such messages.

        @return   [samples, [d, R]] for samples the vectors sampled, d the
                  dimension of the vectors, and R a parameter that controls
                  the standard deviation of the noise. """

    # Setup and start a timer.
    timer = Timer().start();

    # If eta is set to None, pick eta = d + 4.
    if None == eta:
      eta = ceil(sqrt(self.n)) + 4;

    # Note: The below function is copied from Ekerå's simulators.
    def inner_product(a, b):

      """ @brief  Computes the inner product between two vectors a and b.

          @param a  The vector a.
          @param b  The vector b.

          @return   The inner product of a and b. """

      if a.dimensions() != b.dimensions():
        raise Exception("Error: Incompatible dimensions.");

      n = a.dimensions()[1];

      result = 0;
      for j in range(n):
        result += a[0, j] * b[0, j];

      return result;

    # Note: The below function is copied from Ekerå's simulators.
    def babai_cvp(B, t, Bgs = None):

      """ @brief  Uses Babai's nearest plane algorithm to find a vector in the
                  lattice L that is close to a target vector t given a reduce
                  basis B for L.

          @param B    The reduced basis B for L.
          @param t    The target vector t.

          @param Bgs  The Gram–Schmidt orthogonal basis for B, or None in which
                      case this function will call B.gram_schmidt() to compute
                      the Gram–Schmidt orthogonal basis for B.

                      When calling this function repeatedly, it is advantageous
                      to pre-compute the Gram–Schmidt orthogonal basis for B.

          @return   A vector in L close to the target vector t. """

      # Compute the Gram–Schmidt orthogonal basis of B.
      if Bgs == None:
        (Bgs, _) = B.gram_schmidt();

      # Let n be the number of rows in B.
      n = B.dimensions()[0]

      b = copy(t);
      for j in range(n - 1, -1, -1):
        bj = B[j, :];
        bjs = Bgs[j, :];

        cj = round(inner_product(b, bjs) / inner_product(bjs, bjs));
        b = b - cj * bj;

      return t - b;

    # Print status message.
    if verbose:
      print("Sampling vectors...");

    samples = [];

    # Append eta samples to the samples vector.
    for j in range(eta):

      # Simulate failures in the error correction.
      while True:
        pivot = RealField(90).random_element(0, 1);
        if pivot > 0:
          break;

      if failure_rate > pivot:
        bad_sample = \
              matrix([IntegerModRing(self.D).random_element().lift() / self.D \
                        for _ in range(self.d)]);
        samples.append(bad_sample);

        # Print status message.
        if verbose:
          print(" Sampled a bad vector for sample", j + 1, \
                  "of", eta, "samples...");

        continue;

      # Sample a vector from the dual via Babai's algorithm.
      t = matrix([IntegerModRing(self.D).random_element().lift() / self.D
                    for _ in range(self.d)]);
      v = babai_cvp(self.B_dual, t, self.Bgs);

      noise_dist = RealDistribution('gaussian', 1 / (2 * sqrt(pi) * self.R));

      # Sample the noise.
      noise = matrix(RealField(2 * self.D.nbits()), \
                      [noise_dist.get_random_element() for _ in range(self.d)]);

      # Add the noise and discretize.
      def coerce(x):
        return (round(self.D * x) % self.D) / self.D;

      # Note that since we have forced 2 * D.nbits() precision when converting
      # the reals to a noise vector, we avoid losing precision in the below step
      # where we add the noise vector to v (resulting in a vector of reals).

      noisy_v = v + noise;
      noisy_v = matrix([coerce(noisy_v[0, i]) for i in range(self.d)]);

      samples.append(noisy_v);

      # Print status message.
      if verbose:
        print(" Sampled a good vector for sample", j + 1, \
                "of", eta, "samples...");

    if verbose:
      print("");
      print(" Time required to sample:", timer);

    # Return the samples.
    return [samples, [self.d, self.R]];

  def __repr__(self):

    """ @brief  Represents the simulator as a string.

        @return   A string representation of the simulator. """

    return "Simulator for N = " + str(self.N) + " with C = " + str(self.C);

def solve_samples_for_factors(samples, N, d, R, verbose = True):

  """ @brief  Solves a set of vectors sampled for the factors of N by using the
              post-processing from [Regev23].

      [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                              ArXiv 2308.06572v2 (2023).

      A minor improvement is to process all relations found, so as to attempt
      to achieve more than a split. In many cases, this yields the complete
      factorization, depending on how parameter such as C and the number of
      runs eta are selected. Compare to [E21b] for Shor's factoring algorithm.

      [E21b]  Ekerå, M.: "On completely factoring any integer efficiently in
                          a single run of an order-finding algorithm".
                          Quantum Inf. Process. 20(6):205 (2021).

      @param samples  The vectors sampled.

      @param N  The integer N to factor.

      @param d  The dimension of the vectors sampled.

      @param R  A parameter that controls the level of noise in the vectors
                sampled. More specifically, the vectors sampled are perturbed by
                Gaussian noise of standard deviation 1 / (2 sqrt(pi) R).

      @param verbose  A flag that may be set to True to print verbose status
                      messages, or to False not to print such messages.

      @return   The set of factors of N found. """

  # Setup and start a timer.
  timer = Timer().start();

  # Print status.
  if verbose:
    print("Post-processing the sampled vectors to find factors...");

  # Build the post-processing matrix from the noisy discretized vectors.
  if verbose:
    print(" Building the post-processing matrix...");

  # Pick S as described in [Regev23]: In Corollary 4.5 of [Regev23], S = δ^-1,
  # and δ = sqrt(d) / (sqrt(2) R) on p. 7 of [Regev23].
  S = ceil(sqrt(2 / d) * R);
  M = build_M_matrix(samples, S);

  # Run LLL on the post-processing matrix and extract the relevant submatrix.
  if verbose:
    print(" Running LLL on the post-processing matrix...");

  X = M.transpose().LLL();
  if X.row_space() != M.transpose().row_space():
    raise Exception("Error: Failed to run LLL (3).");

  LB = X[:d, :d].transpose();

  [a, b] = build_a_b_vectors(d);

  # Setup the ring of integers modulo N.
  R = IntegerModRing(N);

  if verbose:
    print("");
    print(" Processing the factoring relations found...");

  factors = FactorCollection(N);

  for v in LB.transpose():
    A = prod([R(a[i])^v[i] for i in range(d)]);
    B = prod([R(b[i])^v[i] for i in range(d)]);

    if (A == 1) and (B not in [R(-1), R(1)]):

      # We have:    B^2 = 1 (mod N) but B \not \in {-1, 1}
      #          => B^2 - 1 = 0 (mod N)
      #          => (B - 1) (B + 1) = 0 (mod N),
      #
      # where B ± 1 ≠ 0 (mod N), and so we have a split.

      factor = gcd(B.lift() - 1, N);

      if verbose:
        print("  Found factor:", factor);

      factors.add(factor);

  if verbose:
    print("");
    print(" Time required to post-process:", timer);

  # Return the factors.
  return factors;

def solve_samples_for_factors_exhaust(samples, c, N, d, R, verbose = True):

  """ @brief  Solves all subset of eta, eta - 1, .., eta - c vectors from a set
              of eta vectors sampled for the factors of N by using the
              post-processing from [Regev23], stopping as soon as at least
              one non-trivial factor is returned.

      [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                              ArXiv 2308.06572v2 (2023).

      The idea is that if error correction failures results in some of the
      eta runs being bad, then we can exclude subsets, similar to how the
      classical lattice-based post-processing works in [EH17]. This function
      simply calls the solve_samples_for_factors() function on each subsets.

      [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                        Discrete Logarithms and Factoring RSA
                                        Integers.". In: PQCrypto 2017.
                                        Springer LNCS 10346, pp. 347–363 (2017).

      See the documentation of solve_samples_for_factors() for further details.

      Note that it would be possible to continue to collect factors from other
      subsets instead of aborting as soon as at least one non-trivial factor is
      returned if the goal is to completely factor general integers.

      @param samples  The vectors sampled.

      @param c  An integer on [0, eta).

      @param N  The integer N to factor.

      @param d  The dimension of the vectors sampled.

      @param R  A parameter that controls the level of noise in the vectors
                sampled. More specifically, the vectors sampled are perturbed by
                Gaussian noise of standard deviation 1 / (2 sqrt(pi) R).

      @param verbose  A flag that may be set to True to print verbose status
                      messages, or to False not to print such messages.

      @return   The set of factors of N found. """

  eta = len(samples);

  # Sanity check.
  if (c < 0) or (c >= eta):
    raise Exception("Error: Incorrect parameter c.")

  ts = [t for t in range(eta - c, eta + 1)];
  ts.reverse();

  for t in ts:
    for index_subset in Subsets(Set([i for i in range(eta)]), t):
      if verbose:
        print("** Trying a", len(index_subset), "sample subset:", index_subset);
        print("");

      sample_subset = [samples[i] for i in index_subset];
      factors = solve_samples_for_factors(sample_subset, N, d, R, verbose);

      if factors != FactorCollection(N):
        return factors;

      if verbose:
        print("");

def test(m = 2, l = 512, verbose = True):

  """ @brief  A convenience function for testing the simulator for the quantum
              algorithm in [Regev23] and the classical post-processing.

      [Regev23]   Regev, O.: "An Efficient Quantum Factoring Algorithm".
                              ArXiv 2308.06572v2 (2023).

      This function first picks an integer N with m distinct l-bit prime factors
      (on special form, so as to enable efficient simulation).

      It then sets up the simulator, uses the simulator to sample vectors
      representative of vectors that the quantum computer would output according
      to Regev's analysis of the quantum algorithm in [Regev23].

      Finally, it solves the vectors sampled for factors of N by using Regev's
      lattice-based post-processing from [Regev23].

      @param m  The number of distinct prime factors m of N.
      @param l  The bit length l of each distinct prime factor of m.

      @param verbose  A flag that may be set to True to print verbose status
                      messages, or to False not to print such messages.

      @return   The factors of N found. """

  # Sample the integer.
  [N, F] = sample_integer(m, l, verbose = verbose);

  # Setup a simulator.
  if verbose:
    print("");
    print("** Setting up the simulator...");
    print("");

  simulator = Simulator(N, F, verbose = verbose);

  # Use the simulator to sample vectors.
  if verbose:
    print("");
    print("** Using the simulator to sample vectors...");
    print("");

  [samples, [d, R]] = simulator.sample_vectors(verbose = verbose);

  # Solve the vectors sampled for the factors.
  if verbose:
    print("");
    print("** Post-processing the sampled vectors to find factors...");
    print("");

  factors = solve_samples_for_factors(samples, N, d, R, verbose = verbose);

  if verbose:
    print("");

  return factors;