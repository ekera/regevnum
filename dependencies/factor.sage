# ------------------------------------------------------------------------------
# Note: This file was copied from the factoritall repository by M. Ekerå.
#
# For further details, please see: https://github.com/ekera/factoritall
#
# Reference to other files and documentation below pertain to the above GitHub
# repository; not to the regevnum GitHub repository.

# ------------------------------------------------------------------------------
# This Sage script implements the procedure described in the paper:
#
# [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a single
#                    run of an order-finding algorithm".
#                   Quantum Inf. Process. 20(6):205 (2021).
#
# Use factor_completely(r, N, c = 1) to solve for r the order and N the integer.
#
# Note: This implementation assumes that the random_element() function (in the
# IntegerModRing class that is provided by Sage) is indistinguishable from a
# function that selects an element uniformly at random from the ring.

from dependencies.timer import Timer;

# An enumeration of optimization options for how factors are processed.
#
# For further details, see "optimizations.md" and Section 3.2.1 of [E21b].
class OptProcessCompositeFactors:
  # Select x uniformly at random from Z_N^*, for N the number to be factored,
  # and exponentiate x modulo N to 2^t o.
  #
  # This is as described in the algorithm in Section 3.2 of [E21b].
  JOINTLY_MOD_N = 1;

  # Select x uniformly at random from Z_N'^*, for N' the product of all pairwise
  # coprime composite factors of N currently stored in the collection, and
  # exponentiate x modulo N' to 2^t o.
  #
  # This is as above, but with optimizations from Section 3.2.1 of [E21b].
  JOINTLY_MOD_Np = 2;

  # Select x uniformly at random from Z_N'^*, for N' the product of all pairwise
  # coprime composite factors of N currently stored in the collection.
  # Exponentiate x modulo N' to 2^t o, as N' runs over the pairwise coprime
  # composite factors of N currently stored in the collection.
  #
  # This is as above, but with more optimizations from Section 3.2.1 of [E21b].
  SEPARATELY_MOD_Np = 3;

# Supporting class to collect the non-trivial factors of N on reduced form.
class FactorCollection:
  def __init__(self, N):
    # The integer N to be factored.
    self.N = N;

    # The set of factors found thus far, reduced so that all factors in the set
    # are pairwise coprime to each other. This property is enforced by add().
    self.found_factors = set();

    # The set of prime factors found thus far; a subset of found_factors.
    self.found_primes = set();

    # A timer for measuring the time spent performing primality tests.
    self.timer_test_primality = Timer();

    # A timer for measuring the time spent detecting perfect powers.
    self.timer_test_perfect_power = Timer();

    # The residual; the product of the composite pairwise coprime factors in the
    # collection, or one if there are no composite factors in the collection.
    self.residual = 1;

    # Add N as a factor.
    self.add(N);

  # Checks if all prime factors have been found.
  def is_complete(self):
    return self.residual == 1;

  # Adds a factor to this collection.
  def add(self, d):
    # Check that the factor is non-trivial and has not already been found.
    if (d == 1) or (d in self.found_factors):
      return;

    # Test if d shares a factor with any of the factors found.
    D = 1;

    for f in self.found_factors:
      D = gcd(f, d);

      if D != 1:
        break;

    if D != 1:
      # If so, remove f, split f and d, and add the resulting factors.
      self.found_factors.remove(f);
      if f not in self.found_primes:
        # Also remove f from the residual when removing f from the collection.
        self.residual /= f;

      f /= D;
      d /= D;

      self.add(D);

      if f != 1:
        self.add(f);

      if d != 1:
        self.add(d);
    else:
      # Check if d is a perfect power, and if so reduce d.
      self.timer_test_perfect_power.start();
      (d, _) = ZZ(d).perfect_power();
      self.timer_test_perfect_power.stop();

      # Add d to the factors found.
      self.found_factors.add(d);

      # Check if d is prime, and if so register it.
      self.timer_test_primality.start();
      result = d.is_prime(proof = False);
      self.timer_test_primality.stop();

      if result:
        self.found_primes.add(d);
      else:
        # If d is not prime, multiply d onto the residual.
        self.residual *= d;

  # Prints status information for this collection.
  def print_status(self):
    print("Found factors:", len(self.found_factors));
    print("Found primes:", len(self.found_primes));

    found_factors = list(self.found_factors);
    found_factors.sort();

    for i in range(len(found_factors)):
      print(" Factor " + str(i) + ":", found_factors[i]);
    print("");

  # Compares this collection to another collection.
  def __eq__(self, other):
    if self.N != other.N:
      return False;
    if self.residual != other.residual:
      return False;
    if self.found_factors != other.found_factors:
      return False;
    if self.found_primes != other.found_primes:
      return False;

    return True; # Note: Disregards timer differences.

  # Represents this collection as a string.
  def __repr__(self):
    return str(self.found_factors);

# An exception that is raised to signal an incomplete factorization. This occurs
# only if an iteration or timeout limit has been specified.
class IncompleteFactorizationException(Exception):
  def __init__(self, message, factors):
    super().__init__(message);
    self.factors = factors;

# ------------------------------------------------------------------------------
# Solves a problem instance given by r and N.
#
# The parameter c is as described in [E21b]. The parameter k in [E21b] need not
# be explicitly specified: By default, as many iterations k as are necessary to
# completely factor N will be performed. The algorithm will then stop.
#
# If you wish, you may specify k and/or a timeout in seconds. If the number of
# iterations performed exceeds k, or if the timeout is exceeded, an exception of
# type IncompleteFactorizationException will be raised.
#
# The remaining arguments are optimization flags. They are documented below in
# the code, and in "optimizations.md". It is recommended to use the defaults.
#
# This function returns the set of all distinct prime factors that divide N.
def factor_completely(r, N, c = 1,
  k = None,
  timeout = None,
  opt_split_factors_with_multiplicity = True,
  opt_report_accidental_factors = True,
  opt_abort_early = True,
  opt_square = True,
  opt_exclude_one = True,
  opt_process_composite_factors =
    OptProcessCompositeFactors.SEPARATELY_MOD_Np):

  # Sanity checks.
  if (r < 1) or (N < 2) or (c < 1):
    raise Exception("Error: Incorrect parameters.");

  print("Solving for the complete factorization...")
  print("")

  # Supporting function to build the product of q^e, for q all primes <= B and e
  # the largest exponent such that q^e <= B for B some bound.
  def build_prime_power_product(B):
    factor = 1;

    for q in prime_range(B + 1):
      e = 1;
      while q^(e + 1) <= B:
        e += 1;
      factor *= q^e;

    return factor;

  # Supporting function for computing t such that x = 2^t o for o odd.
  def kappa(x):
    if x == 0:
      return 0;

    t = 0;

    while (x % 2) == 0:
      t += 1;
      x /= 2;

    return t;

  # Note: Step 1 is already completed.
  r = ZZ(r);
  N = ZZ(N);
  m = N.nbits();

  # Setup and start a timer to measure the total time required to solve.
  timer = Timer().start();

  # Setup and reset a timer to measure the time spent exponentiating.
  timer_exponentiation = Timer();

  # Step 2: Build the product of prime factors q^e < cm and multiply onto r.
  rp = build_prime_power_product(c * m) * r;

  # Step 3: Let rp = 2^t o for o odd.
  t = kappa(rp);
  o = rp / 2^t;

  # Define a pairwise coprime set and add in N.
  F = FactorCollection(N);

  # Optimization: Initially split N when factors of N occur with multiplicity.
  #
  # If p^e divides N for e > 1, then p^(e - 1) is likely to divide r. We may use
  # this fact to initially split N when prime factors occur with multiplicity.
  #
  # Note that splitting N in this way is advantageous, as it can be done without
  # exponentiating, and as it may speed up the subsequent exponentiations (when
  # opt_process_composite_factors is not set to JOINTLY_MOD_N).
  #
  # For further details, see "optimizations.md" and [GLMS15].
  if opt_split_factors_with_multiplicity:
    d = gcd(r, N);
    if d != 1:
      print("Note: Splitting N by gcd(r, N) before commencing to iterate...\n");

      F.add(d);

  # Step 4: For j = 1, 2, ... up to k where k is unbounded.
  j = 0;

  while True:
    # Print current status information before proceeding.
    print("Iteration:", j);
    F.print_status();

    # Check if we are done...
    if F.is_complete():
      break;

    # Increment j for the next iteration.
    j += 1;

    # Check if j > k, if k is specified, and if so raise an exception passing
    # along the factors that have been found thus far.
    if (k != None) and (j > k):
      raise IncompleteFactorizationException(
        "Error: The iteration limit has been exceeded.", F.found_factors);

    # Check if the timeout is exceeded, if specified, and if so raise an
    # exception passing along the factors that have been found thus far.
    if (timeout != None) and (timer.peek() > timeout):
      raise IncompleteFactorizationException(\
        "Error: The timeout limit has been exceeded.", F.found_factors);

    # Step 4.1: Select x uniformly at random from Z_N^*.

    # Optimization: Select x uniformly at random from Z_N'^*, for N' the product
    # of all pairwise coprime composite factors of N stored in the collection.
    # This as opposed to selecting x uniformly at random from Z_N'^*, for
    # N' = N, when not applying the optimization.
    #
    # For details, see "optimizations.md" and Section 3.2.1 of [E21b].
    if opt_process_composite_factors == \
      OptProcessCompositeFactors.JOINTLY_MOD_N:
      # Let N' (denoted Np in the code) be N when not applying the optimization.
      Np = N;
    elif opt_process_composite_factors in \
      [OptProcessCompositeFactors.JOINTLY_MOD_Np,
       OptProcessCompositeFactors.SEPARATELY_MOD_Np]:
      # Let N' (denoted Np in the code) be the product of all pairwise coprime
      # composite factors of N stored in the collection.
      Np = F.residual;
    else:
      raise Exception("Error: Invalid option: opt_process_composite_factors.");

    while True:
      # Sample x uniformly at random from Z_N'^*.
      x = IntegerModRing(Np).random_element();
      if x == 0:
        continue; # Not in Z_N'^*, and not a non-trivial factor of N'.

      # Optimization: Sample x uniformly at random from Z_N'^* \ {1}.
      #
      # For details, see "optimizations.md" and Section 3.2.1 of [E21b].
      if (x == 1) and opt_exclude_one:
        print("Note: Sampled x = 1; excluding and sampling again...\n");
        continue;

      d = gcd(x.lift(), Np);
      if d == 1:
        break; # The element is in Z_N'^*.

      # Optimization: Report the non-trivial factor d found "by accident" above.
      #
      # For further details, see "optimizations.md".
      if opt_report_accidental_factors:
        print("Note: Reporting a factor (" + str(d) + ") found \"by " +
              "accident\" when sampling. This is likely to occur only if N " +
              "has small factors.\n");

        # Report the factor.
        F.add(d);

        # Print status again.
        F.print_status();

        # Check again if we are done, since we have updated F.
        if F.is_complete():
          break;

        if opt_process_composite_factors in \
          [OptProcessCompositeFactors.JOINTLY_MOD_Np,
           OptProcessCompositeFactors.SEPARATELY_MOD_Np]:
          Np = F.residual; # Update to the potentially new N'.

    # Check again if we are done, since we may have updated F above.
    if F.is_complete():
      break;

    # Optimization: Exponentiate x modulo N', for N' the product of all pairwise
    # coprime composite factors of N stored in the collection, or as N' runs
    # over the pairwise coprime composite factor of N stored in the collection.
    # This as opposed to exponentiating x modulo N', for N' = N, when not
    # applying the optimization.
    #
    # For further details, see "optimizations.md" and Section 3.2.1 of [E21b].
    if opt_process_composite_factors == \
      OptProcessCompositeFactors.SEPARATELY_MOD_Np:
      # Note: For SEPARATELY_MOD_Np, we compute the set of composite pairwise
      # coprime factors stored in the collection and let N' run over the set.
      factors = F.found_factors.difference(F.found_primes);
    elif opt_process_composite_factors in \
      [OptProcessCompositeFactors.JOINTLY_MOD_Np,
       OptProcessCompositeFactors.JOINTLY_MOD_N]:
      # Note: For JOINTLY_MOD_N we have N' = N (where we recall that N' is
      # denoted Np in the code). For JOINTLY_MOD_Np, we have that N' is the
      # product of all pairwise coprime composite factors of N stored in the
      # collection. (Note: This by the manner in which N' was setup above.)
      factors = set([Np]);
    else:
      raise Exception("Error: Invalid option: opt_process_composite_factors.");

    # Exponentiate x for each factor in the set setup above.
    #
    # Note that when opt_process_composite_factors is set to SEPARATELY_MOD_Np,
    # for each composite factor N' processed, any non-trivial factors of N'
    # reported can only split N', as N' is coprime with all other factors of N
    # stored in the collection. This implies that there is no need to go back
    # and re-examine the factor collection after it has been updated with the
    # non-trivial factors reported.
    for Np in factors:
      Rp = IntegerModRing(Np); # Define the subring Z_N'^*.
      xp = Rp(x); # Coerce x to Z_N'^*.

      # Step 4.2: For i = 0, 1, .., t do:
      timer_exponentiation.start();
      tmp = xp^o;
      timer_exponentiation.stop();

      # Optimization: Test if tmp = 1, and if so abort early.
      #
      # Were we to proceed, we would obtain d = gcd(tmp^{2^i} - 1, Np) = Np in
      # step 4.2.1, as tmp^{2^i} = 1 for all i, and as gcd(0, Np) = Np.
      #
      # For further details, see "optimizations.md" and Section 3.2.1 of [E21b].
      if (tmp == 1) and opt_abort_early:
        continue; # No point in continuing with the below procedure.

      # Step 4.2.1 for i = 0.
      d = gcd((tmp - 1).lift(), Np);
      if 1 < d < Np:
        F.add(d);

      for i in range(1, t + 1):
        # Optimization: To speed up the arithmetic, we may use a temporary
        # variable tmp, that we initially set to xp^o and then square
        # repeatedly, as opposed to computing xp^(2^i o) for each i.
        #
        # For more details, see "optimizations.md" and Section 3.2.1 of [E21b].
        timer_exponentiation.start();
        if opt_square:
          tmp = tmp^2;
        else:
          tmp = xp^((2^i) * o);
        timer_exponentiation.stop();

        if (tmp == 1) and opt_abort_early:
          break; # No point in continuing to iterate, see above.

        # Step 4.2.1 for i = 1, .., t.
        d = gcd((tmp - 1).lift(), Np);
        if 1 < d < Np:
          F.add(d);

  # Stop the timer.
  timer.stop();

  # The complete factorization has been found.
  print("Time required to solve:", timer);
  print(" Time spent exponentiating:", timer_exponentiation);
  print(" Time spent checking primality:", F.timer_test_primality);
  print(" Time spent reducing perfect powers:", F.timer_test_perfect_power);

  # Sanity check to assert that the factorization is correct and complete.
  tmp = N;

  for f in F.found_primes:
    if (tmp % f) != 0:
      raise Exception("Error: Failed to factor N correctly.");

    tmp /= f;
    while (tmp % f) == 0:
      tmp /= f;

  if tmp != 1:
    raise Exception("Error: Failed to completely factor N.");

  # Return the set of prime factors found.
  return F;
