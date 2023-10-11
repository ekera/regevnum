# Note: This file was copied from the factoritall repository by M. Ekerå.
# 
# More specifically, this file contains the class FactorCollection defined in
# the file factor.sage in the factoritall repository.
#
# For further details, please see: https://github.com/ekera/factoritall

from timer import Timer;

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