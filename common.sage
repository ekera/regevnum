DEFAULT_THREADS=8
DEFAULT_BOUND_B=10^5
DEFAULT_CONSTANT_C=2
DEFAULT_PRECISION=0.1
DEFAULT_BLOCK_SIZE=2

ACTUAL_SAFE_PRIMES=False

def generate_new_random_seed():

  """
  @brief  Generates a new random 128-bit seed represented as an integer.

  @return   The random 128-bit seed generated.
  """

  return IntegerModRing(2^128).random_element().lift()