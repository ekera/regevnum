load("logarithm-finding.sage")


def main():
    if test_logarithm_finding_in_safe_prime_group(l=128, verbose=False):
      print("Passed: test_logarithm_finding_in_safe_prime_group()")
    else:
      print("Falied: test_logarithm_finding_in_safe_prime_group()")

    if test_logarithm_finding_in_schnorr_group(l=128, k=32, verbose=False):
       print("Passed: test_logarithm_finding_in_schnorr_group()")
    else:
      print("Failed: test_logarithm_finding_in_schnorr_group()")


    if test_logarithm_finding_in_safe_prime_group_with_precomputation(
          l=128,
          verbose=False):
        print("Passed: test_logarithm_finding_in_safe_prime_group_with_" + \
                  "precomputation()")
    else:
        print("Failed: test_logarithm_finding_in_safe_prime_group_with_" + \
                  "precomputation()")


    if test_logarithm_finding_in_schnorr_group_with_precomputation(
          l=128,
          k=32,
          verbose=False):
        print("Passed: test_logarithm_finding_in_schnorr_group_with_" + \
                  "precomputation()")
    else:
        print("Failed: test_logarithm_finding_in_schnorr_group_with_" + \
                  "precomputation()")


    C = find_minimum_C_for_logarithm_finding_in_safe_prime_group(
            l=128,
            verbose=False,
            log_file_prefix="tests/tmp/logarithm-")
    print("Info: find_minimum_C_for_logarithm_finding_in_safe_prime_" + \
              f"group() yielded C = {C}")

    C = find_minimum_C_for_logarithm_finding_in_safe_prime_group(
            l=128,
            failure_rate=0.3,
            verbose=False,
            log_file_prefix="tests/tmp/logarithm-failures-")
    print("Info: find_minimum_C_for_logarithm_finding_in_safe_prime_" + \
              f"group() yielded C = {C} when the failure rate = 0.3")


if __name__ == "__main__":
    main()
