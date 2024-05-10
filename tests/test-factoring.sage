load("factoring.sage")


def main():
    result = test_factoring(t=2, l=32, verbose=False)
    print("Info: test_factoring() with 2 factors yielded: " + \
              f"{result.is_complete()}")

    result = test_factoring(t=6, l=32, verbose=False)
    print("Info: test_factoring() with 6 factors yielded: " + \
              f"{result.is_complete()}")


    C = find_minimum_C_for_factoring(
            t=2,
            l=32,
            verbose=False,
            log_file_prefix="tests/tmp/factoring-2-factors-")
    print(f"Info: find_minimum_C_for_factoring() for 2 factors yielded C = {C}")

    C = find_minimum_C_for_factoring(
            t=2,
            l=32,
            failure_rate=0.3,
            verbose=False,
            log_file_prefix="tests/tmp/factoring-2-factors-failures-")
    print("Info: find_minimum_C_for_factoring() for 2 factors yielded "
              f"C = {C} when the failure rate = 0.3")

    C = find_minimum_C_for_factoring(
            t=6,
            l=32,
            verbose=False,
            log_file_prefix="tests/tmp/factoring-6-factors-")
    print(f"Info: find_minimum_C_for_factoring() for 6 factors yielded C = {C}")

    C = find_minimum_C_for_factoring(
            t=6,
            l=32,
            failure_rate=0.3,
            verbose=False,
            log_file_prefix="tests/tmp/factoring-6-factors-failures-")
    print("Info: find_minimum_C_for_factoring() for 6 factors yielded "
              f"C = {C} when the failure rate = 0.3")


if __name__ == '__main__':
    main()
