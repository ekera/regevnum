load("order-finding.sage")


def main():
    if test_order_finding(t=2, l=32, verbose=False):
        print("Passed: test_order_finding()")
    else:
        print("Failed: test_order_finding()")

    if test_order_finding_phi(t=2, l=32, verbose=False):
        print("Passed: test_order_finding_phi()")
    else:
        print("Failed: test_order_finding_phi()")


    C = find_minimum_C_for_order_finding(
            t=2,
            l=32,
            verbose=False,
            log_file_prefix="tests/tmp/order-finding-")
    print(f"Info: find_minimum_C_for_order_finding() yielded C = {C}")

    C = find_minimum_C_for_order_finding(
            t=2,
            l=32,
            failure_rate=0.3,
            verbose=False,
            log_file_prefix="tests/tmp/order-finding-failures-")
    print(f"Info: find_minimum_C_for_order_finding() yielded C = {C} " + \
              "when the failure rate = 0.3")

    C = find_minimum_C_for_order_finding_phi(
            t=2,
            l=32,
            verbose=False,
            log_file_prefix="tests/tmp/order-finding-phi-")
    print(f"Info: find_minimum_C_for_order_finding_phi() yielded C = {C}")

    C = find_minimum_C_for_order_finding_phi(
            t=2,
            l=32,
            failure_rate=0.3,
            verbose=False,
            log_file_prefix="tests/tmp/order-finding-phi-failures-")
    print(f"Info: find_minimum_C_for_order_finding_phi() yielded C = {C} " + \
              "when the failure rate = 0.3")


if __name__ == "__main__":
    main()
