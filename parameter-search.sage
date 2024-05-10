from datetime import datetime

def find_minimum_C(
    test_for_C,
    precision=DEFAULT_PRECISION,
    log_file=None,
    verbose=True):

    """
    @brief  A function to search for the minimum value of the constant C for
            which the test function test_for_C() returns True.

    The search is essentially a binary search in C up to the specified level
    of precision. The search starts from C = 1.

    @param test_for_C   Function which tests if a specific value for C is
                        acceptable or not. Function is expected to return
                        a boolean value signifying the result.

    @param precision  The level of precision with which to perform the search.

    @param log_file   A log file to which to log intermediary results. If
                      omitted, no intermediary results are logged.

    @param verbose  A flag that may be set to True to print verbose status
                    messages, or to False not to print such messages.

    @return   The minimum value of the constant C.
    """

    # Setup parameters for the binary search.
    C = 1
    low = 0
    high = None

    step = round(C / 2 / precision)
    val = int(C / precision)

    while high == None or high - low > 1:
        C = val * precision

        if verbose:
            print("")
            print("** Attempting with C =", C)
            print("")

        success = test_for_C(C)

        if success:
            if verbose:
                print("")
                print("** Successful with C =", C)

            if None != log_file:
                print("Successful with C =", C, file=log_file)
                log_file.flush()

            high = val
            step = ceil((high - low) / 2)
            val -= step
        else:
            if verbose:
                print("")
                print("** Failed with C =", C)

            if None != log_file:
                print("Failed with C =", C, file=log_file)
                log_file.flush()

            low = val
            if high != None:
                step = ceil((high - low) / 2)
            else:
                step = ceil(low / 2)
            val += step

    if verbose:
        print("")
        print("** Picked C =", high * precision)
        print("")

    if None != log_file:
        print("\nPicked C =", high * precision, file=log_file)

        print(
            "\nStop time:",
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            file=log_file,
        )

        log_file.close()

    return high * precision
