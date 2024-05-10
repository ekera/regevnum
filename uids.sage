def uid():

    """
    @brief  Returns a random 32-bit UID.

    @return   The UID generated.
    """

    ID = IntegerModRing(2^32).random_element().lift().hex()
    while len(ID) < 8:
        ID = "0" + ID
    return ID


def uid_for_N(N):

    """
    @brief  Returns a 32-bit UID deterministically generated from an integer N.

    @param N  The integer N for which to generate a UID.

    @return   The UID generated from the integer N.
    """

    ID = (hash(N) % 2^32).hex()
    while len(ID) < 8:
        ID = "0" + ID
    return str(ZZ(N).nbits()) + "-" + ID