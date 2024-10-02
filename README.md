# Simulating Regev's quantum factoring algorithm and Ekerå–Gärtner's extensions to discrete logarithm finding, order finding and factoring via order finding
The repository contains [Sage](https://www.sagemath.org) scripts that implement a simulator for the quantum part of Regev's multi-dimensional variation [[Regev23]](https://doi.org/10.48550/arXiv.2308.06572) of Shor's factoring algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), and of Ekerå–Gärtner's extensions [[EG23p]](https://doi.org/10.48550/arXiv.2311.05545) of Regev's algorithm to discrete logarithm finding, order finding, and factoring via order finding.
The scripts furthermore implement the lattice-based post-processing algorithms from the aforementioned works.

The simulator works by constructing a basis for the lattice for the problem instance to be simulated.
For the simulator to be able to construct such a basis, it requires the modulus that defines the problem instance to be on special form.
More specifically,
- when the modulus is a prime $p$, the simulator requires $p - 1$ to be smooth, and
- when the modulus is a composite $`N = {\prod}_{i=1}^t \, p_i`$, for the $p_i$ distinct prime factors, the simulator requires each $p_i - 1$ to be smooth and to not share any factor with $p_j - 1$ for $j \neq i$ except for a factor of two.

Imposing the above requirements enables the simulator to efficiently compute discrete logarithms in $`\mathbb Z_N^*`$ and $`\mathbb Z_p^*`$, respectively, which in turn enables it to construct a basis for the lattice.
Given this basis, a basis for the dual can be trivially computed, and noisy vectors sampled from the dual.

Note that imposing the above requirements makes the problem instance classically tractable:
The simulator can hence not be used to break classically hard problem instances.
This having been said, we expect its performance (in terms for the success probability of the post-processing for a given problem instance size and set of parameters) to be representative of that for Regev's algorithm, and of Ekerå–Gärtner's extensions thereof, also for classically hard problem instances.

For factoring, the simulator is sufficiently efficient to simulate Regev's algorithm, and our extensions of it to factoring via order finding, for 2048-bit RSA integers.
For discrete logarithms, the simulator is sufficiently efficient to simulate computing discrete logarithms in safe-prime groups and Schnorr groups, as used in Diffie–Hellman and DSA, with 2048-bit moduli.

The high-level functionality for factoring, logarithm finding, and order finding (including factoring via order finding), respectively, is implemented in the scripts
- [<code>factoring.sage</code>](factoring.sage),
- [<code>logarithm-finding.sage</code>](logarithm-finding.sage), and
- [<code>order-finding.sage</code>](order-finding.sage).

These scripts also contain convenience test functions, and functions for finding the minimum constant $C$ that suffices for the post-processing to be successful for a given problem instance and parameterization.

The aforementioned high-level scripts in turn depend on a number of other scripts, such as
- [<code>simulator.sage</code>](simulator.sage) that implements the simulator,
- [<code>problem-instances.sage</code>](problem-instances.sage) that samples problem instances of special form,
- [<code>parameter-search.sage</code>](parameter-search.sage) that implements the search for $C$, and
- [<code>common.sage</code>](simulator.sage) that defines default parameters.

To support the [<code>factoring.sage</code>](factoring.sage) and [<code>order-finding.sage</code>](order-finding.sage) scripts, this repository also contains a [Sage](https://www.sagemath.org) script [<code>dependencies/factor.sage</code>](dependencies/factor.sage), alongside supporting [Python](https://www.python.org) scripts, that are copied directly from the [factoritall](https://www.github.com/ekera/factoritall) repository.
These scripts implement procedures from [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).
They are used by to improve the post-processing for Regev's factoring algorithm, so as to use all available factoring relations to attempt to factor completely.
Furthermore, they are used to factor completely via order finding via Ekerå–Gärtner's extensions.

The aforementioned scripts were developed for academic research purposes.
They grew out of our research project in an organic manner as research questions were posed and answered.
They are distributed "as is" without warranty of any kind, either expressed or implied.
For further details, see the [license](LICENSE.md).

## Prerequisites
To install [Sage](https://www.sagemath.org) under [Ubuntu 22.04 LTS](https://releases.ubuntu.com/22.04), simply execute:

```console
$ sudo apt install sagemath
```
For other Linux and Unix distributions, or operating systems, you may need to [download Sage](https://www.sagemath.org/download) and install it manually.
These scripts were developed for Sage 9.5.

### Loading the scripts
Launch [Sage](https://www.sagemath.org) in this directory and load the main scripts by executing:

```console
$ sage
(..)
sage: load("factoring.sage")
sage: load("logarithm-finding.sage")
sage: load("order-finding.sage")
```

## Examples
For [examples](examples) that illustrate how to use the simulator, please see
- [<code>examples/factoring.md</code>](examples/factoring.md) for factoring,
- [<code>examples/logarithm-finding.md</code>](examples/logarithm-finding.md) for finding discrete logarithms, and
- [<code>examples/order-finding.md</code>](examples/order-finding.md) for finding orders, and for factoring completely via order finding.

## About and acknowledgments
These scripts were developed by [Martin Ekerå](mailto:ekera@kth.se) and Joel Gärtner, in part at [KTH, the Royal Institute of Technology](https://www.kth.se/en), in Stockholm, [Sweden](https://www.sweden.se).
Martin Ekerå thanks Sam Jaques for useful discussions.

Funding and support was provided by the Swedish NCSA that is a part of the [Swedish Armed Forces](https://www.mil.se).
