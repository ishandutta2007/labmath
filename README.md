labmath 2.2.0
=============

pip install labmath Copy PIP instructions

[Latest version](https://pypi.org/project/labmath/)

Released: Nov 27, 2021

Module for basic math in the general vicinity of computational number theory


### Navigation

*   [Project description](https://pypi.org/project/labmath/#description)
*   [Release history](https://pypi.org/project/labmath/#history)
*   [Download files](https://pypi.org/project/labmath/#files)

### Project links

*   [Homepage](https://pypi.org/manage/project/labmath)

### Statistics

View statistics for this project via [Libraries.io](https://libraries.io/pypi/labmath "External link"), or by using [our public dataset on Google BigQuery](https://packaging.python.org/guides/analyzing-pypi-package-downloads/)

### Meta

**License:** MIT License

**Author:** [lucasbrown.cit](mailto:lucasbrown.cit@gmail.com)

Tags math, mathematics, computational, number, theory, integer, factoring, factorization, primes, prime, numbers, semiprime, almost, prime, almost-prime, legendre, symbol, jacobi, symbol, kronecker, symbol, elliptic, curve, method, bpsw, miller, rabin, quadratic, frobenius, prp, sprp, lprp, slprp, xslprp, primality, testing, linear, recurrences, lucas, sequences, modular, square, root, generalized, Pell, equations, divisor, counting, function, euler's, totient, function, mobius, function, möbius, function, continued, fractions, partitions, stormer's, theorem, størmer's, theorem, smooth, numbers, Dirichlet, convolution

**Requires:** Python >=3.6


### Maintainers

 [![Avatar for lucasbrown.cit from gravatar.com](./labmath · PyPI_files/68747470733a2f2f7365637572652e67726176617461722e636f6d2f6176617461722f30333861633230343932316334323235656330333533653161323064383637363f73697a653d3530 "Avatar for lucasbrown.cit from gravatar.com")lucasbrown.cit](https://pypi.org/user/lucasbrown.cit/)

### Classifiers

*   **Development Status**
    *   [5 - Production/Stable](https://pypi.org/search/?c=Development+Status+%3A%3A+5+-+Production%2FStable)
*   **Environment**
    *   [Console](https://pypi.org/search/?c=Environment+%3A%3A+Console)
*   **Intended Audience**
    *   [Education](https://pypi.org/search/?c=Intended+Audience+%3A%3A+Education)
    *   [Science/Research](https://pypi.org/search/?c=Intended+Audience+%3A%3A+Science%2FResearch)
*   **License**
    *   [OSI Approved :: MIT License](https://pypi.org/search/?c=License+%3A%3A+OSI+Approved+%3A%3A+MIT+License)
*   **Operating System**
    *   [OS Independent](https://pypi.org/search/?c=Operating+System+%3A%3A+OS+Independent)
*   **Programming Language**
    *   [Python :: 3](https://pypi.org/search/?c=Programming+Language+%3A%3A+Python+%3A%3A+3)
*   **Topic**
    *   [Scientific/Engineering](https://pypi.org/search/?c=Topic+%3A%3A+Scientific%2FEngineering)
    *   [Scientific/Engineering :: Mathematics](https://pypi.org/search/?c=Topic+%3A%3A+Scientific%2FEngineering+%3A%3A+Mathematics)
    *   [Software Development :: Libraries](https://pypi.org/search/?c=Topic+%3A%3A+Software+Development+%3A%3A+Libraries)
    *   [Software Development :: Libraries :: Python Modules](https://pypi.org/search/?c=Topic+%3A%3A+Software+Development+%3A%3A+Libraries+%3A%3A+Python+Modules)

[![Sponsored: Python Software Foundation](./labmath · PyPI_files/Ansys_Logo_WhiteBackground_PNG_qjqAtuZ.png)](https://server.ethicalads.io/proxy/click/3341/92ddc541-2815-4c45-b302-bc56042e1abd/)

Ansys is a Maintaining sponsor of the Python Software Foundation.

_[PSF Sponsor](https://www.python.org/psf/sponsorship/?ref=ethicalads-placement) · [Served ethically](https://www.ethicalads.io/sponsorship-platform/?ref=psf)_

*   [Project description](https://pypi.org/project/labmath/#description)
*   [Project details](https://pypi.org/project/labmath/#data)
*   [Release history](https://pypi.org/project/labmath/#history)
*   [Download files](https://pypi.org/project/labmath/#files)



Project description
-------------------

labmath version 2.2.0
---------------------

This is a module for basic math in the general vicinity of computational number theory. It includes functions associated with primality testing, integer factoring, prime counting, linear recurrences, modular square roots, generalized Pell equations, the classic arithmetical functions, continued fractions, partitions, Størmer’s theorem, smooth numbers, and Dirichlet convolution. Integer arithmetic is used wherever feasible.

Functions
---------

We make a few imports:

    from multiprocessing import Process, Queue as mpQueue
    from itertools import chain, count, groupby, islice, tee, cycle, takewhile, compress, product, zip_longest
    from fractions import Fraction
    from random import randrange
    from math import log, log2, ceil, sqrt, factorial; inf = float('inf')
    from heapq import merge
    
    try: from gmpy2 import mpz; mpzv, inttypes = 2, (int, type(mpz(1)))
    except ImportError: mpz, mpzv, inttypes = int, 0, (int,)
    
    labmathversion = "2.2.0"

The _new_ functions provided by this module are as follows. Further details, including examples and input details, are available in docstrings and accessible via the built-in help function.

    primegen(limit=inf)

Generates primes less than the given limit (which may be infinite) lazily via a segmented sieve of Eratosthenes. Memory usage depends on the sequence of prime gaps; on Cramér’s conjecture, it is O(sqrt(_p_) · log(_p_)2).

    rpn(instr)

Evaluates a string in reverse Polish notation. The acceptable binary operators are \+ - \* / // % \*\* and correspond directly to those same operators in Python3 source code. The acceptable unary operators are ! #. They are the factorial and primorial, respectively. There are four aliases: x for \* , xx for \*\* , f for !, and p! for #.

    iterprod(l)

Product of the elements of any iterable. The product of an empty iterable is 1. DEPRECATED: now that Python 3.8 has math.prod, this function definition will in a future version be replaced with the statement from math import prod. That will probably happen when PyPy3 supports Python 3.8.

    listprod(l)

Product of the elements of a list. The product of the empty list is 1. We use a binary algorithm because this can easily generate huge numbers, and calling reduce(lambda x,y: x\*y, a) in such situations is quite a bit slower due to the time-complexity of multiplication. However, the size of the problem required to make this superior to iterprod() is quite large, so iterprod() should usually be used instead.

    polyval(f, x, m=None)

Evaluates a polynomial at a particular point, optionally modulo something.

    binomial(n,k)

The binomial coefficient nCr(n, k). DEPRECATED: now that Python 3.8 has math.comb, this function definition will in a future version be replaced with the statement from math import comb. That will probably happen when PyPy3 supports Python 3.8.

    powerset(l)

Generates the powerset of a list, tuple, or string. The yielded objects are always lists.

    primephi(x, a, ps, phicache={})

Legendre’s phi function. Helper function for primepi.

    primepi(x, ps=[], picache={}, phicache={}, sqrts={})

Computes the number of primes ≤ x via the Meissel-Lehmer method. The arguments ps, pichache, phicache, and sqrts are for internal use only.

    primesum(n)

Sum of primes ≤ n.

    altseriesaccel(a, n)

Convergence acceleration for alternating series. This is algorithm 1 from _Convergence Acceleration of Alternating Series_ by Cohen, Villegas, and Zagier [(pdf)](https://people.mpim-bonn.mpg.de/zagier/files/exp-math-9/fulltext.pdf), with a minor tweak so that the _d_\-value isn’t computed via floating point.

    riemannzeta(n, k=24)

Computes the Riemann zeta function by applying altseriesaccel to the [Dirichlet eta function](https://en.wikipedia.org/wiki/Dirichlet_eta_function). Should be rather accurate throughout the complex plane except near n such that 1 = 2n-1.

    zetam1(n, k=24)

Computes riemannzeta(n, k) - 1 by applying altseriesaccel to the Dirichlet eta function. Designed to be accurate even when riemannzeta(n) is machine-indistinguishable from 1.0 — in particular, when n is a large real number.

    riemannR(x, n=None, zc={})

Uses the [Gram series](http://mathworld.wolfram.com/GramSeries.html) to compute [Riemann’s R function](http://mathworld.wolfram.com/RiemannPrimeCountingFunction.html), which is a rather good approximation to primepi. The argument zc is a cache of zeta values.

    nthprimeapprox(n)

Produces an integer that should be rather close to the nth prime by using binary splitting on Riemann’s R function.

    nthprime(n)

Returns the nth prime (counting 2 as #1). This is done with some efficiency by using nthprimeapprox as an initial estimate, computing primepi of that, and then sieving to remove the error.

    gcd(a, *r)

Greatest common divisor of any number of values. Now that math.gcd supports any number of arguments, this function definition will in some future version be replaced with from math import gcd.

    xgcd(a, b)

Extended Euclidean altorithm: returns a tuple (g, _x_, _y_) such that g = gcd(a, b) and g = a·_x_ + b·_y_.

    modinv(a, m)

Returns the inverse of a modulo m, normalized to lie between 0 and m-1. If a is not coprime to m, returns 1. DEPRECATED: as of version 3.8, this can be computed using Python’s built-in pow function as pow(a, \-1, m). As such, a future version of this library will remove this function. That will probably happen once PyPy3 supports Python 3.8.

    crt(rems, mods)

Returns the unique integer _c_ in range(iterprod(mods)) such that _c_ ≡ _x_ (mod _y_) for (_x_, _y_) in zip(rems, mods). All elements of mods must be pairwise coprime.

    lcm(a, *r)

The least common multiple of any number of values. Now that math.lcm supports any number of arguments, a future version of this library will replace this function definition with from math import lcm.

    isqrt(n)

Greatest integer whose square is ≤ n. Now that Python 3.9 has the math.isqrt function, a future version of this library will remove this function definition in favor of the line from math import isqrt.

    introot(n, r=2)

For non-negative n, returns the greatest integer ≤ the rth root of n. For negative n, returns the least integer ≥ the rth root of n, or None if r is even.

    semiprimegen()

Generates the semiprimes. This is done by filtering the primes out of the output of pspgen.

    pspgen()

Generates the primes and semiprimes. This is done using a segmented sieve based on the sieve of Eratosthenes and the fact that these are precisely the numbers not divisible by any smaller semiprimes.

    almostprimegen(k)

Generates the k\-almost-primes, which are the numbers that have precisely k prime factors, counted with multiplicity. This is done by filtering nearlyprimegen(k-1) out of the output of nearlyprimegen(k).

    nearlyprimegen(k)

Generates the numbers (other than 1) that have k or fewer prime factors, counted with multipicity. This is done via a segmented sieve based on the sieve of Eratosthenes and the fact that these are precisely the numbers not divisible by any smaller k\-almost-primes.

    ispower(n, r=0)

If r = 0: If n is a perfect power, return a tuple containing the largest integer that, when squares/cubed/etc, yields n as the first component and the relevant power as the second component. If n is not a perfect power, return None.

If r > 0: We check whether n is a perfect rth power; we return its rth root if it is and None if it isn’t.

    ilog(x, b)

Greatest integer _k_ such that bk ≤ x.

    fibogen()

Generates the Fibonacci numbers, starting with 0 and 1.

    fibo(n, f={0:0, 1:1, 2:1})

Efficiently extracts the nth Fibonacci number, indexed so that fibo(0) = 0 and fibo(1) = fibo(2) = 1. The argument f is used for memoization. We compute O(log(n)) earlier Fibonaccis along the way. This is, in the big-O sense, just about as fast as possible.

    fibomod(n, m, f={0:0, 1:1, 2:1})

Efficiently extracts the nth Fibonacci number modulo m, indexed so that fibo(0) = 0 and fibo(1) == fibo(2) = 1. The argument f is used for memoization. We compute O(log(n)) earlier Fibonaccis along the way. This is the asymptotically fastest algorithm.

    lucaschain(n, x0, x1, op1, op2)

Algorithm 3.6.7 from _Prime Numbers: A Computational Perspective_ by Crandall & Pomerance (2nd edition): Evaluation of a binary Lucas chain. To quote their description:

> For a sequence _x_0, _x_1, … with a rule for computing _x_2j from _x_j and a rule for computing _x_2j+1 from _x_j and _x_j+1, this algorithm computes (_x_n, _x_n+1) for a given positive integer _n_. We have _n_ in binary as (_n_0, _n_1, …, _n_b-1) with _n_0 being the low-order bit. We write the rules as follows: _x_2j = op1(_x_j) and _x_2j+1 = op2(_x_j, _x_j+1).

    lucasgen(P, Q):

Generates the Lucas U- and V-sequences with parameters (P, Q).

    lucas(k, P, Q)

Efficiently computes the kth terms in the Lucas U- and V-sequences U(P, Q) and V(P, Q). More explicitly, if

> U0, U1, V0, V1 = 0, 1, 2, P

and we have the recursions

> Un = P · Un-1 - Q · Un-2
> 
> Vn = P · Vn-1 - Q · Vn-2

then we compute Uk and Vk in O(ln(k)) arithmetic operations. If P2 ≠ 4·Q, then these sequences grow exponentially, so the number of bit operations is anywhere from O(k · ln(k)2 · ln(ln(k))) to O(k2) depending on how multiplication is handled. We recommend using MPZs when k > 100 or so. We divide by P2 - 4·Q at the end, so we handle separately the case where this is zero.

    binlinrecgen(P, Q, a, b)

The general binary linear recursion. Exactly like lucasgen, except we only compute one sequence, and we supply the seeds.

    binlinrec(k, P, Q, a, b)

The general binary linear recursion. Exactly like lucas, except we compute only one sequence, and we supply the seeds.

    linrecgen(a, b, m=None)

The general homogenous linear recursion: we generate in order the sequence defined by

> _x_n+1 = ak · _x_n + ak-1 · _x_n-1 + … + a0 · _x_n-k,

where the initial values are \[_x_0, …, _x_k\] = b. If m is supplied, then we compute the sequence modulo m. The terms of this sequence usually grow exponentially, so computing a distant term incrementally by plucking it out of this generator takes O(n2) bit operations. Extractions of distant terms should therefore be done via linrec, which takes anywhere from O(n · ln(n)2 · ln(ln(n))) to O(n2) bit operations depending on how multiplication is handled.

    linrec(n, a, b, m=None)

The general homogeneous linear recursion. If m is supplied, terms are computed modulo m. We use matrix methods to efficiently compute the nth term of the recursion

> _x_n+1 = ak · _x_n + ak-1 · _x_n-1 + … + a0 · _x_n-k,

where the initial values are \[_x_0, …, _x_k\] = b.

    legendre(a, p)

Legendre symbol (a | p): 1 if a is a quadratic residue mod p, -1 if it isn’t, and 0 if a ≡ 0 (mod p). Not meaningful if p isn’t prime.

    jacobi(a, n)

The Jacobi symbol (a | n).

    kronecker(a, n)

The Kronecker symbol (a | n). Note that this is the generalization of the Jacobi symbol, _not_ the Dirac-delta analogue.

    sprp(n, b)

The strong probable primality test (aka single-round Miller-Rabin).

    mrab(n, basis)

Miller-Rabin probable primality test.

    miller(n)

Miller’s primality test. If the extended Riemann hypothesis (the one about Dirichlet L-functions) is true, then this test is deterministic.

    lprp(n, a, b)

Lucas probable primality test as described in _Prime Numbers: A Computational Perspective_ by Crandall & Pomerance (2nd edition).

    lucasmod(k, P, Q, m)

Efficiently computes the kth terms of Lucas U- and V-sequences modulo m with parameters (P, Q). Currently just a helper function for slprp and xslprp. Will be upgraded to full status when the case gcd(D,m)!=1 is handled properly.

    slprp(n, a, b)

Strong lucas probable primality test as described on Wikipedia. Its false positives are a strict subset of those for lprp with the same parameters.

    xslprp(n, a)

Extra strong Lucas probable primality test as described on Wikipedia. Its false positives are a strict subset of those for slprp (and therefore lprp) with parameters (a, 1).

    bpsw(n)

The Baille-Pomerance-Selfridge-Wagstaff probable primality test. Infinitely many false positives are conjectured to exist, but none are known, and the test is known to be deterministic below 264.

    qfprp(n, a, b)

Quadratic Frobenius probable primality test as described in _Prime Numbers: A Computational Perspective_ by Crandall & Pomerance (2nd edition).

    polyaddmodp(a, b, p)

Adds two polynomials and reduces their coefficients mod p. Polynomials are written as lists of integers with the constant terms first. If the high-degree coefficients are zero, those terms will be deleted from the answer so that the highest-degree term is nonzero. We assume that the inputs also satisfy this property. The zero polynomial is represented by the empty list. If one of the input polynomials is None, we return None.

    polysubmodp(a, b, p)

Subtracts the polynomial b from a and reduces their coefficients mod p. Polynomials are written as lists of integers with the constant terms first. If the high-degree coefficients are zero, those terms will be deleted from the answer so that the highest-degree term is nonzero. We assume that the inputs also satisfy this property. The zero polynomial is represented by the empty list. If one of the input polynomials is None, we return None.

    polymulmodp(a, b, p)

Multiplies the polynomials a and b and reduces their coefficients mod p. Polynomials are written as lists of integers with the constant terms first. If the high-degree coefficients are zero, those terms will be deleted from the answer so that the highest-degree term is nonzero. We assume that the inputs also satisfy this property. The zero polynomial is represented by the empty list. If one of the input polynomials is None, we return None.

    polydivmodmodp(a, b, p)

Divides the polynomial a by the polynomial b and returns the quotient and remainder. The coefficients are interpreted mod p. Polynomials are written as lists of integers with the constant terms first. If the high-degree coefficients are zero, those terms will be deleted from the answer so that the highest-degree term is nonzero. We assume that the inputs also satisfy this property. The zero polynomial is represented by the empty list. If one of the input polynomials is None, we return None. The result is not guaranteed to exist; in such cases we return (None, None).

    gcmd(f, g, p)

Computes the greatest common monic divisor of the polynomials f and g. The coefficients are interpreted mod p. Polynomials are written as lists of integers with the constant terms first. If the high-degree coefficients are zero, those terms will be deleted from the answer so that the highest-degree term is nonzero. We assume that the inputs also satisfy this property. The zero polynomial is represented by the empty list. If one of the input polynomials is None, or if both input polynomials are \[\], we return None. The result is not guaranteed to exist; in such cases, we return None. Coded after algorithm 2.2.1 from _Prime Numbers: A Computational Perspective_ by Crandall & Pomerance (2nd edition).

    polypowmodpmodpoly(a, e, p, f)

Computes the remainder when the polynomial a is raised to the eth power and reduced modulo f. The coefficients are interpreted mod p. Polynomials are written as lists of integers with the constant terms first. If the high-degree coefficients are zero, those terms will be deleted from the answer so that the highest-degree term is nonzero. We assume that the inputs also satisfy this property. The zero polynomial is represented by the empty list. If one of the input polynomials is None, or if f == \[\], we return None. The answer is not guaranteed to exist. In such cases, we return None.

    frobenius_prp(n, poly, strong=False)

Grantham’s general Frobenius probable primality test, in both the strong and weak versions, as described in [his paper introducing the test](https://doi.org/10.1090/S0025-5718-00-01197-2).

    isprime(n, tb=(3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59))

The workhorse primality test. It is a BPSW primality test variant: we use the strong Lucas PRP test and preface the computation with trial division for speed. No composites are known to pass the test, though it is suspected that infinitely many will do so. There are definitely no such errors below 264. This function is mainly a streamlined version of bpsw.

    isprimepower(n)

Determines whether n is of the form _p_e for some prime number _p_ and positive integer _e_, and returns (_p_, _e_) if so.

    isprime_mersenne(p)

The Lucas-Lehmer test. Deterministically and efficiently checks whether the Mersenne number 2p\-1 is prime.

    nextprime(n, primetest=isprime)

Smallest prime strictly greater than n.

    prevprime(n, primetest=isprime)

Largest prime strictly less than n, or None if no such prime exists.

    randprime(digits, base=10, primetest=isprime)

Returns a random prime with the specified number of digits when rendered in the specified base.

    randomfactored(n, primetest=isprime)

Efficiently generates an integer selected uniformly from the range \[1, n\] with its factorization. Uses Adam Kalai’s algorithm, which uses in the average case O(log(n)2) primality tests. When called with the default primality test, this then uses O(log(n)3) arithmetic operations, which in turn results in just over O(log(n)4) to O(log(n)5) bit operations, depending on how multiplication is handled.

    sqrtmod_prime(a, p)

Finds _x_ such that _x_2 ≡ a (mod p). We assume that p is a prime and a is a quadratic residue modulo p. If any of these conditions is false, then the return value is meaningless.

    cbrtmod_prime(a, p)

Returns in a sorted list all cube roots of a mod p. There are a bunch of easily-computed special formulae for various cases with p != 1 (mod 9); we do those first, and then if p == 1 (mod 9) we use Algorithm 4.2 in [Taking Cube Roots in Zm](https://doi.org/10.1016/S0893-9659(02)00031-9) by Padro and Saez, which is essentially a variation on the Tonelli-Shanks algorithm for modular square roots.

    pollardrho_brent(n)

Factors integers using Brent’s variation of Pollard’s rho algorithm. If n is prime, we immediately return n; if not, we keep chugging until a nontrivial factor is found.

    pollard_pm1(n, B1=100, B2=1000)

Integer factoring function. Uses Pollard’s p-1 algorithm. Note that this is only efficient if the number to be factored has a prime factor _p_ such that _p_\-1’s largest prime factor is “small”.

    mlucas(v, a, n)

Helper function for williams\_pp1. Multiplies along a Lucas sequence modulo n.

    williams_pp1(n)

Integer factoring function. Uses Williams’ p+1 algorithm, single-stage variant. Note that this is only efficient when the number to be factored has a prime factor _p_ such that _p_+1’s largest prime factor is “small”.

    ecadd(p1, p2, p0, n)

Helper function for ecm. Adds two points on a Montgomery curve modulo n.

    ecdub(p, A, n)

Helper function for ecm. Doubles a point on a Montgomery curve modulo n.

    ecmul(m, p, A, n)

Helper function for ecm. Multiplies a point on Montgomery curve by an integer modulo n.

    secm(n, B1, B2, seed)

Seeded elliptic curve factoring using the two-phase algorithm on Montgomery curves. Helper function for ecm. Returns a possibly-trivial divisor of n given two bounds and a seed.

    ecmparams(n)

Generator of parameters to use for secm.

    ecm(n, paramseq=ecmparams, nprocs=1)

Integer factoring via elliptic curves using the two-phase algorithm on Montgomery curves, and optionally uses multiple processes. This is a shell function that repeatedly calls secm using parameters provided by ecmparams; the actual factoring work is done there. Multiprocessing incurs relatively significant overhead, so when nprocs==1 (default), we don’t call the multiprocessing functions.

    siqs(n)

Factors an integer via the self-initializing quadratic sieve. Most of this function is copied verbatim from [https://github.com/skollmann/PyFactorise](https://github.com/skollmann/PyFactorise).

    multifactor(n, methods)

Integer factoring function. Uses several methods in parallel. Waits for one of them to return, kills the rest, and reports.

    primefac(n, trial=1000, rho=42000, primetest=isprime, methods=(pollardrho_brent,))

The workhorse integer factorizer. Generates the prime factors of the input. Factors that appear _x_ times are yielded _x_ times.

    factorint(n, trial=1000, rho=42000, primetest=isprime, methods=(pollardrho_brent,))

Compiles the output of primefac into a dictionary with primes as keys and multiplicities as values.

    factorsieve(stop)

Uses a sieve to compute the factorizations of all whole numbers strictly less than the input. This uses a lot of memory; if you aren’t after the factors directly, it’s usually better to write a dedicated function for whatever it is that you actually want.

    divisors(n)

Generates all natural numbers that evenly divide n. The output is not necessarily sorted.

    divisors_factored(n)

Generates the divisors of n, written as their prime factorizations in factorint format.

    divcount(n)

Counts the number of divisors of n.

    divsigma(n, x=1)

Sum of divisors of a natural number, raised to the _x_th power. The conventional notation for this in mathematical literature is σx(n), hence the name of this function.

    divcountsieve(stop)

Uses a sieve to compute the number of divisors of all whole numbers strictly less than the input.

    totient(n, k=1)

Jordan’s totient function: the number of k\-tuples of positive integers all ≤ n that form a coprime (k+1)-tuple together with n. When k = 1, this is Euler’s totient: the number of numbers less than a number that are relatively prime to that number.

    totientsieve(n)

Uses a sieve to compute the totients up to (and including) n.

    totientsum(n)

Computes sum(totient(n) for n in range(1, n+1)) efficiently.

    mobius(n)

The Möbius function of n: 1 if n is squarefree with an even number of prime factors, -1 if n is squarefree with an odd number of prime factors, and 0 if n has a repeated prime factor.

    mobiussieve(stop)

Uses a sieve to compute the Möbius function of all whole numbers strictly less than the input.

    liouville(n)

The Liouville lambda function of n: the strongly multiplicative function that is -1 on the primes.

    polyroots_prime(g, p, sqfr=False)

Generates with some efficiency and without multiplicity the zeros of a polynomial modulo a prime. Coded after algorithm 2.3.10 from _Prime Numbers: A Computational Perspective_ by Crandall & Pomerance (2nd edition), which is essentially Cantor-Zassenhaus.

    hensel(f, p, k, given=None)

Uses Hensel lifting to generate with some efficiency all zeros of a polynomial modulo a prime power.

    sqrtmod(a, n)

Computes all square roots of a modulo n and returns them in a sorted list.

    polyrootsmod(pol, n)

Computes the zeros of a polynomial modulo an integer. We do this by factoring the modulus, solving modulo the prime power factors, and putting the results back together via the Chinese Remainder Theorem.

    PQa(P, Q, D)

Generates some sequences related to simple continued fractions of certain quadratic surds. A helper function for pell. Let P, Q, and D be integers such that Q ≠ 0, D > 0 is a nonsquare, and P2 ≡ D (mod Q). We yield a sequence of tuples (_B_i, _G_i, _P_i, _Q_i) where _i_ is an index counting up from 0, _x_ = (P+√D)/Q = \[_a_0; _a_1, _a_2, …\], (_P_i+√D))/_Q_i is the _i_th complete quotient of _x_, and _B_i is the denominator of the _i_th convergent to _x_. For full details, see [https://web.archive.org/web/20180831180333/http://www.jpr2718.org/pell.pdf](https://web.archive.org/web/20180831180333/http://www.jpr2718.org/pell.pdf).

    pell(D, N)

This function solves the generalized Pell equation: we find all non-negative integers (_x_, _y_) such that _x_2 - D · _y_2 = N. We have several cases:

Case 1: N = 0. We solve _x_2 = D · _y_2. (0,0) is always a solution.

> Case 1a: If D is a nonsquare, then there are no further solutions.
> 
> Case 1b: If D is a square, then there are infinitely many solutions, parametrized by (_t_·√D, _t_).

Case 2: N ≠ 0 = D. We solve _x_2 = N.

> Case 2a: If N is a nonsquare, then there are no solutions.
> 
> Case 2b: If N is a square, then there are infinitely many solutions, parametrized by (√N, _t_).

Case 3: N ≠ 0 > D. We solve _x_2 + |D| · _y_2 = N. The number of solutions will be finite.

Case 4: N ≠ 0 < D. We find lattice points on a hyperbola.

> Case 4a: If D is a square, then the number of solutions will be at most finite. This case is solved by factoring.
> 
> Case 4b: If D is a nonsquare, then we run the PQa/LMM algorithms: we produce a set of primitive solutions; if this set is empty, there are no solutions; if this set has members, an ininite set of solutions can be produced by repeatedly composing them with the fundamental solution of _x_2 - D · _y_2 = 1.

References:

*   [https://web.archive.org/web/20180831180333/https://www.jpr2718.org/pell.pdf](https://web.archive.org/web/20180831180333/https://www.jpr2718.org/pell.pdf)
    
*   [http://www.offtonic.com/blog/?p=12](http://www.offtonic.com/blog/?p=12)
    
*   [http://www.offtonic.com/blog/?p=18](http://www.offtonic.com/blog/?p=18)
    

Input: D, N – integers

Output:

> A 3-tuple.
> 
> If the number of solutions is finite, it is (None, z, None), where z is the sorted list of all solutions.
> 
> If the number of solutions is infinite and the equation is degenerate, it’s (gen, None, None), where gen yields all solutions.
> 
> If the number of solutions if infinite and the equation is nondegenerate, it is (gen, z, f), where z is the set of primitive solutions, represented as a sorted list, and f is the fundamental solution — i.e., f is the primitive solution of _x_2 - D · _y_2 = 1.
> 
> Note that we can check the infinitude of solutions by calling bool(pell(D,N)\[0\]).

    simplepell(D, bail=inf)

Generates the positive solutions of _x_2 - D · _y_2 = 1. We use some optimizations specific to this case of the Pell equation that makes this more efficient than calling pell(D,1)\[0\]. Note that this function is not equivalent to calling pell(D,1)\[0\]: pell is concerned with the general equation, which may or may not have trivial solutions, and as such yields all non-negative solutions, whereas this function is concerned only with the simple Pell equation, which always has an infinite family of positive solutions generated from a single primitive solution and always has the trivial solution (1,0).

We yield only those solutions with _x_ ≤ bail.

    carmichael(n)

The Carmichael lambda function: the smallest positive integer _m_ such that _a_m ≡ 1 (mod n) for all _a_ such that gcd(_a_, n) = 1. Also called the reduced totient or least universal exponent.

    multord(b, n)

Computes the multiplicative order of b modulo n; i.e., finds the smallest _k_ such that bk ≡ 1 (mod n).

    pythags_by_perimeter(p)

Generates all Pythagorean triples of a given perimeter by examining the perimeter’s factors.

    collatz(n)

Generates the Collatz sequence initiated by n. Stops after yielding 1.

    sqrtcfrac(n)

Computes the simple continued fraction for √n. We return the answer as (isqrt(n), \[a,b,c,...,d\]), where \[a,b,c,...,d\] is the minimal reptend.

    convergents(a)

Generates the convergents of a simple continued fraction.

    contfrac_rat(n, d)

Returns the simple continued fraction of the rational number n/d.

    quadratic_scf(P,Q,D)

Computes the simple continued fraction of the expression (P + sqrt(D)) / Q, for any integers P, Q, and D with D ≥ 0 ≠ Q. Note that D can be a square or a nonsquare.

    ngonal(x, n)

Returns the xth n\-gonal number. Indexing begins with 1 so that ngonal(1, n) = 1 for all applicable n.

    is_ngonal(p, n)

Checks whether p is an n\-gonal number.

    partitions(n, parts=[1])

Computes with some semblance of efficiency the number of additive partitions of an integer. The parts argument is for memoization.

    partgen(n)

Generates partitions of integers in ascending order via an iterative algorithm. It is the fastest known algorithm as of June 2014.

    partconj(p)

Computes the conjugate of a partition.

    farey(n)

Generates the Farey sequence of maximum denominator n. Includes 0/1 and 1/1.

    fareyneighbors(n, p, q)

Returns the neighbors of p/q in the Farey sequence of maximum denominator n.

    ispractical(n)

Tests whether n is a practical number – i.e., whether every integer from 1 through n (inclusive) can be written as a sum of divisors of n. These are also called panarithmic numbers.

    hamming(ps, *ps2)

Generates all ps\-smooth numbers, where ps is a list of primes.

    arithmeticderivative(n)

The arithmetic derivative of n: if n is prime, then n’ = 1; if -2 < n < 2, then n’ = 0; if n < 0, then n’ = -(-n)’; and (_ab_)’ = _a_’·_b_ + _b_’·_a_.

    perfectpowers()

Generates the sequence of perfect powers without multiplicity.

    sqfrgen(ps)

Generates the squarefree products of the elements of ps.

    sqfrgenb(ps, b, k=0, m=1)

Generates the squarefree products of elements of ps. Does not yield anything > b. For best performance, ps should be sorted in decreasing order.

    stormer(ps, *ps2, abc=None)

Størmer’s theorem asserts that for any given set ps of prime numbers, there are only finitely many pairs of consecutive integers that are both ps\-smooth; the theorem also gives an effective algorithm for finding them. We implement Lenstra’s improvement to this theorem.

The abc argument indicates that we are to assume an effective abc conjecture of the form _c_ < abc\[0\] · rad(_a_·_b_·_c_)abc\[1\]. This enables major speedups. If abc is None, then we make no such assumptions.

    quadintroots(a, b, c)

Given integers a, b, and c, we return in a tuple all distinct integers _x_ such that a·_x_2 + b·_x_ + c = 0. This is primarily a helper function for cubicintrootsgiven and cubicintroots.

    cubicintrootsgiven(a, b, c, d, r)

Given integers a, b, c, d, and r such that a·r3 + b·r2 + c·r + d = 0, we find the cubic’s other two roots and return in a tuple all distinct integer roots (including r). This is primarily a helper function for cubicintroots.

    cubicintroots(a, b, c, d)

Given integers a, b, c, d, we return in a tuple all distinct integer roots of a·_x_3 + b·_x_2 + c·_x_ + d. This is primarily a helper function for isprime\_nm1.

    isprime_nm1(n, fac=None)

The _n_\-1 primality test: given an odd integer n > 214 and a fully-factored integer _F_ such that _F_ divides n\-1 and _F_ > n0.3, we quickly determine without error whether n is prime. If the provided (partial) factorization of n\-1 is insufficient, we compute the factorization ourselves.

    isprime_np1(n, fac=None)

The _n_+1 primality test: given an odd integer n > 214 and a fully-factored integer _F_ such that _F_ divides n+1 and _F_ > n0.3, we quickly determine without error whether n is prime. If the provided (partial) factorization of n+1 is insufficient, we compute the factorization ourselves.

    mulparts(n, r=None, nfac=None)

Generates all ordered r\-tuples of positive integers whose product is n. If r is None, then we generate all such tuples (regardless of size) that do not contain 1.

    dirconv(f, g, ffac=False, gfac=False)

This returns a function that is the Dirichlet convolution of f and g. When called with the keyword arguments at their default values, this is equivalent to the expression lambda n: sum(f(d) \* g(n//d) for d in divisors(n)). If f or g needs to factor its argument, such as f == totient or g == mobius or something like that, then that lambda expression calls the factorizer a lot more than it needs to — we’re already factoring n, so instead of feeding those functions the integer forms of n’s factors, we can instead pass ffac=True or gfac=True when dirconv is called and we will call divisors\_factored(n) instead and feed those factored divisors into f or g as appropriate. This optimization becomes more noticeable as the factoring becomes more difficult.

    dirichletinverse(f)

Computes the Dirichlet inverse of the input function f. Mathematically, functions _f_ such that _f_(1) = 0 have no Dirichlet inverses due to a division by zero. This is reflected in this implementation by raising a ZeroDivisionError when attempting to evaluate dirichletinverse(f)(n) for any such f and any n. If f(1) is neither 1 nor -1, then dirichletinverse(f) will return Fraction objects (as imported from the fractions module).

    dirichletroot(f, r, val1)

Computes the rth Dirichlet root of the input function f whose value at 1 is val1. More precisely, let f be a function on the positive integers, let r be a positive integer, and let val1r = f(1). Then we return the unique function g such that f = g \* g \* … \* g, where g appears r times and \* represents Dirichlet convolution. The values returned will be Fraction objects (as imported from the fractions module).

    determinant(M)

Computes the determinant of a matrix via the Schur determinant identity.

    discriminant(coefs)

Computes the discriminant of a polynomial. The input list is ordered from lowest degree to highest — i.e., coefs\[k\] is the coefficient of the _x_k term. For low-degree polynomials, explicit formulae are used; for degrees 5 and higher, we compute it by taking the determinant (using this package’s determinant function) of the Sylvester matrix of the input and its derivative. This in turn is calculated by the Schur determinant identity. Note that this has the effect of setting the discriminant of a linear polynomial to 1 (which is conventional) and that of a constant to 0.

    egypt_short(n, d, terms=0, minden=1)

Generates all shortest Egyptian fractions for n/d using at least the indicated number of terms and whose denominators are all ≥ minden. No algorithm is known for this problem that significantly improves upon brute force, so this can take impractically long times on even modest-seeming inputs.

    egypt_greedy(n, d)

The greedy algorithm for Egyptian fraction expansion; also called the Fibonacci-Sylvester algorithm.

### Dependencies

This package imports items from multiprocessing, itertools, fractions, random, math, and heapq. These are all in the Python standard library.

We attempt to import mpz from gmpy2, but this is purely for efficiency: if this import fails, we simply set mpz = int.

Project details
---------------

### Project links

*   [Homepage](https://pypi.org/manage/project/labmath)

### Statistics

View statistics for this project via [Libraries.io](https://libraries.io/pypi/labmath "External link"), or by using [our public dataset on Google BigQuery](https://packaging.python.org/guides/analyzing-pypi-package-downloads/)

### Meta

**License:** MIT License

**Author:** [lucasbrown.cit](mailto:lucasbrown.cit@gmail.com)

Tags math, mathematics, computational, number, theory, integer, factoring, factorization, primes, prime, numbers, semiprime, almost, prime, almost-prime, legendre, symbol, jacobi, symbol, kronecker, symbol, elliptic, curve, method, bpsw, miller, rabin, quadratic, frobenius, prp, sprp, lprp, slprp, xslprp, primality, testing, linear, recurrences, lucas, sequences, modular, square, root, generalized, Pell, equations, divisor, counting, function, euler's, totient, function, mobius, function, möbius, function, continued, fractions, partitions, stormer's, theorem, størmer's, theorem, smooth, numbers, Dirichlet, convolution

**Requires:** Python >=3.6

### Maintainers

 [![Avatar for lucasbrown.cit from gravatar.com](./labmath · PyPI_files/68747470733a2f2f7365637572652e67726176617461722e636f6d2f6176617461722f30333861633230343932316334323235656330333533653161323064383637363f73697a653d3530 "Avatar for lucasbrown.cit from gravatar.com")lucasbrown.cit](https://pypi.org/user/lucasbrown.cit/)

### Classifiers

*   **Development Status**
    *   [5 - Production/Stable](https://pypi.org/search/?c=Development+Status+%3A%3A+5+-+Production%2FStable)
*   **Environment**
    *   [Console](https://pypi.org/search/?c=Environment+%3A%3A+Console)
*   **Intended Audience**
    *   [Education](https://pypi.org/search/?c=Intended+Audience+%3A%3A+Education)
    *   [Science/Research](https://pypi.org/search/?c=Intended+Audience+%3A%3A+Science%2FResearch)
*   **License**
    *   [OSI Approved :: MIT License](https://pypi.org/search/?c=License+%3A%3A+OSI+Approved+%3A%3A+MIT+License)
*   **Operating System**
    *   [OS Independent](https://pypi.org/search/?c=Operating+System+%3A%3A+OS+Independent)
*   **Programming Language**
    *   [Python :: 3](https://pypi.org/search/?c=Programming+Language+%3A%3A+Python+%3A%3A+3)
*   **Topic**
    *   [Scientific/Engineering](https://pypi.org/search/?c=Topic+%3A%3A+Scientific%2FEngineering)
    *   [Scientific/Engineering :: Mathematics](https://pypi.org/search/?c=Topic+%3A%3A+Scientific%2FEngineering+%3A%3A+Mathematics)
    *   [Software Development :: Libraries](https://pypi.org/search/?c=Topic+%3A%3A+Software+Development+%3A%3A+Libraries)
    *   [Software Development :: Libraries :: Python Modules](https://pypi.org/search/?c=Topic+%3A%3A+Software+Development+%3A%3A+Libraries+%3A%3A+Python+Modules)

  

Release history [Release notifications](https://pypi.org/help/#project-release-notifications) | [RSS feed](https://pypi.org/rss/project/labmath/releases.xml)
-------------------------------------------------------------------------------------------------------------------------------------------------------------

This version

![](./labmath · PyPI_files/blue-cube.fff94ebc.svg)

[

2.2.0

Nov 27, 2021

](https://pypi.org/project/labmath/2.2.0/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.1.3

Nov 14, 2021

](https://pypi.org/project/labmath/2.1.3/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.1.2

Nov 14, 2021

](https://pypi.org/project/labmath/2.1.2/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.1.1

Nov 13, 2021

](https://pypi.org/project/labmath/2.1.1/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.1.0

Nov 13, 2021

](https://pypi.org/project/labmath/2.1.0/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.0.2

Nov 10, 2021

](https://pypi.org/project/labmath/2.0.2/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.0.1

Nov 10, 2021

](https://pypi.org/project/labmath/2.0.1/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

2.0.0

Nov 10, 2021

](https://pypi.org/project/labmath/2.0.0/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

1.2.0

Nov 19, 2018

](https://pypi.org/project/labmath/1.2.0/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

1.1.0

Apr 28, 2018

](https://pypi.org/project/labmath/1.1.0/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

1.0.8

Mar 29, 2018

](https://pypi.org/project/labmath/1.0.8/)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

[

1.0.7

Mar 2, 2018

](https://pypi.org/project/labmath/1.0.7/)

Download files
--------------

Download the file for your platform. If you're not sure which to choose, learn more about [installing packages](https://packaging.python.org/installing/ "External link").

### Source Distribution

[labmath-2.2.0.tar.gz](https://files.pythonhosted.org/packages/76/02/93aaa577f88aae2c4d4b25bfd482928b6bf601bdcb9272fd3a5d959e4104/labmath-2.2.0.tar.gz) (94.4 kB [view hashes](https://pypi.org/project/labmath/#copy-hash-modal-901ec201-a223-470e-8d40-d2d69294e8d3))

Uploaded Nov 27, 2021 `source`

### Built Distribution

[labmath-2.2.0-py3-none-any.whl](https://files.pythonhosted.org/packages/14/da/cc4e782ac107401fb47b7399d0ba3277a34045948e07b8ca2760b825552d/labmath-2.2.0-py3-none-any.whl) (69.8 kB [view hashes](https://pypi.org/project/labmath/#copy-hash-modal-3e534cb4-c4d2-4b00-86e0-1a45caf6f4ff))

Uploaded Nov 27, 2021 `py3`

[Close](https://pypi.org/project/labmath/#modal-close "Close")

### [Hashes](https://pip.pypa.io/en/stable/cli/pip_install/#hash-checking-mode "External link") for labmath-2.2.0.tar.gz

Hashes for labmath-2.2.0.tar.gz

Algorithm

Hash digest

SHA256

`773278b333f19dc9340a0739204a79141987bc8c800b4aea29156b92ddb49ed4`

Copy

MD5

`0ba90b84fcd30b6db73a773adb9b58db`

Copy

BLAKE2b-256

`760293aaa577f88aae2c4d4b25bfd482928b6bf601bdcb9272fd3a5d959e4104`

Copy

[Close](https://pypi.org/project/labmath/#modal-close)

[Close](https://pypi.org/project/labmath/#modal-close "Close")

### [Hashes](https://pip.pypa.io/en/stable/cli/pip_install/#hash-checking-mode "External link") for labmath-2.2.0-py3-none-any.whl

Hashes for labmath-2.2.0-py3-none-any.whl

Algorithm

Hash digest

SHA256

`ac00538bab16346a7a56406022d9a00a899cc1ee2deb5b882d88e55c4f85345b`

Copy

MD5

`1996daaee8c1c78eeb4c18e6eaeec74f`

Copy

BLAKE2b-256

`14dacc4e782ac107401fb47b7399d0ba3277a34045948e07b8ca2760b825552d`

Copy

[Close](https://pypi.org/project/labmath/#modal-close)

![](./labmath · PyPI_files/white-cube.7e91782c.svg)

Help
----

*   [Installing packages](https://packaging.python.org/installing/ "External link")
*   [Uploading packages](https://packaging.python.org/tutorials/packaging-projects/ "External link")
*   [User guide](https://packaging.python.org/ "External link")
*   [Project name retention](https://www.python.org/dev/peps/pep-0541/ "External link")
*   [FAQs](https://pypi.org/help/)

About PyPI
----------

*   [PyPI on Twitter](https://twitter.com/PyPI "External link")
*   [Infrastructure dashboard](https://dtdg.co/pypi "External link")
*   [Statistics](https://pypi.org/stats/)
*   [Logos & trademarks](https://pypi.org/trademarks/)
*   [Our sponsors](https://pypi.org/sponsors/)

Contributing to PyPI
--------------------

*   [Bugs and feedback](https://pypi.org/help/#feedback)
*   [Contribute on GitHub](https://github.com/pypi/warehouse "External link")
*   [Translate PyPI](https://hosted.weblate.org/projects/pypa/warehouse/ "External link")
*   [Sponsor PyPI](https://pypi.org/sponsors/)
*   [Development credits](https://github.com/pypi/warehouse/graphs/contributors "External link")

Using PyPI
----------

*   [Code of conduct](https://github.com/pypa/.github/blob/main/CODE_OF_CONDUCT.md "External link")
*   [Report security issue](https://pypi.org/security/)
*   [Privacy policy](https://www.python.org/privacy/ "External link")
*   [Terms of use](https://pypi.org/policy/terms-of-use/)
*   [Acceptable Use Policy](https://pypi.org/policy/acceptable-use-policy/)

* * *

Status: [Service Under Maintenance](https://status.python.org/ "External link")

Developed and maintained by the Python community, for the Python community.  
[Donate today!](https://donate.pypi.org/)

"PyPI", "Python Package Index", and the blocks logos are registered [trademarks](https://pypi.org/trademarks/) of the [Python Software Foundation](https://python.org/psf-landing).  

© 2023 [Python Software Foundation](https://www.python.org/psf-landing/ "External link")  
[Site map](https://pypi.org/sitemap/)

Switch to desktop version

*   English
*   español
*   français
*   日本語
*   português (Brasil)
*   українська
*   Ελληνικά
*   Deutsch
*   中文 (简体)
*   中文 (繁體)
*   русский
*   עברית
*   esperanto
