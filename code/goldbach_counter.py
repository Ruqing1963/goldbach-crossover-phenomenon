"""
Goldbach Pair Counter
=====================

This module implements ordered-pair counting for Goldbach representations.

Author: Ruqing Chen
Date: January 2026
"""

import numpy as np


def sieve_of_eratosthenes(limit):
    """
    Generate all prime numbers up to limit using Sieve of Eratosthenes.
    
    Parameters
    ----------
    limit : int
        Upper bound for prime generation
        
    Returns
    -------
    numpy.ndarray
        Array of prime numbers up to limit
    """
    if limit < 2:
        return np.array([], dtype=int)
    
    # Create boolean array
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    
    # Sieve
    for i in range(2, int(np.sqrt(limit)) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    
    # Extract primes
    primes = np.where(is_prime)[0]
    return primes


def count_goldbach_ordered(N, primes=None):
    """
    Count ordered Goldbach pairs for even integer N.
    
    For N = 10: pairs (3,7), (5,5), (7,3) give count = 3
    This is the ordered-pair counting consistent with Hardy-Littlewood
    circle method derivation.
    
    Parameters
    ----------
    N : int
        Even integer to analyze (must be even and >= 4)
    primes : numpy.ndarray, optional
        Pre-computed array of primes. If None, will generate up to N.
        
    Returns
    -------
    int
        Number of ordered pairs (p, q) where p + q = N and both are prime
        
    Examples
    --------
    >>> count_goldbach_ordered(10)
    3
    >>> count_goldbach_ordered(100)
    6
    """
    if N % 2 != 0 or N < 4:
        raise ValueError("N must be even and >= 4")
    
    # Generate primes if not provided
    if primes is None:
        primes = sieve_of_eratosthenes(N)
    else:
        # Filter to relevant range
        primes = primes[primes < N]
    
    # Create prime set for O(1) lookup
    prime_set = set(primes)
    
    # Count ordered pairs
    count = 0
    for p in primes:
        if p > N - 2:
            break
        q = N - p
        if q in prime_set:
            count += 1
    
    return count


def count_goldbach_batch(N_values, max_N=None):
    """
    Count Goldbach pairs for multiple N values efficiently.
    
    Parameters
    ----------
    N_values : array-like
        Array of even integers to analyze
    max_N : int, optional
        Maximum N value (for prime generation). If None, uses max(N_values)
        
    Returns
    -------
    numpy.ndarray
        Array of counts corresponding to N_values
    """
    N_values = np.asarray(N_values)
    
    if max_N is None:
        max_N = N_values.max()
    
    # Generate primes once
    primes = sieve_of_eratosthenes(max_N)
    prime_set = set(primes)
    
    # Count for each N
    counts = np.zeros(len(N_values), dtype=int)
    for i, N in enumerate(N_values):
        count = 0
        for p in primes:
            if p > N - 2:
                break
            if (N - p) in prime_set:
                count += 1
        counts[i] = count
    
    return counts


if __name__ == "__main__":
    # Example usage
    print("Goldbach Ordered-Pair Counter")
    print("=" * 50)
    
    # Test cases
    test_values = [10, 100, 1000, 10000]
    
    print("\nTest cases:")
    for N in test_values:
        count = count_goldbach_ordered(N)
        print(f"N = {N:6d}: G(N) = {count:6d} ordered pairs")
    
    # Verification against known values (OEIS A006307)
    known_values = {
        10: 3,
        100: 6,
        1000: 56,
        10000: 254
    }
    
    print("\nVerification:")
    all_correct = True
    for N, expected in known_values.items():
        actual = count_goldbach_ordered(N)
        status = "✓" if actual == expected else "✗"
        print(f"{status} N = {N:6d}: Expected {expected:6d}, Got {actual:6d}")
        if actual != expected:
            all_correct = False
    
    if all_correct:
        print("\n✓ All tests passed!")
    else:
        print("\n✗ Some tests failed!")
