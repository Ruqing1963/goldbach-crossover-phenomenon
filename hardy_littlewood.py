"""
Hardy-Littlewood Goldbach Formula Implementations
==================================================

This module implements both series expansion and logarithmic integral
formulations of the Hardy-Littlewood asymptotic formula.

Author: Ruqing Chen
Date: January 2026
"""

import numpy as np
from scipy import integrate
from scipy.special import gamma


# Twin prime constant
C2 = 0.6601618158468695739278121100145557784326233602847330039854


def singular_series(N):
    """
    Compute the singular series S(N) for the Hardy-Littlewood formula.
    
    S(N) = ∏ (1 - 1/(p-1)²) / (1 - 2/(p-1)²)
           p|N, p odd prime
    
    Parameters
    ----------
    N : int
        Even integer
        
    Returns
    -------
    float
        Singular series value (typically between 1.0 and 3.0)
    """
    if N % 2 != 0:
        raise ValueError("N must be even")
    
    # Factor N/2
    n = N // 2
    
    # Find odd prime factors
    S = 1.0
    p = 3
    temp_n = n
    
    while p * p <= temp_n:
        if temp_n % p == 0:
            # p divides N/2
            factor = (1 - 1/(p-1)**2) / (1 - 2/(p-1)**2)
            S *= factor
            # Remove all factors of p
            while temp_n % p == 0:
                temp_n //= p
        p += 2
    
    # Check if remaining number > 1 (it's a prime factor)
    if temp_n > 1:
        p = temp_n
        factor = (1 - 1/(p-1)**2) / (1 - 2/(p-1)**2)
        S *= factor
    
    return S


def hardy_littlewood_series(N, order=4):
    """
    Hardy-Littlewood formula using series expansion.
    
    r(N) ≈ 2C₂S(N) × N/log²(N) × [1 + a₁/log(N) + a₂/log²(N) + ...]
    
    Parameters
    ----------
    N : int
        Even integer
    order : int, optional
        Expansion order (default: 4)
        
    Returns
    -------
    float
        Predicted number of Goldbach representations
    """
    if N < 4 or N % 2 != 0:
        raise ValueError("N must be even and >= 4")
    
    log_N = np.log(N)
    S_N = singular_series(N)
    
    # Leading term
    result = 2 * C2 * S_N * N / (log_N ** 2)
    
    # Series correction terms (empirically fitted)
    # These coefficients come from asymptotic expansion
    if order >= 1:
        a1 = 2.0
        result *= (1 + a1 / log_N)
    
    if order >= 2:
        a2 = 2.0
        result *= (1 + a2 / (log_N ** 2))
    
    if order >= 3:
        a3 = 4.0 / 3.0
        result *= (1 + a3 / (log_N ** 3))
    
    if order >= 4:
        a4 = 2.0 / 3.0
        result *= (1 + a4 / (log_N ** 4))
    
    return result


def hardy_littlewood_integral(N, limit=200):
    """
    Hardy-Littlewood formula using logarithmic integral.
    
    r(N) ≈ 2C₂S(N) × ∫₂^(N-2) dt / [log(t)log(N-t)]
    
    Parameters
    ----------
    N : int
        Even integer
    limit : int, optional
        Integration limit for adaptive quadrature (default: 200)
        
    Returns
    -------
    float
        Predicted number of Goldbach representations
    """
    if N < 4 or N % 2 != 0:
        raise ValueError("N must be even and >= 4")
    
    S_N = singular_series(N)
    
    # Define integrand
    def integrand(t):
        if t <= 2 or t >= N - 2:
            return 0.0
        return 1.0 / (np.log(t) * np.log(N - t))
    
    # Numerical integration using adaptive Gaussian quadrature
    integral_value, error = integrate.quad(
        integrand, 
        2, 
        N - 2,
        limit=limit
    )
    
    result = 2 * C2 * S_N * integral_value
    
    return result


def compare_methods(N):
    """
    Compare series and integral methods for a given N.
    
    Parameters
    ----------
    N : int
        Even integer
        
    Returns
    -------
    dict
        Dictionary with 'series', 'integral', 'singular_series' values
    """
    return {
        'N': N,
        'series': hardy_littlewood_series(N),
        'integral': hardy_littlewood_integral(N),
        'singular_series': singular_series(N)
    }


if __name__ == "__main__":
    # Example usage
    print("Hardy-Littlewood Formula Implementations")
    print("=" * 60)
    
    # Test values
    test_N = [100, 1000, 10000, 100000, 1000000]
    
    print("\nComparison of Methods:\n")
    print(f"{'N':>10} {'Series':>12} {'Integral':>12} {'S(N)':>8}")
    print("-" * 60)
    
    for N in test_N:
        result = compare_methods(N)
        print(f"{N:10d} {result['series']:12.2f} {result['integral']:12.2f} "
              f"{result['singular_series']:8.4f}")
    
    # Detailed example for N = 10^6
    print("\n" + "=" * 60)
    print("Detailed Analysis for N = 1,000,000:")
    print("=" * 60)
    
    N = 1000000
    S_N = singular_series(N)
    series = hardy_littlewood_series(N)
    integral = hardy_littlewood_integral(N)
    
    print(f"\nSingular series S(N) = {S_N:.6f}")
    print(f"\nSeries expansion:")
    print(f"  Prediction: {series:.2f}")
    print(f"\nLogarithmic integral:")
    print(f"  Prediction: {integral:.2f}")
    print(f"\nDifference: {abs(series - integral):.2f}")
    print(f"Relative difference: {abs(series - integral) / integral * 100:.2f}%")
