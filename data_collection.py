"""
Data Collection Script for Goldbach Crossover Analysis
=======================================================

This script generates the complete dataset used in the paper.

WARNING: Full data collection (N up to 10^8) takes 4-6 hours!

Author: Ruqing Chen
Date: January 2026
"""

import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import time

# Import our modules
from goldbach_counter import sieve_of_eratosthenes, count_goldbach_ordered
from hardy_littlewood import hardy_littlewood_series, hardy_littlewood_integral


def generate_sampling_points(min_n=1000, max_n=100000000):
    """
    Generate stratified sampling points across logarithmic scales.
    
    Dense sampling near crossover region (N ≈ 10^5),
    sparse sampling at extremes.
    
    Parameters
    ----------
    min_n : int
        Minimum N value
    max_n : int
        Maximum N value
        
    Returns
    -------
    numpy.ndarray
        Array of even integers to sample
    """
    points = []
    
    # Very dense near crossover (10^4 to 10^6)
    points.extend(range(10000, 100000, 1000))      # Every 1k
    points.extend(range(100000, 1000000, 5000))    # Every 5k
    
    # Medium density in other ranges
    points.extend(range(1000, 10000, 500))         # Every 500
    points.extend(range(1000000, 10000000, 50000)) # Every 50k
    
    # Sparse at large N
    points.extend(range(10000000, 100000000, 500000)) # Every 500k
    
    # Add key milestones
    milestones = [1000, 10000, 100000, 1000000, 10000000, 
                  54950000, 70000000, 97000000, 100000000]
    points.extend(milestones)
    
    # Convert to array, ensure even, remove duplicates
    points = np.array(points)
    points = points[points % 2 == 0]  # Keep only even
    points = np.unique(points)
    
    # Filter to range
    points = points[(points >= min_n) & (points <= max_n)]
    
    return np.sort(points)


def collect_data(N_values, primes):
    """
    Collect Goldbach data for array of N values.
    
    Parameters
    ----------
    N_values : numpy.ndarray
        Array of even integers
    primes : numpy.ndarray
        Pre-computed prime array
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: N, G_N, Series, Integral, Series_Bias, Integral_Bias
    """
    results = []
    
    total = len(N_values)
    start_time = time.time()
    
    for i, N in enumerate(N_values):
        # Actual count
        G_N = count_goldbach_ordered(N, primes)
        
        # Predictions
        series_pred = hardy_littlewood_series(N)
        integral_pred = hardy_littlewood_integral(N)
        
        # Biases
        series_bias = (series_pred - G_N) / G_N
        integral_bias = (integral_pred - G_N) / G_N
        
        results.append({
            'N': N,
            'G_N': G_N,
            'Series_Pred': series_pred,
            'Integral_Pred': integral_pred,
            'Series_Bias': series_bias,
            'Integral_Bias': integral_bias,
            'Abs_Bias_Series': abs(series_bias),
            'Abs_Bias_Integral': abs(integral_bias)
        })
        
        # Progress update
        if (i + 1) % 100 == 0 or i == total - 1:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            remaining = (total - i - 1) / rate if rate > 0 else 0
            print(f"Progress: {i+1}/{total} ({100*(i+1)/total:.1f}%) "
                  f"- {elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining")
    
    return pd.DataFrame(results)


def main():
    """Main data collection function."""
    parser = argparse.ArgumentParser(
        description='Generate Goldbach crossover dataset'
    )
    parser.add_argument('--min-n', type=int, default=1000,
                        help='Minimum N value (default: 1000)')
    parser.add_argument('--max-n', type=int, default=100000,
                        help='Maximum N value (default: 100000)')
    parser.add_argument('--output', type=str, default='goldbach_data.csv',
                        help='Output CSV filename')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Goldbach Crossover Data Collection")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Min N: {args.min_n:,}")
    print(f"  Max N: {args.max_n:,}")
    print(f"  Output: {args.output}")
    
    # Generate sampling points
    print(f"\nGenerating sampling points...")
    N_values = generate_sampling_points(args.min_n, args.max_n)
    print(f"  Total points: {len(N_values):,}")
    
    # Generate primes
    print(f"\nGenerating primes up to {args.max_n:,}...")
    start_time = time.time()
    primes = sieve_of_eratosthenes(args.max_n)
    elapsed = time.time() - start_time
    print(f"  Generated {len(primes):,} primes in {elapsed:.1f}s")
    
    # Collect data
    print(f"\nCollecting data...")
    print(f"  (This may take a while for large max_n)")
    df = collect_data(N_values, primes)
    
    # Save results
    print(f"\nSaving results to {args.output}...")
    df.to_csv(args.output, index=False)
    
    # Summary statistics
    print(f"\n" + "=" * 70)
    print("Summary Statistics:")
    print("=" * 70)
    print(f"\nTotal points collected: {len(df):,}")
    print(f"\nSeries expansion:")
    print(f"  Mean |bias|: {df['Abs_Bias_Series'].mean()*100:.2f}%")
    print(f"  Median |bias|: {df['Abs_Bias_Series'].median()*100:.2f}%")
    print(f"\nLogarithmic integral:")
    print(f"  Mean |bias|: {df['Abs_Bias_Integral'].mean()*100:.2f}%")
    print(f"  Median |bias|: {df['Abs_Bias_Integral'].median()*100:.2f}%")
    
    # Find crossover point
    series_wins = df['Abs_Bias_Series'] < df['Abs_Bias_Integral']
    integral_wins = df['Abs_Bias_Integral'] < df['Abs_Bias_Series']
    
    if series_wins.any() and integral_wins.any():
        crossover_idx = np.where(integral_wins)[0][0]
        crossover_N = df.iloc[crossover_idx]['N']
        print(f"\nCrossover point: N ≈ {crossover_N:,}")
    
    # Find peak advantage
    df['Advantage_Ratio'] = df['Abs_Bias_Series'] / df['Abs_Bias_Integral']
    peak_idx = df['Advantage_Ratio'].idxmax()
    peak_N = df.iloc[peak_idx]['N']
    peak_ratio = df.iloc[peak_idx]['Advantage_Ratio']
    print(f"Peak advantage: {peak_ratio:.1f}× at N = {peak_N:,}")
    
    print(f"\n✓ Data collection complete!")
    print(f"  Output saved to: {args.output}")


if __name__ == "__main__":
    main()
