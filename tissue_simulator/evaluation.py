"""
Evaluation metrics for comparing graph statistics.

This module provides various metrics for comparing the statistical
properties of colored graphs, including divergence measures and
distance metrics.
"""

import numpy as np
from scipy.stats import entropy
from typing import Dict, List, Tuple


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """
    Calculate the cosine similarity between two vectors.

    Args:
        a: First vector (1D numpy array)
        b: Second vector (1D numpy array)

    Returns:
        Cosine similarity between a and b (0 to 1, where 1 is identical)
    """
    # Convert inputs to numpy arrays
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if a.sum() == 0 and b.sum() == 0:
        return 1.0  # Both zero vectors are considered identical
    if a.sum() == 0 or b.sum() == 0:
        return 0.0  # One is zero, the other isn't
    
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    
    return dot_product / (norm_a * norm_b)


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    """
    Calculate the cosine distance between two vectors.

    Args:
        a: First vector (1D numpy array)
        b: Second vector (1D numpy array)

    Returns:
        Cosine distance between a and b (0 to 2, where 0 is identical)
    """
    return 1 - cosine_similarity(a, b)


def js_divergence(p: np.ndarray, q: np.ndarray, base: int = 2) -> float:
    """
    Calculate the Jensen-Shannon Divergence between two probability distributions.

    The JS Divergence is a symmetric and smoothed measure of the similarity between
    two probability distributions. It is based on the Kullback-Leibler (KL) Divergence.
    The result is in bits if the base is 2. A score of 0 means the distributions
    are identical, while a larger value indicates greater divergence.

    Args:
        p: First distribution (1D numpy array). Will be normalized to sum to 1.
        q: Second distribution (1D numpy array). Will be normalized to sum to 1.
        base: Logarithmic base to use for calculation (default: 2)

    Returns:
        Jensen-Shannon Divergence between p and q (0 to 1 for base 2)
    """
    # Convert inputs to numpy arrays
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    if p.sum() == 0 and q.sum() == 0:
        return 0.0
    if p.sum() == 0 or q.sum() == 0:
        # If one distribution is all zeros, the divergence is maximal
        # For base 2, the max is 1.0
        return 1.0

    # Normalize the vectors to get probability distributions
    p /= p.sum()
    q /= q.sum()

    # Calculate the mixture distribution (M)
    m = 0.5 * (p + q)

    # Calculate the JS Divergence using the KL Divergence (entropy)
    # JSD(P||Q) = 0.5 * KLD(P||M) + 0.5 * KLD(Q||M)
    jsd = 0.5 * entropy(p, m, base=base) + 0.5 * entropy(q, m, base=base)
    
    return jsd


def calculate_percent_difference(source_val: float, target_val: float) -> float:
    """
    Calculate percent difference between two values.
    
    Args:
        source_val: Source/target value
        target_val: Generated/actual value
        
    Returns:
        Percent difference
    """
    if source_val == 0 and target_val == 0:
        return 0.0
    elif source_val == 0:
        return float('inf')
    else:
        return (abs(target_val - source_val) / source_val) * 100


def mean_absolute_error(source: np.ndarray, target: np.ndarray) -> float:
    """
    Calculate mean absolute error between two arrays.
    
    Args:
        source: Source values
        target: Target values
        
    Returns:
        Mean absolute error
    """
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    
    return np.mean(np.abs(source - target))


def root_mean_squared_error(source: np.ndarray, target: np.ndarray) -> float:
    """
    Calculate root mean squared error between two arrays.
    
    Args:
        source: Source values
        target: Target values
        
    Returns:
        Root mean squared error
    """
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    
    return np.sqrt(np.mean((source - target) ** 2))


def evaluate_graph_coloring(source_stats: Dict, target_stats: Dict) -> Dict:
    """
    Comprehensive evaluation of graph coloring quality.
    
    Args:
        source_stats: Target/desired statistics
        target_stats: Generated statistics
        
    Returns:
        Dictionary with multiple evaluation metrics
    """
    # Extract edge statistics for comparison
    edge_keys = [k for k in source_stats if k.startswith('edges_')]
    source_edges = [source_stats[k] for k in edge_keys]
    target_edges = [target_stats.get(k, 0) for k in edge_keys]
    
    # Calculate various metrics
    results = {
        'mean_absolute_error': mean_absolute_error(source_edges, target_edges),
        'root_mean_squared_error': root_mean_squared_error(source_edges, target_edges),
        'cosine_similarity': cosine_similarity(source_edges, target_edges),
        'cosine_distance': cosine_distance(source_edges, target_edges),
        'js_divergence': js_divergence(source_edges, target_edges),
    }
    
    # Calculate percent differences for each edge type
    percent_diffs = []
    for key in edge_keys:
        source_val = source_stats[key]
        target_val = target_stats.get(key, 0)
        diff = calculate_percent_difference(source_val, target_val)
        if diff != float('inf'):
            percent_diffs.append(diff)
    
    results['avg_percent_difference'] = np.mean(percent_diffs) if percent_diffs else 0.0
    results['max_percent_difference'] = np.max(percent_diffs) if percent_diffs else 0.0
    
    return results


def print_evaluation_report(source_stats: Dict, target_stats: Dict, 
                           evaluation: Dict = None):
    """
    Print a comprehensive evaluation report.
    
    Args:
        source_stats: Target/desired statistics
        target_stats: Generated statistics
        evaluation: Pre-computed evaluation metrics (optional)
    """
    if evaluation is None:
        evaluation = evaluate_graph_coloring(source_stats, target_stats)
    
    print("\n" + "="*60)
    print("GRAPH COLORING EVALUATION REPORT")
    print("="*60)
    
    print("\n--- Overall Metrics ---")
    print(f"Mean Absolute Error:       {evaluation['mean_absolute_error']:.4f}")
    print(f"Root Mean Squared Error:   {evaluation['root_mean_squared_error']:.4f}")
    print(f"Cosine Similarity:         {evaluation['cosine_similarity']:.4f}")
    print(f"Cosine Distance:           {evaluation['cosine_distance']:.4f}")
    print(f"Jensen-Shannon Divergence: {evaluation['js_divergence']:.4f}")
    print(f"Avg Percent Difference:    {evaluation['avg_percent_difference']:.2f}%")
    print(f"Max Percent Difference:    {evaluation['max_percent_difference']:.2f}%")
    
    print("\n--- Node Counts ---")
    node_keys = [k for k in source_stats if k.startswith('nodes_')]
    for key in node_keys:
        source_val = source_stats.get(key, 0)
        target_val = target_stats.get(key, 0)
        diff = calculate_percent_difference(source_val, target_val)
        print(f"{key:20s}: Source={source_val:4d}, Target={target_val:4d}, Diff={diff:6.2f}%")
    
    print("\n--- Edge Counts ---")
    edge_keys = [k for k in source_stats if k.startswith('edges_')]
    for key in edge_keys:
        source_val = source_stats.get(key, 0)
        target_val = target_stats.get(key, 0)
        diff = calculate_percent_difference(source_val, target_val)
        print(f"{key:20s}: Source={source_val:4d}, Target={target_val:4d}, Diff={diff:6.2f}%")
    
    print("\n" + "="*60)
