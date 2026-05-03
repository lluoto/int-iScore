"""
Analysis utilities for generating int_iScore from CSV data.
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler


def load_csv_data(csv_path):
    """
    Load CSV data file.
    
    Args:
        csv_path: Path to CSV file
        
    Returns:
        DataFrame with loaded data
    """
    df = pd.read_csv(csv_path, encoding='ISO-8859-1')
    return df


DEFAULT_WEIGHTS = {
    'cpscore': 0.0281,
    'scaled_soap': 0.1502,
    'sc': 0.3052,
    'ec': 0.0931,
    'dope': 0.0954,
    'scaled_frustrate': 0.2567,
    'scaled_bsa': 0.0712,
}

DEFAULT_FEATURE_COLUMNS = [
    'cpscore', 'scaled_soap', 'sc', 'ec', 'dope', 'scaled_frustrate', 'scaled_bsa'
]


def calculate_correlation_weights(df, feature_columns, target_column='dockq'):
    """
    Calculate Pearson correlation weights between features and target.
    
    Args:
        df: DataFrame containing data
        feature_columns: List of feature column names
        target_column: Target column name (default: 'dockq')
        
    Returns:
        Tuple of (weights, correlations) where both are lists aligned with feature_columns
    """
    weights = []
    correlations = []
    
    for feature in feature_columns:
        if feature in df.columns and target_column in df.columns:
            corr_coefficient, p_value = pearsonr(df[feature], df[target_column])
            weights.append(abs(corr_coefficient))
            correlations.append(corr_coefficient)
        else:
            weights.append(0.0)
            correlations.append(0.0)
    
    return weights, correlations


def calculate_int_score(df, feature_columns, weights=None, scaler=None):
    """
    Calculate int_iScore using weighted combination of features.
    
    Args:
        df: DataFrame containing data
        feature_columns: List of feature column names
        weights: Optional list of weights (if None, uniform weights used)
        scaler: Optional StandardScaler instance
        
    Returns:
        Series with int_iScore values
    """
    if weights is None:
        weights = [1.0 / len(feature_columns)] * len(feature_columns)
    
    if scaler is None:
        scaler = StandardScaler()
    
    scaled_features = {}
    for feature in feature_columns:
        if feature in df.columns:
            scaled = scaler.fit_transform(df[feature].values.reshape(-1, 1)).flatten()
            scaled_features[feature] = scaled
        else:
            scaled_features[feature] = np.zeros(len(df))
    
    int_score = sum(
        scaled_features[feature] * weight 
        for feature, weight in zip(feature_columns, weights)
    )
    
    return int_score, scaler


def generate_weights_from_csv(csv_path, output_path=None):
    """
    Generate weights from CSV data using correlation with DockQ.
    
    Args:
        csv_path: Path to input CSV file
        output_path: Optional path to save weights CSV
        
    Returns:
        Tuple of (weight_dict, corr_dict)
    """
    df = load_csv_data(csv_path)
    
    df['frustrate'] = -df['frustrate']
    df['dope'] = -df['dope']
    df['soap'] = -df['soap']
    
    scaler = StandardScaler()
    for f in ['bsa', 'soap', 'frustrate', 'average_plddt', 'interface_plddt']:
        df[f'scaled_{f}'] = scaler.fit_transform(df[f].values.reshape(-1, 1)).flatten()
    
    feature_columns = ['cpscore','scaled_soap','sc','ec','dope','scaled_frustrate','scaled_bsa']
    available = [f for f in feature_columns if f in df.columns]
    
    weights, correlations = calculate_correlation_weights(df, available, 'dockq')
    
    weight_dict = dict(zip(available, weights))
    corr_dict = dict(zip(available, correlations))
    
    weight_sum = sum(weights)
    normalized = {k: v/weight_sum for k, v in weight_dict.items()}
    
    if output_path:
        result_df = pd.DataFrame({
            'feature': list(normalized.keys()),
            'weight': list(normalized.values())
        })
        result_df.to_csv(output_path, index=False)
    
    return weight_dict, corr_dict


def analyze_csv_correlations(csv_path):
    """
    Analyze all features in CSV and return correlation with DockQ.
    
    Args:
        csv_path: Path to input CSV file
        
    Returns:
        Dictionary mapping feature names to correlation coefficients
    """
    df = load_csv_data(csv_path)
    
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    numeric_columns = [c for c in numeric_columns if c != 'dockq']
    
    correlations = {}
    for col in numeric_columns:
        try:
            corr, pval = pearsonr(df[col], df['dockq'])
            correlations[col] = {'correlation': corr, 'p_value': pval, 'abs_correlation': abs(corr)}
        except Exception:
            correlations[col] = {'correlation': 0, 'p_value': 1, 'abs_correlation': 0}
    
    return correlations


def create_int_score_pipeline(csv_path, output_csv=None, include_iptm=True):
    """
    Complete pipeline to generate int_iScore from CSV.
    
    Args:
        csv_path: Path to input CSV file
        output_csv: Optional path to save result CSV
        include_iptm: Whether to include iPTM in the score
        
    Returns:
        DataFrame with int_iScore added
    """
    df = load_csv_data(csv_path)
    
    df['frustrate'] = -df['frustrate']
    df['dope'] = -df['dope']
    df['soap'] = -df['soap']
    
    scaler = StandardScaler()
    for f in ['bsa', 'soap', 'frustrate', 'average_plddt', 'interface_plddt']:
        if f in df.columns:
            df[f'scaled_{f}'] = scaler.fit_transform(df[f].values.reshape(-1, 1)).flatten()
    
    feature_columns = ['cpscore','scaled_soap','sc','ec','dope','scaled_frustrate','scaled_bsa']
    if include_iptm and 'iptm' in df.columns:
        feature_columns.append('iptm')
    
    available = [f for f in feature_columns if f in df.columns]
    
    weights, _ = calculate_correlation_weights(df, available, 'dockq')
    weight_sum = sum(weights)
    normalized_weights = [w/weight_sum for w in weights]
    
    df['int_iscore'] = sum(
        df[col] * weight 
        for col, weight in zip(available, normalized_weights)
    )
    
    if output_csv:
        df.to_csv(output_csv, index=False)
    
    return df


if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python analysis.py <csv_path> [output_csv]")
        print("Or import as module for CSV analysis")
        sys.exit(1)
    
    csv_path = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None
    
    df = create_int_score_pipeline(csv_path, output_csv)
    corr, pval = pearsonr(df['int_iscore'], df['dockq'])
    print(f"int_iscore vs DockQ correlation: {corr:.4f} (p={pval:.2e})")
