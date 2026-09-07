#!/usr/bin/env python3
"""
PON-Del Prediction Script

This script predicts pathogenicity of deletion variants using trained models.

IMPORTANT: Deletion start positions (where start_pos == 1) are always labeled as 
"P" (Pathogenic) and do not require model predictions, as these positions represent 
the start of deletions which are inherently pathogenic.

For all other positions, the model provides:
- Prediction probabilities (0-1 scale)
- Cross-validation predictions with uncertainty measures
- Classification labels (P=Pathogenic, B=Benign, U=Uncertain)
"""

import sys
import os
import pandas as pd
import numpy as np

# Set base directory
BASE_DIR = "/srv/shiny-server/pon_del"
TEMP_DIR = os.path.join(BASE_DIR, "temp")

# Add the code directory to Python path
sys.path.append('./code')

from train_00_function import load_model, process_X


def load_model_multiSplit(model_name, fold, version):
    """Load model with version prefix from multiSplit directory"""
    import lightgbm as lgb
    
    # Use Docker container path (mapped from host)
    # Host: /srv/data/website/structure_static/lu_psb_rshiny_server/pon_del
    # Container: /srv/shiny-server/pon_del
    base_dir = f'{BASE_DIR}/model/multiSplit'
    
    if model_name == "PON_Del":
        if fold == 'all':
            model_path = f'{base_dir}/{version}_{model_name}_train.txt'
        else:
            model_path = f'{base_dir}/{version}_{model_name}_fold{fold}.txt'
        if os.path.exists(model_path):
            return lgb.Booster(model_file=model_path)
        return None
    return None


def load_and_validate_data(task_id):
    """Load and validate input data for prediction."""
    input_file = os.path.join(TEMP_DIR, f'{task_id}_withfeat.csv')
    print(f"Reading input file: {input_file}")
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    df = pd.read_csv(input_file)
    print(f"Input file has {len(df)} rows")
    
    # Check if input column exists
    if 'input' in df.columns:
        print(f"Input column found with {len(df['input'])} values")
    else:
        print("Input column not found - this is normal for protein inputs")
    
    return df


def prepare_features(df):
    """Extract and prepare features for prediction."""
    # Extract key columns
    np_ids = df['NP_id']
    start_pos = df['start']
    end_pos = df['end']
    
    # Check if input column exists and extract it
    input_col = None
    if 'input' in df.columns:
        input_col = df['input'].copy()
        # print(f"Found input column with {len(input_col)} values")
        # print(f"Sample input values: {input_col.head(3).tolist()}")
    
    # Prepare features (drop non-feature columns)
    X_new = df.drop(['NP_id', 'end', 'input'], axis=1)
    
    # Load training data for imputation
    print("Loading training data...")
    X_train = pd.read_csv('/srv/shiny-server/db/pon_del/data/X_train.csv')
    X_train = X_train[X_new.columns]
    
    # Process features
    print("Processing features...")
    X_train_transformed, X_new_transformed = process_X(X_train, X_new, "PON_Del")
    
    return np_ids, start_pos, end_pos, input_col, X_new_transformed


def get_cv_predictions(X_new_transformed, model, start_pos, n_folds=5):
    """Get cross-validation predictions from 25 models (5 original + 20 multi-split) using bootstrap.
    
    Note: Positions with start_pos == 1 are always labeled as "P" (Pathogenic)
    and do not require model predictions as they represent deletion start positions.
    """
    all_predictions = []  # List of arrays, each shape (n_variants,)
    
    # 1. Get predictions from 5 folds of original split
    for fold in range(n_folds):
        model_cv = load_model("PON_Del", fold)
        if model_cv is not None:
            model_feats = model_cv.feature_name()
            X_new_model = X_new_transformed[model_feats]
            fold_pred = model_cv.predict(X_new_model)
            all_predictions.append(fold_pred)
        else:
            print(f"    Warning: Model not found for original split, fold {fold}")
    
    # 2. Get predictions from 4 versions (v2-v5) * 5 folds each
    versions = ['v2', 'v3', 'v4', 'v5']
    for version in versions:
        for fold in range(n_folds):
            model_cv = load_model_multiSplit("PON_Del", fold, version)
            if model_cv is not None:
                model_feats = model_cv.feature_name()
                X_new_model = X_new_transformed[model_feats]
                fold_pred = model_cv.predict(X_new_model)
                all_predictions.append(fold_pred)
            else:
                print(f"    Warning: Model not found for {version}, fold {fold}")
    
    if len(all_predictions) > 0:
        # Convert to numpy array: shape (n_models, n_variants)
        all_predictions_array = np.array(all_predictions)
        n_variants = all_predictions_array.shape[1]
        n_models = all_predictions_array.shape[0]
        
        # Calculate mean and std
        pred_cv_mean = np.mean(all_predictions_array, axis=0)
        pred_cv_std = np.std(all_predictions_array, axis=0)
        
        # Bootstrap to calculate p-values
        print("  Running bootstrap to get p-value...")
        n_bootstrap = 1000
        pvalues = np.zeros(n_variants)
        
        for i in range(n_variants):
            # Get all predictions for this variant
            variant_predictions = all_predictions_array[:, i]
            
            # Calculate mean prediction
            mean_pred = pred_cv_mean[i]
            
            # Bootstrap: resample with replacement 1000 times
            bootstrap_means = []
            for _ in range(n_bootstrap):
                # Resample from the predictions
                bootstrap_sample = np.random.choice(variant_predictions, size=n_models, replace=True)
                bootstrap_means.append(np.mean(bootstrap_sample))
            
            bootstrap_means = np.array(bootstrap_means)
            
            # Calculate p-value: test if mean is significantly different from 0.5
            # Two-tailed test using bootstrap
            if mean_pred > 0.5:
                # Test if significantly greater than 0.5
                pvalue = np.mean(bootstrap_means <= 0.5) * 2  # Two-tailed
            elif mean_pred < 0.5:
                # Test if significantly less than 0.5
                pvalue = np.mean(bootstrap_means >= 0.5) * 2  # Two-tailed
            else:
                # observed_mean == 0.5, p-value = 1.0
                pvalue = 1.0
            
            # Ensure p-value is between 0 and 1
            pvalue = min(1.0, max(0.0, pvalue))
            pvalues[i] = pvalue
        
        pred_cv_p_value = pvalues
        
        # Assign labels: P (Pathogenic), B (Benign), U (Uncertain)
        # Deletion start positions (start_pos == 1) are always "P"
        # For others: U (Uncertain) if p-value > 0.05 (not statistically significant)
        # Otherwise: P if mean > 0.5, B if mean <= 0.5
        p_threshold = 0.05
        pred_cv_label = np.where(start_pos == 1, "P", 
                                np.where(pred_cv_p_value > p_threshold, "U", 
                                        np.where(pred_cv_mean > 0.5, "P", "B")))
    else:
        print("Warning: No cross-validation models found, using single model predictions")
        pred_cv_mean = None
        pred_cv_std = None
        pred_cv_p_value = None
        pred_cv_label = None
    
    return pred_cv_mean, pred_cv_std, pred_cv_p_value, pred_cv_label


def create_prediction_labels(pred, pred_cv_mean, start_pos):
    """Create prediction labels based on probabilities.
    
    Note: Deletion start positions (start_pos == 1) are always labeled as "P" (Pathogenic)
    regardless of model predictions, as these positions represent the start of deletions.
    """
    # Single model labels
    # Deletion start positions are always "P", others based on prediction threshold
    single_label = np.where(start_pos == 1, "P", np.where(pred > 0.5, "P", "B"))
    
    # Cross-validation labels (fallback to single model if CV not available)
    if pred_cv_mean is not None:
        cv_label = np.where(start_pos == 1, "P", np.where(pred_cv_mean > 0.5, "P", "B"))
    else:
        cv_label = single_label
    
    return single_label, cv_label


def create_result_dataframe(np_ids, start_pos, end_pos, pred, pred_cv_mean, 
                          pred_cv_std, pred_cv_p_value, single_label, cv_label, input_col=None):
    """Create the final result dataframe with all predictions.
    
    Note: For deletion start positions (start_pos == 1):
    - All prediction values are set to 'P' (no numerical predictions)
    - P-values are set to 'P' (no statistical testing)
    - Labels are always 'P' (Pathogenic)
    """
    # Handle p-value: create array if None
    if pred_cv_p_value is None:
        pvalue_array = np.full(len(start_pos), 1.0)
    else:
        pvalue_array = pred_cv_p_value
    
    result_df = pd.DataFrame({
        'NP_id': np_ids,
        'start': start_pos,
        'end': end_pos,
        # For start_pos == 1: set to 'P' (no prediction value needed)
        # For others: use actual prediction values
        'pondel_pred': np.where(start_pos == 1, 'P', np.round(pred, 2)),
        'pondel_label': single_label,
        # P-value from bootstrap (25 models): 'N/A' for start positions, actual values for others
        'pondel_pred_pvalue': np.where(start_pos == 1, 'N/A', np.round(pvalue_array, 4)),
        'pondel_pred_label': cv_label
    })
    
    # Add input column as first column if it exists
    if input_col is not None:
        result_df.insert(0, 'Input', input_col)
    
    return result_df


def save_predictions(result_df, task_id):
    """Save predictions to CSV file."""
    output_file = os.path.join(TEMP_DIR, f'{task_id}_predictions.csv')
    print(f"Saving predictions to: {output_file}")
    result_df.to_csv(output_file, index=False)
    print("Predictions saved successfully")


def print_summary(results):
    """Print summary statistics of predictions."""
    n_total = len(results)
    n_P = np.sum(results['pondel_label'] == "P")
    n_B = np.sum(results['pondel_label'] == "B")
    n_U = np.sum(results['pondel_pred_label'] == "U")
    
    print("\nSummary:")
    print(f"Total positions: {n_total}")
    print(f"Pathogenic (P): {n_P}")
    print(f"Benign (B): {n_B}")
    print(f"Uncertain (U): {n_U}")


def predict_file(task_id):
    """Main prediction function."""
    print(f"Processing task ID: {task_id}")
    
    # Load the main model
    print("Loading model...")
    model = load_model("PON_Del", "all")
    if model is None:
        raise FileNotFoundError("Main model file not found. Please ensure the model has been trained.")
    
    # Load and validate data
    df = load_and_validate_data(task_id)
    
    # Prepare features
    np_ids, start_pos, end_pos, input_col, X_new_transformed = prepare_features(df)
    
    # Make single model predictions
    print("Making predictions...")
    pred = model.predict(X_new_transformed[model.feature_name()])
    
    # Get cross-validation predictions
    pred_cv_mean, pred_cv_std, pred_cv_p_value, pred_cv_label = get_cv_predictions(
        X_new_transformed, model, start_pos
    )
    
    # Create prediction labels
    single_label, cv_label = create_prediction_labels(pred, pred_cv_mean, start_pos)
    
    # Create result dataframe
    result_df = create_result_dataframe(
        np_ids, start_pos, end_pos, pred, pred_cv_mean, 
        pred_cv_std, pred_cv_p_value, single_label, cv_label, input_col
    )
    
    # Save predictions
    save_predictions(result_df, task_id)
    
    return result_df


def main():    
    if len(sys.argv) != 2:
        print("Usage: python 04_predict.py <task_id>")
        sys.exit(1)
        
    task_id = sys.argv[1]
    try:
        results = predict_file(task_id)
        print_summary(results)
        
    except FileNotFoundError as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        print(f"Error details: {type(e).__name__}")
        sys.exit(1)


if __name__ == "__main__":
    main() 
