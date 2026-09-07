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
from scipy.stats import t as pt

# Set base directory
BASE_DIR = "/srv/shiny-server/pon_del"
TEMP_DIR = os.path.join(BASE_DIR, "temp")

# Add the code directory to Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

from train_00_function import load_model, process_X


def load_and_validate_data(task_id):
    """Load and validate input data for prediction."""
    input_file = os.path.join(TEMP_DIR, f'{task_id}_withfeat.csv')
    print(f"Reading input file: {input_file}")
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    df = pd.read_csv(input_file)
    print(f"Input file has {len(df)} rows")
    return df


def prepare_features(df):
    """Extract and prepare features for prediction."""
    # Extract key columns
    np_ids = df['NP_id']
    start_pos = df['start']
    end_pos = df['end']
    
    # Prepare features (drop non-feature columns)
    X_new = df.drop(['NP_id', 'end'], axis=1)
    
    # Load training data for imputation
    print("Loading training data...")
    X_train = pd.read_csv('/srv/shiny-server/db/pon_del/data/X_train.csv')
    X_train = X_train[X_new.columns]
    
    # Process features
    print("Processing features...")
    X_train_transformed, X_new_transformed = process_X(X_train, X_new, "PON_Del")
    
    return np_ids, start_pos, end_pos, X_new_transformed


def get_cv_predictions(X_new_transformed, model, start_pos, n_folds=5):
    """Get cross-validation predictions from multiple models.
    
    Note: Positions with start_pos == 1 are always labeled as "P" (Pathogenic)
    and do not require model predictions as they represent deletion start positions.
    """
    pred_cv = []
    
    for i in range(n_folds):
        model_cv = load_model("PON_Del", i)
        if model_cv is not None:
            pred_cv.append(model_cv.predict(X_new_transformed[model.feature_name()]))
        else:
            print(f"Warning: Model for fold {i} not found")
    
    if len(pred_cv) > 0:
        pred_cv = np.array(pred_cv)
        pred_cv_mean = np.mean(pred_cv, axis=0)
        pred_cv_std = np.std(pred_cv, axis=0)
        
        # Calculate t-statistic and p-value for uncertainty assessment
        # Only for positions where start_pos != 1 (not deletion start positions)
        pred_cv_t_stat = (pred_cv_mean - 0.5) / (pred_cv_std / np.sqrt(len(pred_cv)))
        pred_cv_p_value = 2 * pt.cdf(-abs(pred_cv_t_stat), df=len(pred_cv)-1)
        
        # Assign labels: P (Pathogenic), B (Benign), U (Uncertain)
        # Deletion start positions (start_pos == 1) are always "P"
        # For others: U (Uncertain) if p-value > 0.05 (not statistically significant)
        # Otherwise: P if mean > 0.5, B if mean <= 0.5
        pred_cv_label = np.where(start_pos == 1, "P", 
                                np.where(pred_cv_p_value > 0.05, "U", 
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
                          pred_cv_std, pred_cv_p_value, single_label, cv_label):
    """Create the final result dataframe with all predictions.
    
    Note: For deletion start positions (start_pos == 1):
    - All prediction values are set to 'P' (no numerical predictions)
    - Standard deviations are set to 'P' (no uncertainty measures)
    - P-values are set to 'P' (no statistical testing)
    - Labels are always 'P' (Pathogenic)
    """
    return pd.DataFrame({
        'NP_id': np_ids,
        'start': start_pos,
        'end': end_pos,
        # For start_pos == 1: set to 'P' (no prediction value needed)
        # For others: use actual prediction values
        'pondel_pred': np.where(start_pos == 1, 'P', np.round(pred, 2)),
        'pondel_label': single_label,
        # Cross-validation predictions: 'N/A' for start positions, actual values for others
        'pondel_pred_cv': np.where(start_pos == 1, 'N/A', 
                                  np.round(pred_cv_mean if pred_cv_mean is not None else pred, 2)),
        # Standard deviations: 'N/A' for start positions, actual values for others
        'pondel_pred_cv_std': np.where(start_pos == 1, 'N/A', 
                                      np.round(pred_cv_std if pred_cv_std is not None else 0, 2)),
        # P-values: 'N/A' for start positions, actual values for others
        'pondel_pred_cv_pvalue': np.where(start_pos == 1, 'N/A', 
                                         np.round(pred_cv_p_value if pred_cv_p_value is not None else 1.0, 3)),
        'pondel_pred_cv_label': cv_label
    })


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
    n_U = np.sum(results['pondel_pred_cv_label'] == "U")
    
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
    np_ids, start_pos, end_pos, X_new_transformed = prepare_features(df)
    
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
        pred_cv_std, pred_cv_p_value, single_label, cv_label
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
