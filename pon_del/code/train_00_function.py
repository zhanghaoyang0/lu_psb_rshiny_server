# =================================================
# prepare, now not used bc put in 01 for better readability
# =================================================
# conda activate py38
import os
# import shap
import numpy as np
import pickle
import pandas as pd
# import matplotlib.pyplot as plt
# import torch
# import torch.nn as nn
# import torch.nn.functional as F
# from torch.nn.utils.rnn import pad_sequence
# from torch.utils.data import DataLoader, TensorDataset
from collections import Counter
import lightgbm as lgb
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold
from sklearn.feature_selection import RFE
from lightgbm import LGBMClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score, precision_recall_curve
from sklearn.metrics import (
    accuracy_score, matthews_corrcoef, precision_score, recall_score, confusion_matrix,
    roc_auc_score
)
from collections import defaultdict

# =================================================
# Settings
# =================================================
seed = 2025

# LightGBM parameters
lgb_params = {
    'objective': 'binary', 
    'metric': ['binary_logloss'], 
    'boosting_type': 'gbdt',
    'learning_rate': 0.01,
    'num_leaves': 10,
    'feature_fraction': 0.8,
    'bagging_fraction': 0.8,
    'bagging_freq': 5,
    'seed': seed
}

# Update with exact parameters from optimization
lgbm_best_params = {'boosting_type': 'gbdt', 'learning_rate': 0.004825553202546152, 'num_leaves': 20, 'feature_fraction': 0.9, 'bagging_fraction': 0.9, 'bagging_freq': 1, 'min_child_samples': 10, 'min_child_weight': 0.006618664867480504, 'min_split_gain': 0.1, 'reg_alpha': 0.4, 'reg_lambda': 0.7}
# Initialize KFold
kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)


# =================================================
# Functions
# =================================================
# cal_freq(df, ['class', 'set'])
def cal_freq(df, group_by_columns):
    # Calculate frequency counts
    freq = df.groupby(group_by_columns).size()
    
    # Calculate total counts
    total_counts = freq.sum()
    
    # Calculate proportions
    proportions = freq / total_counts
    
    # Combine into a DataFrame
    result = pd.DataFrame({'count': freq, 'proportion': proportions}).reset_index()
    
    return result


# Fill missing values
def fill_na(dataset_train, dataset_val=None):
    dataset_train = dataset_train.copy()
    if dataset_val is not None:
        dataset_val = dataset_val.copy()
    
    # Separate binary and continuous features
    binary_cols = [col for col in dataset_train.columns if dataset_train[col].dtype in [int, float] and set(dataset_train[col].dropna().unique()) <= {0, 1}]
    continuous_cols = [col for col in dataset_train.columns if dataset_train[col].dtype in [int, float] and col not in binary_cols]
    
    # Create imputers
    binary_imputer = SimpleImputer(strategy="most_frequent")  # For binary features
    continuous_imputer = SimpleImputer(strategy="mean")  # For continuous features
    
    # Fit and transform training data
    dataset_train[binary_cols] = binary_imputer.fit_transform(dataset_train[binary_cols])
    dataset_train[continuous_cols] = continuous_imputer.fit_transform(dataset_train[continuous_cols])
    
    # Transform validation data if provided
    if dataset_val is not None:
        dataset_val[binary_cols] = binary_imputer.transform(dataset_val[binary_cols])
        dataset_val[continuous_cols] = continuous_imputer.transform(dataset_val[continuous_cols])
        return dataset_train, dataset_val
    
    return dataset_train

# Save model
def save_model(model, model_name, fold):
    os.makedirs("./model", exist_ok=True)  # Ensure directory exists
    if model_name == "LightGBM" or model_name == "PON_Del":
        if fold == 'all':
            model.save_model(f'./model/{model_name}_train.txt')
        else:
            model.save_model(f'./model/{model_name}_fold{fold}.txt')
    elif model_name in ["DL_MLP", "DL_CNN", "DL_GRU"]:
        if fold == 'all':
            torch.save(model.state_dict(), f'./model/{model_name}_train.pt')
        else:
            torch.save(model.state_dict(), f'./model/{model_name}_fold{fold}.pt')
    else:
        if fold == 'all':
            with open(f'./model/{model_name}_train.pkl', 'wb') as f:
                pickle.dump(model, f)
        else:
            with open(f'./model/{model_name}_fold{fold}.pkl', 'wb') as f:
                pickle.dump(model, f)

# Load model
def load_model(model_name, fold, input_dim=None):
    if model_name == "LightGBM" or model_name == "PON_Del":
        if fold == 'all':
            model_path = f'./model/{model_name}_train.txt'
        else:
            model_path = f'./model/{model_name}_fold{fold}.txt'
        if os.path.exists(model_path):
            return lgb.Booster(model_file=model_path)
        return None # Return None
    elif model_name in ["DL_MLP", "DL_CNN", "DL_GRU"]:
        if fold == 'all':
            model_path = f'./model/{model_name}_train.pt'
        else:
            model_path = f'./model/{model_name}_fold{fold}.pt'
        if os.path.exists(model_path):
            if input_dim is None:
                raise ValueError(f"{model_name} model requires input_dim to load model.")
            
            # Create the appropriate model based on model_name
            if model_name == "DL_MLP":
                model = DL_MLP(input_dim)
            elif model_name == "DL_CNN":
                model = DL_CNN(input_dim=input_dim)
            elif model_name == "DL_GRU":
                model = DL_GRU(input_dim=input_dim)
                
            model.load_state_dict(torch.load(model_path, weights_only=True))
            model.eval()  # Set to evaluation mode
            return model
        return None # Return None
    else:
        if fold == 'all':
            model_path = f'./model/{model_name}_train.pkl'
        else:
            model_path = f'./model/{model_name}_fold{fold}.pkl'
        if os.path.exists(model_path):
            with open(model_path, 'rb') as f:
                return pickle.load(f)
        return None  # Return None


# Calculate metrics
def calculate_metrics(y_true, y_prob, y_label=None):
    """Calculate various classification metrics from true labels and predicted probabilities.
    Args:
        y_true: True labels
        y_prob: Predicted probabilities
        y_label: Optional binary labels. If provided, use these instead of calculating from y_prob
    """
    y_true = np.asarray(y_true, dtype=np.int64)
    y_prob = np.asarray(y_prob, dtype=np.float64)
    
    # Use provided labels or calculate from probabilities
    if y_label is not None:
        y_pred = np.asarray(y_label, dtype=np.int64)
    else:
        y_pred = (y_prob > 0.5).astype(np.int64)
    
    # Get confusion matrix values and calculate denominators
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    denoms = {'acc': tp + tn + fp + fn, 'ppv': tp + fp, 'npv': tn + fn, 'sens': tp + fn, 'spec': tn + fp}
    
    # Calculate basic metrics
    metrics = {
        "Accuracy": (tp + tn) / denoms['acc'],
        "MCC": ((tp * tn) - (fp * fn)) / np.sqrt(denoms['sens'] * denoms['ppv'] * denoms['npv'] * denoms['spec']),
        "PPV": tp / denoms['ppv'] if denoms['ppv'] > 0 else 0,
        "NPV": tn / denoms['npv'] if denoms['npv'] > 0 else 0,
        "Sensitivity": tp / denoms['sens'] if denoms['sens'] > 0 else 0,
        "Specificity": tn / denoms['spec'] if denoms['spec'] > 0 else 0,
        "AUC": roc_auc_score(y_true, y_prob)
    }
    
    # Calculate OPM
    metrics["OPM"] = ((metrics["PPV"] + metrics["NPV"]) * (metrics["Sensitivity"] + metrics["Specificity"]) * (metrics["Accuracy"] + (1 + metrics["MCC"])/2)) / 8
    
    # Calculate normalized metrics
    factor =  (tp + fn) / (tn + fp)
    n_tp, n_fn = tp, fn
    n_tn, n_fp = round(tn * factor), round(fp * factor)
    n_denoms = {'acc': n_tp + n_tn + n_fp + n_fn, 'ppv': n_tp + n_fp, 'npv': n_tn + n_fn, 'sens': n_tp + n_fn, 'spec': n_tn + n_fp}
    
    # Add normalized metrics
    metrics.update({
        "n_Accuracy": (n_tp + n_tn) / n_denoms['acc'],
        "n_MCC": ((n_tp * n_tn) - (n_fp * n_fn)) / np.sqrt(n_denoms['sens'] * n_denoms['ppv'] * n_denoms['npv'] * n_denoms['spec']),
        "n_PPV": n_tp / n_denoms['ppv'] if n_denoms['ppv'] > 0 else 0,
        "n_NPV": n_tn / n_denoms['npv'] if n_denoms['npv'] > 0 else 0,
        "n_Sensitivity": n_tp / n_denoms['sens'] if n_denoms['sens'] > 0 else 0,
        "n_Specificity": n_tn / n_denoms['spec'] if n_denoms['spec'] > 0 else 0
    })
    
    # Calculate normalized OPM
    metrics["n_OPM"] = ((metrics["n_PPV"] + metrics["n_NPV"]) * (metrics["n_Sensitivity"] + metrics["n_Specificity"]) * (metrics["n_Accuracy"] + (1 + metrics["n_MCC"])/2)) / 8
    
    # Round all metrics and add confusion matrix values
    metrics = {k: round(v, 3) for k, v in metrics.items()}
    metrics.update({"tp": tp, "tn": tn, "fp": fp, "fn": fn, "n_tp": n_tp, "n_tn": n_tn, "n_fp": n_fp, "n_fn": n_fn})
    
    return metrics

# Process data
def process_X(X_train_sub, X_val_sub, model_name):
    # Get the list of sequence columns that actually exist in the DataFrame
    seq_cols = ["del_seq", "up5seq", "down5seq"]
    
    # Only drop columns if they exist in both DataFrames
    common_cols = list(set(X_train_sub.columns) & set(X_val_sub.columns))
    existing_seq_cols = [col for col in seq_cols if col in common_cols]
    
    # Only drop columns if they exist
    if existing_seq_cols:
        X_train_sub = X_train_sub.drop(columns=existing_seq_cols)
        X_val_sub = X_val_sub.drop(columns=existing_seq_cols)
    
    # Fill missing values
    X_train_sub, X_val_sub = fill_na(X_train_sub, X_val_sub)
    
    # Transform data for SVM, LogisticRegression, and DL_MLP
    if model_name in ["SVM", "LR", "DL_MLP"]:
        scaler = MinMaxScaler()
        return scaler.fit_transform(X_train_sub), scaler.transform(X_val_sub)
    return X_train_sub, X_val_sub  # Return unchanged data for other models


# Sequence encoding
def seq_encode(seqs):
    standard_aa = list("ACDEFGHIKLMNPQRSTVWY")
    char_to_index = {char: i + 1 for i, char in enumerate(standard_aa)}
    alternative_encoding = 0  # if not in 20 aa, set to 0   
    # Convert series to list if pandas Series is passed
    if isinstance(seqs, pd.Series):
        seqs = seqs.tolist()
    # Encode sequences
    encoded_seqs = [
        torch.tensor([ char_to_index.get(c, alternative_encoding) for c in seq])
        for seq in seqs
    ]
    # Pad the sequences
    padded_seqs = pad_sequence(encoded_seqs, batch_first=True, padding_value=0)
    return padded_seqs


def train_DL(model_name, X_train_transformed, y_train, input_dim=None, seq_data=None, fold=None, epochs=300, 
             is_final=False, patience=30, val_data=None, learning_rate=0.001):
    """Train deep learning models (DL_MLP, DL_CNN, or DL_GRU) with the given data
    Args:
        model_name: Name of the model ('DL_MLP', 'DL_CNN', or 'DL_GRU')
        X_train_transformed: transformed training features (numpy array or DataFrame)
        y_train: training labels
        input_dim: input dimension for the model
        seq_data: dictionary containing sequence data for DL_CNN/DL_GRU {'del_seq': tensor}
        fold: fold number or 'final' for final training
        epochs: number of training epochs
        is_final: whether this is final training
        patience: number of epochs to wait for improvement before early stopping
        val_data: validation data for early stopping (X_val, y_val, val_seq_data)
        learning_rate: learning rate for optimizer (default: 0.001)
    Returns:
        trained model, best epoch
    """
    # Set random seeds for reproducibility
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
    # Convert features to numpy if DataFrame
    if isinstance(X_train_transformed, pd.DataFrame):
        X_train_transformed = X_train_transformed.values
    
    # Convert regular features to PyTorch tensors
    X_train_tensor = torch.tensor(X_train_transformed, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train.values if isinstance(y_train, pd.Series) else y_train, dtype=torch.float32).reshape(-1, 1)
    
    # Prepare validation data if provided
    if val_data is not None:
        X_val, y_val, val_seq_data = val_data
        if isinstance(X_val, pd.DataFrame):
            X_val = X_val.values
        X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
        y_val_tensor = torch.tensor(y_val.values if isinstance(y_val, pd.Series) else y_val, dtype=torch.float32).reshape(-1, 1)
    
    # Load saved model or create a new one
    model = load_model(model_name, fold, input_dim)
    if model is None:
        print(f"Training new {model_name} model (fold {fold})")
        # Initialize model based on type
        if model_name == "DL_MLP":
            model = DL_MLP(input_dim)
        elif model_name == "DL_CNN":    
            model = DL_CNN(input_dim=input_dim)
        elif model_name == "DL_GRU":
            model = DL_GRU(input_dim=input_dim)
        else:
            raise ValueError(f"Unknown model type: {model_name}")
            
        criterion = nn.BCELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        
        # Create DataLoader with fixed random seed for worker
        generator = torch.Generator()
        generator.manual_seed(seed)
        
        if model_name == "DL_MLP":
            train_dataset = torch.utils.data.TensorDataset(X_train_tensor, y_train_tensor)
        else:  # DL_CNN or DL_GRU
            if seq_data is None:
                raise ValueError("seq_data is required for DL_CNN and DL_GRU")
            train_dataset = torch.utils.data.TensorDataset(
                X_train_tensor, 
                seq_data['del_seq'],
                seq_data['up5seq'],
                seq_data['down5seq'],
                y_train_tensor
            )
            
        # Use entire dataset as one batch
        batch_size = len(X_train_tensor)  # Use full dataset size
        train_loader = torch.utils.data.DataLoader(
            train_dataset, 
            batch_size=batch_size,  # Full batch training
            shuffle=True,
            generator=generator,
            worker_init_fn=lambda worker_id: np.random.seed(seed)
        )
        
        # Early stopping variables
        best_val_loss = float('inf')
        best_epoch = 0 if not is_final else epochs  # For final training, use total epochs
        no_improve = 0
        best_model_state = None
        
        # Train model
        model.train()
        for epoch in range(epochs):
            total_loss = 0
            batch_count = 0
            
            if model_name == "DL_MLP":
                for batch_X, batch_y in train_loader:
                    optimizer.zero_grad()
                    y_pred = model(batch_X)
                    loss = criterion(y_pred, batch_y)
                    loss.backward()
                    optimizer.step()
                    total_loss += loss.item()
                    batch_count += 1
            else:  # DL_CNN or DL_GRU
                for batch_X, batch_del, batch_up5, batch_down5, batch_y in train_loader:
                    optimizer.zero_grad()
                    y_pred = model(batch_X, batch_del, batch_up5, batch_down5)
                    loss = criterion(y_pred, batch_y)
                    loss.backward()
                    optimizer.step()
                    total_loss += loss.item()
                    batch_count += 1
            
            # Calculate validation loss if validation data is provided
            if val_data is not None and not is_final:
                model.eval()
                with torch.no_grad():
                    if model_name == "DL_MLP":
                        val_pred = model(X_val_tensor)
                        val_loss = criterion(val_pred, y_val_tensor).item()
                    else:  # DL_CNN or DL_GRU
                        val_pred = model(
                            X_val_tensor,
                            val_seq_data['del_seq'],
                            val_seq_data['up5seq'],
                            val_seq_data['down5seq']
                        )
                        val_loss = criterion(val_pred, y_val_tensor).item()
                model.train()
                
                # Early stopping check
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_epoch = epoch + 1
                    no_improve = 0
                    best_model_state = model.state_dict().copy()
                else:
                    no_improve += 1
                    
                if no_improve >= patience:
                    print(f"Early stopping at epoch {epoch+1}. Best epoch was {best_epoch}")
                    model.load_state_dict(best_model_state)
                    break
                    
            if (epoch + 1) % 10 == 0:  # Print every 10 epochs
                avg_loss = total_loss / batch_count
                status = f"Fold {fold}, Epoch {epoch+1}, Train Loss: {avg_loss:.3f}"
                if val_data is not None and not is_final:
                    status += f", Val Loss: {val_loss:.3f}"
                print(status)
        
        # Save model
        save_model(model, model_name, fold)
        
        # Remove storing best_epochs in defaultdict since we'll use metrics
        if not is_final and val_data is not None:
            pass  # We'll handle this in the main training loop
    else:
        print(f"Loaded existing {model_name} model (fold {fold})")
        best_epoch = None if not is_final else epochs  # For final training, use total epochs
    
    return model, best_epoch


# # =================================================
# # Prepare functions
# # =================================================
# import sys
# sys.path.append('/home/ha0214zh/github/lupsb/pon_del/code')
# from train_00_function import *