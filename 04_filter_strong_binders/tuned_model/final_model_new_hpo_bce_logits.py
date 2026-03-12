"""
Final Standalone Transformer+AAIndex Model for TCR-Peptide Binding Prediction
(New HPO Results - REVISED)

This model uses the best hyperparameters found during the REVISED hyperparameter optimization (HPO).
The REVISED HPO includes class ratio (pos_ratio) and loss weighting strategy as hyperparameters.
This variant sets the loss function to BCEWithLogitsLoss (logits output, no sigmoid in model).

Model Configuration (Hard-coded from NEW HPO results):
- Architecture: Transformer
- Encoding: AAIndex (5-dimensional)
- Best hyperparameters from NEW HPO (REVISED):
  * Learning Rate: 1.0224517407745427e-05
  * Batch Size: 64
  * Weight Decay: 5.618925735679741e-06
  * Dropout: 0.12507655751802993
  * Transformer Hidden Dimension: 256
  * Transformer Number of Layers: 2
  * Number of Attention Heads: 16
  * Feedforward Dimension: 256
  * Positive Ratio (pos_ratio): 0.5
  * Loss Strategy: sqrt_inv (square root of inverse frequency)
  * Max TCR Length: 50
  * Max Peptide Length: 30
  * Encoding Dimension: 5

Performance (from NEW HPO):
- Best Validation AUROC: 0.8215
- Final Validation AUROC: 0.7829
- Final Test AUROC: 0.4905
- Final Test AUPR: 0.6428

Comparison with Previous Model:
- Previous Model Validation AUROC: 0.829
- Previous Model Test AUROC: 0.509
- New Model Validation AUROC: 0.7829 (lower)
- New Model Test AUROC: 0.4905 (lower)
- Note: New model uses different class ratio and loss weighting strategy
"""

import os
import json
import numpy as np
import torch
import torch.nn as nn
from typing import Dict, List, Optional, Union, Tuple
from pathlib import Path


# AAIndex feature vectors (5-dimensional)
AAINDEX_PROPS = {
    'A': np.array([1.8, 0.0, 0.0, 89.1, 15.0], dtype=np.float32),
    'C': np.array([2.5, 0.0, 0.0, 121.2, 13.0], dtype=np.float32),
    'D': np.array([-3.5, 1.0, -1.0, 133.1, 59.0], dtype=np.float32),
    'E': np.array([-3.5, 1.0, -1.0, 147.1, 73.0], dtype=np.float32),
    'F': np.array([2.8, 0.0, 0.0, 165.2, 2.0], dtype=np.float32),
    'G': np.array([-0.4, 0.0, 0.0, 75.1, 66.0], dtype=np.float32),
    'H': np.array([-3.2, 1.0, 0.5, 155.2, 10.0], dtype=np.float32),
    'I': np.array([4.5, 0.0, 0.0, 131.2, 3.0], dtype=np.float32),
    'K': np.array([-3.9, 1.0, 1.0, 146.2, 100.0], dtype=np.float32),
    'L': np.array([3.8, 0.0, 0.0, 131.2, 5.0], dtype=np.float32),
    'M': np.array([1.9, 0.0, 0.0, 149.2, 1.0], dtype=np.float32),
    'N': np.array([-3.5, 1.0, 0.0, 132.1, 33.0], dtype=np.float32),
    'P': np.array([-1.6, 0.0, 0.0, 115.1, 95.0], dtype=np.float32),
    'Q': np.array([-3.5, 1.0, 0.0, 146.2, 33.0], dtype=np.float32),
    'R': np.array([-4.5, 1.0, 1.0, 174.2, 100.0], dtype=np.float32),
    'S': np.array([-0.8, 1.0, 0.0, 105.1, 95.0], dtype=np.float32),
    'T': np.array([-0.7, 1.0, 0.0, 119.1, 60.0], dtype=np.float32),
    'V': np.array([4.2, 0.0, 0.0, 117.1, 3.0], dtype=np.float32),
    'W': np.array([-0.9, 0.0, 0.0, 204.2, 5.0], dtype=np.float32),
    'Y': np.array([-1.3, 1.0, 0.0, 181.2, 2.0], dtype=np.float32)
}


# ============================================================================
# HARD-CODED HYPERPARAMETERS (from NEW HPO results - REVISED)
# ============================================================================
LEARNING_RATE = 1.0224517407745427e-05
BATCH_SIZE = 64
WEIGHT_DECAY = 5.618925735679741e-06
DROPOUT = 0.12507655751802993
TRANSFORMER_HIDDEN_DIM = 256
TRANSFORMER_NUM_LAYERS = 2
NHEAD = 16
DIM_FEEDFORWARD = 256
MAX_TCR_LEN = 50
MAX_PEPTIDE_LEN = 30
ENCODING_DIM = 5

# NEW HPO parameters (training-time, not model architecture)
POS_RATIO = 0.5  # Class ratio for training data resampling
LOSS_STRATEGY = "bce_logits"  # BCEWithLogitsLoss (model outputs logits)


def get_loss() -> nn.Module:
    """Return the configured BCE-with-logits loss."""
    return nn.BCEWithLogitsLoss()


class AAIndexEncoder:
    """AAIndex feature vectors encoder for amino acid sequences."""
    
    def __init__(self):
        self.encoding_dim = 5
        self.props = AAINDEX_PROPS.copy()
        all_vals = np.array(list(self.props.values()))
        self.min_vals = all_vals.min(axis=0)
        self.max_vals = all_vals.max(axis=0)
        self.range_vals = self.max_vals - self.min_vals + 1e-8
    
    def encode(self, sequence: str, max_length: int) -> np.ndarray:
        """Encode a sequence to AAIndex representation.
        
        Args:
            sequence: Amino acid sequence string
            max_length: Maximum sequence length (padded/truncated)
        
        Returns:
            Encoded sequence array of shape (max_length, 5)
        """
        seq = str(sequence).upper().replace('|', '').strip()
        encoded = np.zeros((max_length, 5), dtype=np.float32)
        for i, aa in enumerate(seq[:max_length]):
            if aa in self.props:
                encoded[i] = (self.props[aa] - self.min_vals) / self.range_vals
        return encoded
    
    def create_mask(self, sequence: str, max_length: int) -> np.ndarray:
        """Create padding mask: 1 for padding positions, 0 for valid.
        
        Args:
            sequence: Amino acid sequence string
            max_length: Maximum sequence length
        
        Returns:
            Mask array of shape (max_length,) where 1 = padding, 0 = valid
        """
        seq = str(sequence).strip()
        seq_len = len(seq) if seq else 0
        mask = np.zeros(max_length, dtype=np.float32)
        if seq_len < max_length:
            mask[seq_len:] = 1.0
        return mask


class TransformerModel(nn.Module):
    """Transformer model for TCR-peptide binding prediction."""
    
    def __init__(self, encoding_dim: int = 5, max_tcr_len: int = 50, max_peptide_len: int = 30,
                 transformer_hidden_dim: int = 256, transformer_num_layers: int = 2, 
                 dropout: float = 0.125, dim_feedforward: int = 256, nhead: int = 16,
                 delta_feature_dim: int = 3):
        super().__init__()
        self.encoding_dim = encoding_dim
        d_model = transformer_hidden_dim
        self.delta_feature_dim = delta_feature_dim
        
        # Positional encoding
        self.tcr_pos_embed = nn.Embedding(max_tcr_len, d_model)
        self.pep_pos_embed = nn.Embedding(max_peptide_len, d_model)
        
        # Project encoding to d_model
        self.tcr_proj = nn.Linear(encoding_dim, d_model)
        self.pep_proj = nn.Linear(encoding_dim, d_model)
        
        # Transformer encoders
        encoder_layer = nn.TransformerEncoderLayer(
            d_model, nhead, dim_feedforward=dim_feedforward, 
            dropout=dropout, batch_first=True
        )
        self.tcr_transformer = nn.TransformerEncoder(encoder_layer, transformer_num_layers)
        self.pep_transformer = nn.TransformerEncoder(encoder_layer, transformer_num_layers)
        
        # Combined classifier
        self.fc1 = nn.Linear(d_model * 2, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64 + self.delta_feature_dim, 1)
        self.dropout = nn.Dropout(dropout)
        self.relu = nn.ReLU()
    
    def forward(self, tcr, peptide, tcr_mask=None, peptide_mask=None, delta_features=None):
        """Forward pass.
        
        Args:
            tcr: TCR sequences (batch, seq_len, encoding_dim)
            peptide: Peptide sequences (batch, seq_len, encoding_dim)
            tcr_mask: TCR padding mask (batch, seq_len), 1 = padding, 0 = valid
            peptide_mask: Peptide padding mask (batch, seq_len), 1 = padding, 0 = valid
            delta_features: Optional dissimilarity features (batch, delta_feature_dim)
        
        Returns:
            Binding probability (batch,)
        """
        batch_size = tcr.size(0)
        
        # Project to d_model
        tcr = self.tcr_proj(tcr)
        peptide = self.pep_proj(peptide)
        
        # Add positional encoding
        tcr_pos = torch.arange(tcr.size(1), device=tcr.device).unsqueeze(0).expand(batch_size, -1)
        tcr = tcr + self.tcr_pos_embed(tcr_pos)
        
        pep_pos = torch.arange(peptide.size(1), device=peptide.device).unsqueeze(0).expand(batch_size, -1)
        peptide = peptide + self.pep_pos_embed(pep_pos)
        
        # Transformer encoding
        # Convert mask to bool: True = padding (to be ignored), False = valid
        tcr_key_padding_mask = tcr_mask.bool() if tcr_mask is not None else None
        pep_key_padding_mask = peptide_mask.bool() if peptide_mask is not None else None
        
        tcr_out = self.tcr_transformer(tcr, src_key_padding_mask=tcr_key_padding_mask)
        pep_out = self.pep_transformer(peptide, src_key_padding_mask=pep_key_padding_mask)
        
        # Global average pooling with mask
        if tcr_mask is not None:
            tcr_mask_expanded = tcr_mask.unsqueeze(-1).expand_as(tcr_out)
            tcr_out = tcr_out.masked_fill(tcr_mask_expanded.bool(), 0)
            tcr_final = tcr_out.sum(dim=1) / (~tcr_mask.bool()).sum(dim=1, keepdim=True).float()
        else:
            tcr_final = tcr_out.mean(dim=1)
        
        if peptide_mask is not None:
            pep_mask_expanded = peptide_mask.unsqueeze(-1).expand_as(pep_out)
            pep_out = pep_out.masked_fill(pep_mask_expanded.bool(), 0)
            pep_final = pep_out.sum(dim=1) / (~peptide_mask.bool()).sum(dim=1, keepdim=True).float()
        else:
            pep_final = pep_out.mean(dim=1)
        
        # Concatenate
        combined = torch.cat([tcr_final, pep_final], dim=1)
        
        # Classifier
        x = self.relu(self.fc1(combined))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        if delta_features is None:
            delta_features = torch.zeros(
                (batch_size, self.delta_feature_dim), device=x.device, dtype=x.dtype
            )
        x = torch.cat([x, delta_features], dim=1)
        x = self.fc3(x)
        
        return x.squeeze(-1)


class TCRPeptidePredictor:
    """Complete predictor class for TCR-peptide binding prediction (NEW HPO)."""
    
    def __init__(self, hpo_results_path: Optional[str] = None, model_path: Optional[str] = None,
                 device: Optional[str] = None):
        """Initialize the predictor with NEW HPO hyperparameters.
        
        Args:
            hpo_results_path: Path to new_hpo40_results.json file (optional, uses hardcoded values by default)
            model_path: Path to saved model weights (.pt file). If None, model will not be loaded.
            device: Device to use ('cuda', 'cpu', or None for auto-detection)
        """
        # Set device
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = torch.device(device)
        
        # Use hardcoded hyperparameters from NEW HPO (default)
        self.best_params = {
            'learning_rate': LEARNING_RATE,
            'batch_size': BATCH_SIZE,
            'weight_decay': WEIGHT_DECAY,
            'dropout': DROPOUT,
            'transformer_hidden_dim': TRANSFORMER_HIDDEN_DIM,
            'transformer_num_layers': TRANSFORMER_NUM_LAYERS,
            'nhead': NHEAD,
            'dim_feedforward': DIM_FEEDFORWARD,
            'pos_ratio': POS_RATIO,
            'loss_strategy': LOSS_STRATEGY
        }
        
        # Override with HPO results if provided
        if hpo_results_path is not None and os.path.exists(hpo_results_path):
            with open(hpo_results_path, 'r') as f:
                hpo_results = json.load(f)
            self.hpo_config = hpo_results.get('config', {})
            best_trial_params = hpo_results.get('best_trial', {}).get('params', {})
            if best_trial_params:
                self.best_params.update(best_trial_params)
            self.performance = hpo_results.get('final_model', {}).get('test_metrics', {})
        else:
            self.hpo_config = {
                'max_tcr_len': MAX_TCR_LEN,
                'max_peptide_len': MAX_PEPTIDE_LEN,
                'encoding_dim': ENCODING_DIM
            }
            self.performance = {
                'auroc': 0.4905,
                'aupr': 0.6428,
                'accuracy': 0.3041,
                'precision': 0.3167,
                'recall': 0.0030,
                'f1': 0.0059
            }
        
        # Extract hyperparameters
        self.max_tcr_len = self.hpo_config.get('max_tcr_len', MAX_TCR_LEN)
        self.max_peptide_len = self.hpo_config.get('max_peptide_len', MAX_PEPTIDE_LEN)
        self.encoding_dim = self.hpo_config.get('encoding_dim', ENCODING_DIM)
        
        # Initialize encoder
        self.encoder = AAIndexEncoder()
        
        # Initialize model with NEW HPO hyperparameters
        self.model = TransformerModel(
            encoding_dim=self.encoding_dim,
            max_tcr_len=self.max_tcr_len,
            max_peptide_len=self.max_peptide_len,
            transformer_hidden_dim=self.best_params.get('transformer_hidden_dim', TRANSFORMER_HIDDEN_DIM),
            transformer_num_layers=self.best_params.get('transformer_num_layers', TRANSFORMER_NUM_LAYERS),
            dropout=self.best_params.get('dropout', DROPOUT),
            dim_feedforward=self.best_params.get('dim_feedforward', DIM_FEEDFORWARD),
            nhead=self.best_params.get('nhead', NHEAD)
        )
        
        self.model.to(self.device)
        
        # Load model weights if provided
        if model_path is not None and os.path.exists(model_path):
            self.load_model(model_path)
            print(f"Model weights loaded from {model_path}")
        else:
            print("No model weights loaded. Model is initialized with random weights.")
            if model_path is not None:
                print(f"Warning: Model path {model_path} does not exist.")
    
    def encode_sequences(self, tcr_seq: str, peptide_seq: str) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Encode TCR and peptide sequences.
        
        Args:
            tcr_seq: TCR CDR3β sequence
            peptide_seq: Peptide sequence
        
        Returns:
            Tuple of (tcr_encoded, peptide_encoded, tcr_mask, peptide_mask) as tensors
        """
        tcr_encoded = self.encoder.encode(tcr_seq, self.max_tcr_len)
        peptide_encoded = self.encoder.encode(peptide_seq, self.max_peptide_len)
        
        tcr_mask = self.encoder.create_mask(tcr_seq, self.max_tcr_len)
        peptide_mask = self.encoder.create_mask(peptide_seq, self.max_peptide_len)
        
        tcr_tensor = torch.from_numpy(tcr_encoded).float().unsqueeze(0).to(self.device)
        peptide_tensor = torch.from_numpy(peptide_encoded).float().unsqueeze(0).to(self.device)
        tcr_mask_tensor = torch.from_numpy(tcr_mask).float().unsqueeze(0).to(self.device)
        peptide_mask_tensor = torch.from_numpy(peptide_mask).float().unsqueeze(0).to(self.device)
        
        return tcr_tensor, peptide_tensor, tcr_mask_tensor, peptide_mask_tensor
    
    def predict(self, tcr_seq: str, peptide_seq: str, delta_features: Optional[np.ndarray] = None) -> float:
        """Predict binding probability for a single TCR-peptide pair.
        
        Args:
            tcr_seq: TCR CDR3β sequence
            peptide_seq: Peptide sequence
            delta_features: Optional array-like with 3 dissimilarity features
        
        Returns:
            Binding probability (0-1)
        """
        self.model.eval()
        with torch.no_grad():
            tcr, peptide, tcr_mask, peptide_mask = self.encode_sequences(tcr_seq, peptide_seq)
            delta_tensor = None
            if delta_features is not None:
                delta_tensor = torch.tensor(delta_features, dtype=torch.float32, device=self.device).unsqueeze(0)
            logits = self.model(tcr, peptide, tcr_mask, peptide_mask, delta_tensor)
            prob = torch.sigmoid(logits)
            return prob.item()
    
    def predict_batch(self, tcr_seqs: List[str], peptide_seqs: List[str], 
                     batch_size: int = 64, delta_features: Optional[np.ndarray] = None) -> np.ndarray:
        """Predict binding probabilities for a batch of TCR-peptide pairs.
        
        Args:
            tcr_seqs: List of TCR CDR3β sequences
            peptide_seqs: List of peptide sequences
            batch_size: Batch size for prediction (default: 64 from NEW HPO)
            delta_features: Optional array-like (N, 3) of dissimilarity features
        
        Returns:
            Array of binding probabilities
        """
        if len(tcr_seqs) != len(peptide_seqs):
            raise ValueError("tcr_seqs and peptide_seqs must have the same length")
        
        self.model.eval()
        all_probs = []
        
        with torch.no_grad():
            for i in range(0, len(tcr_seqs), batch_size):
                batch_tcr_seqs = tcr_seqs[i:i+batch_size]
                batch_peptide_seqs = peptide_seqs[i:i+batch_size]
                
                # Encode batch
                batch_tcr = []
                batch_peptide = []
                batch_tcr_mask = []
                batch_peptide_mask = []
                
                for tcr_seq, pep_seq in zip(batch_tcr_seqs, batch_peptide_seqs):
                    tcr_enc = self.encoder.encode(tcr_seq, self.max_tcr_len)
                    pep_enc = self.encoder.encode(pep_seq, self.max_peptide_len)
                    tcr_m = self.encoder.create_mask(tcr_seq, self.max_tcr_len)
                    pep_m = self.encoder.create_mask(pep_seq, self.max_peptide_len)
                    
                    batch_tcr.append(tcr_enc)
                    batch_peptide.append(pep_enc)
                    batch_tcr_mask.append(tcr_m)
                    batch_peptide_mask.append(pep_m)
                
                # Convert to tensors
                tcr_tensor = torch.from_numpy(np.array(batch_tcr)).float().to(self.device)
                peptide_tensor = torch.from_numpy(np.array(batch_peptide)).float().to(self.device)
                tcr_mask_tensor = torch.from_numpy(np.array(batch_tcr_mask)).float().to(self.device)
                peptide_mask_tensor = torch.from_numpy(np.array(batch_peptide_mask)).float().to(self.device)
                
                # Predict
                delta_tensor = None
                if delta_features is not None:
                    delta_batch = delta_features[i:i+batch_size]
                    delta_tensor = torch.tensor(delta_batch, dtype=torch.float32, device=self.device)
                logits = self.model(tcr_tensor, peptide_tensor, tcr_mask_tensor, peptide_mask_tensor, delta_tensor)
                probs = torch.sigmoid(logits)
                all_probs.extend(probs.cpu().numpy())
        
        return np.array(all_probs)
    
    def save_model(self, save_path: str):
        """Save model weights to a file.
        
        Args:
            save_path: Path to save the model (.pt file)
        """
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'hyperparameters': self.best_params,
            'config': {
                'max_tcr_len': self.max_tcr_len,
                'max_peptide_len': self.max_peptide_len,
                'encoding_dim': self.encoding_dim,
                'transformer_hidden_dim': self.best_params.get('transformer_hidden_dim', TRANSFORMER_HIDDEN_DIM),
                'transformer_num_layers': self.best_params.get('transformer_num_layers', TRANSFORMER_NUM_LAYERS),
                'dropout': self.best_params.get('dropout', DROPOUT),
                'dim_feedforward': self.best_params.get('dim_feedforward', DIM_FEEDFORWARD),
                'nhead': self.best_params.get('nhead', NHEAD),
                'delta_feature_dim': self.model.delta_feature_dim,
                'pos_ratio': self.best_params.get('pos_ratio', POS_RATIO),
                'loss_strategy': self.best_params.get('loss_strategy', LOSS_STRATEGY)
            },
            'performance': self.performance
        }, save_path)
        print(f"Model saved to {save_path}")
    
    def load_model(self, model_path: str):
        """Load model weights from a file.
        
        Args:
            model_path: Path to the saved model (.pt file)
        """
        checkpoint = torch.load(model_path, map_location=self.device)
        state_dict = checkpoint.get('model_state_dict', checkpoint)
        model_state = self.model.state_dict()
        compatible_state = {
            k: v for k, v in state_dict.items()
            if k in model_state and v.shape == model_state[k].shape
        }
        missing_keys = [k for k in model_state.keys() if k not in compatible_state]
        self.model.load_state_dict(compatible_state, strict=False)
        print(f"Model loaded from {model_path}")
        if missing_keys:
            print(f"Warning: Skipped loading {len(missing_keys)} keys due to size mismatch or absence.")
    
    def get_model_info(self) -> Dict:
        """Get model information including hyperparameters and performance.
        
        Returns:
            Dictionary with model information
        """
        return {
            'hyperparameters': self.best_params,
            'config': {
                'max_tcr_len': self.max_tcr_len,
                'max_peptide_len': self.max_peptide_len,
                'encoding_dim': self.encoding_dim,
                'delta_feature_dim': self.model.delta_feature_dim
            },
            'performance': self.performance,
            'device': str(self.device),
            'model_parameters': sum(p.numel() for p in self.model.parameters())
        }


# Example usage
if __name__ == '__main__':
    # Example 1: Initialize predictor with NEW HPO hyperparameters
    print("Initializing predictor with NEW HPO hyperparameters...")
    predictor = TCRPeptidePredictor(
        hpo_results_path='new_hpo_results/new_hpo40_results.json',
        device='cuda' if torch.cuda.is_available() else 'cpu'
    )
    
    # Print model information
    print("\nModel Information (NEW HPO):")
    info = predictor.get_model_info()
    print(f"  Device: {info['device']}")
    print(f"  Parameters: {info['model_parameters']:,}")
    print(f"  Max TCR Length: {info['config']['max_tcr_len']}")
    print(f"  Max Peptide Length: {info['config']['max_peptide_len']}")
    print(f"  Encoding Dimension: {info['config']['encoding_dim']}")
    print(f"\nHyperparameters (NEW HPO):")
    for key, value in info['hyperparameters'].items():
        print(f"  {key}: {value}")
    if info['performance']:
        print(f"\nPerformance (Test Set):")
        for key, value in info['performance'].items():
            print(f"  {key}: {value:.4f}")
    
    # Comparison with previous model
    print("\n" + "="*50)
    print("Comparison with Previous Model:")
    print("="*50)
    print("Previous Model (Old HPO):")
    print("  Validation AUROC: 0.829")
    print("  Test AUROC: 0.509")
    print("\nNew Model (REVISED HPO):")
    print(f"  Validation AUROC: 0.7829")
    print(f"  Test AUROC: {info['performance'].get('auroc', 0.4905):.4f}")
    print("\nNote: New model uses pos_ratio=0.5 and loss_strategy='sqrt_inv'")
    
    # Example 2: Single prediction
    print("\n" + "="*50)
    print("Example: Single Prediction")
    print("="*50)
    tcr_seq = "CASSQGVGPLYEQYF"
    peptide_seq = "NLVPMVATV"
    prob = predictor.predict(tcr_seq, peptide_seq)
    print(f"TCR: {tcr_seq}")
    print(f"Peptide: {peptide_seq}")
    print(f"Binding Probability: {prob:.4f}")
    
    # Example 3: Batch prediction
    print("\n" + "="*50)
    print("Example: Batch Prediction")
    print("="*50)
    tcr_seqs = ["CASSQGVGPLYEQYF", "CASSDVGPLYEQYF", "CASSQGVGPLYEQYF"]
    peptide_seqs = ["NLVPMVATV", "NLVPMVATV", "ILKEPVHGV"]
    probs = predictor.predict_batch(tcr_seqs, peptide_seqs)
    for tcr, pep, prob in zip(tcr_seqs, peptide_seqs, probs):
        print(f"TCR: {tcr}, Peptide: {pep}, Probability: {prob:.4f}")
    
    # Example 4: Save model
    print("\n" + "="*50)
    print("Example: Save Model")
    print("="*50)
    save_path = 'final_model_new_hpo_weights.pt'
    predictor.save_model(save_path)
    print(f"Model saved to {save_path}")
    
    print("\n" + "="*50)
    print("Example Usage Complete!")
    print("="*50)
