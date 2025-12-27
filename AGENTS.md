# AGENTS.md - Antimicrobial Prediction System

## Project Overview

This project provides an AI-powered system for predicting the antimicrobial activity of small molecules based on their molecular structure. It uses a pre-trained MolE representation model combined with XGBoost classifiers to predict antibacterial activity against multiple bacterial strains.

### Key Capabilities

- **Molecular Representation**: Convert SMILES strings to high-dimensional molecular embeddings using MolE
- **Antimicrobial Prediction**: Predict growth inhibition probabilities for 100+ bacterial strains
- **Broad-Spectrum Detection**: Identify compounds with broad-spectrum antibacterial activity
- **Batch Processing**: Support for single or multiple molecule predictions
- **Customizable Thresholds**: Adjustable inhibition and broad-spectrum criteria

## System Architecture

### Core Components

1. **MolE Representation Model** (`pretrained_model/`)
   - Pre-trained graph neural network for molecular embeddings
   - Converts SMILES → 8000-dimensional vectors
   - Download required: [Zenodo link](https://zenodo.org/records/10803099)

2. **XGBoost Classifier** (`data/03.model_evaluation/`)
   - Trained on Maier et al. (2018) dataset
   - Predicts antimicrobial activity across 100+ strains
   - Models: `MolE-XGBoost-08.03.2024_14.20.pkl`

3. **Prediction Pipeline** (`mole_antimicrobial_prediction.py`)
   - Command-line interface for batch predictions
   - Supports SMILES input or pre-computed representations
   - Outputs detailed or aggregated results

4. **MCP/API Server** (`mcp_server.py`, `predict_api.py`)
   - FastMCP server for LLM agent integration
   - FastAPI REST endpoint with SSE streaming
   - JSON schema for tool discovery

### Data Flow

```
SMILES → MolE Model → XGBoost → Predictions → Aggregation
                    ↓
            Molecular Embeddings
                    ↓
         Strain-Specific Probabilities
                    ↓
         Antimicrobial Potential Scores
```

## API Endpoints

### MCP Server (Recommended for Agents)

**Endpoint**: `POST http://<host>:8000/mcp`

**Request Format**:
```json
{
  "id": "request-id",
  "method": "predict",
  "params": {
    "smiles": "CCO" | ["CCO", "CCN"],
    "chem_id": "ethanol" | ["ethanol", "ethylamine"],
    "molecules": [{"smiles": "CCO", "chem_id": "ethanol"}],
    "aggregate_scores": true,
    "app_threshold": 0.04374,
    "min_nkill": 10
  }
}
```

**Response Format** (SSE Stream):
```
event: start
data: {"id": "request-id"}

event: data
data: {"id": "request-id", "result": [...]}

event: end
data: {"id": "request-id"}
```

### REST API (Legacy)

**Endpoint**: `POST http://<host>:8000/mcp` (same as MCP)

**Alternative**: Direct function calls via `predict_api.py`

## Usage Examples

### For LLM Agents

#### 1. Single Molecule - Detailed Prediction
```json
{
  "method": "predict",
  "params": {
    "smiles": "CCO",
    "aggregate_scores": false
  }
}
```
**Returns**: 100+ predictions per strain with probabilities

#### 2. Single Molecule - Aggregated Summary
```json
{
  "method": "predict",
  "params": {
    "smiles": "CCO",
    "aggregate_scores": true
  }
}
```
**Returns**: Antimicrobial scores, inhibition counts, broad-spectrum flag

#### 3. Batch Processing with Custom IDs
```json
{
  "method": "predict",
  "params": {
    "molecules": [
      {"smiles": "CCO", "chem_id": "ethanol"},
      {"smiles": "CCN", "chem_id": "ethylamine"}
    ],
    "aggregate_scores": true
  }
}
```

#### 4. Custom Thresholds
```json
{
  "method": "predict",
  "params": {
    "smiles": "CCO",
    "aggregate_scores": true,
    "app_threshold": 0.1,
    "min_nkill": 5
  }
}
```

### Command Line

```bash
# Get molecular representations
python mole_representation.py input.tsv output.tsv --chemid_colname chem_id

# Predict antimicrobial activity
python mole_antimicrobial_prediction.py input.tsv output.tsv --smiles_input --chemid_colname chem_id

# Aggregate results
python mole_antimicrobial_prediction.py input.tsv output.tsv --smiles_input --aggregate_scores
```

## Response Field Reference

### Non-Aggregated Results (per strain)
- `pred_id`: "compound:strain" identifier
- `antimicrobial_predictive_probability`: 0-1 probability of inhibition
- `growth_inhibition`: 1 if probability ≥ threshold, else 0

### Aggregated Results (per molecule)
- `chem_id`: Molecule identifier
- `apscore_total`: Log geometric mean across all strains (lower = stronger)
- `apscore_gnegative`: Gram-negative specific score
- `apscore_gpositive`: Gram-positive specific score
- `ginhib_total`: Total strains inhibited
- `ginhib_gnegative`: Gram-negative strains inhibited
- `ginhib_gpositive`: Gram-positive strains inhibited
- `broad_spectrum`: 1 if ginhib_total ≥ min_nkill, else 0

## Interpretation Guide

### Antimicrobial Potential Scores
- **Lower values** indicate stronger antibacterial activity
- Scores are log-transformed geometric means
- Compare relative values between compounds

### Growth Inhibition
- Binary classification based on probability threshold
- Default: 0.04374 (from original publication)
- Customizable based on sensitivity requirements

### Broad-Spectrum Classification
- Compounds inhibiting ≥10 strains (default) are broad-spectrum
- Adjust `min_nkill` for different criteria
- Useful for identifying multi-target antibiotics

## Integration Guidelines for Agents

### Best Practices

1. **Input Validation**
   - Validate SMILES format before submission
   - Ensure chem_id uniqueness for batch processing
   - Use structured `molecules` input for complex scenarios

2. **Result Processing**
   - For drug discovery: Use aggregated results
   - For mechanism studies: Use detailed strain-level data
   - Consider both Gram-positive and Gram-negative activity

3. **Error Handling**
   - Check for "event: error" in SSE stream
   - Validate model files exist before deployment
   - Handle GPU/CPU fallback gracefully

4. **Performance Optimization**
   - Batch multiple molecules in single requests
   - Use aggregated results for high-throughput screening
   - Cache molecular representations if repeating queries

### Common Patterns

#### Drug Discovery Workflow
```json
{
  "method": "predict",
  "params": {
    "molecules": [{"smiles": "...", "chem_id": "candidate_1"}],
    "aggregate_scores": true,
    "min_nkill": 10
  }
}
```
→ Filter for broad-spectrum candidates

#### Lead Optimization
```json
{
  "method": "predict",
  "params": {
    "smiles": ["SMILES_A", "SMILES_B", "SMILES_C"],
    "aggregate_scores": false
  }
}
```
→ Compare strain-specific activity profiles

#### High-Throughput Screening
```json
{
  "method": "predict",
  "params": {
    "molecules": [...],  // 1000+ compounds
    "aggregate_scores": true,
    "app_threshold": 0.05
  }
}
```
→ Identify top candidates by broad-spectrum score

## Deployment

### Docker Compose
```bash
docker compose build
docker compose up -d
```

### Environment Configuration
```yaml
environment:
  - MODE=mcp              # mcp or api
  - HOST=0.0.0.0
  - PORT=8000
  - TRANSPORT=http        # http, sse, or stdio
  - WORKERS=1
```

### Manual Setup
```bash
conda env create -f environment.yaml
conda activate mole
python mcp_server.py --host 0.0.0.0 --port 8000 --transport http
```

## Troubleshooting

### Model Files Missing
- Download MolE model from Zenodo
- Place in `pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001/`
- Ensure `model.pth` and `config.yaml` exist

### GPU Memory Issues
- Set device to "cpu" in requests
- Process molecules in smaller batches
- Use aggregated results to reduce memory

### Invalid SMILES
- Validate SMILES format before submission
- Use RDKit or similar library for validation
- Check for special characters or encoding issues

### Performance Slow
- Batch molecules together
- Use aggregated results when possible
- Ensure GPU is available for MolE model

## References

- **Paper**: [Pre-trained molecular representations enable antimicrobial discovery](https://www.nature.com/articles/s41467-025-58804-4)
- **MolE Framework**: [GitHub Repository](https://github.com/rolayoalarcon/MolE)
- **Training Data**: Maier et al. (2018), Nature

## Contact

For technical support or API questions, refer to the project README or contact the repository maintainer.
