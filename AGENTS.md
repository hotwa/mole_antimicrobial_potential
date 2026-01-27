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
   - Model: `MolE-XGBoost-08.03.2024_14.20.pkl`

3. **Prediction Core** (`src/predictor.py`, `src/models.py`)
   - Lazy, concurrency-safe loading
   - Structured input normalization (Pydantic v2)
   - Returns pure items (no protocol metadata)

4. **Predictor Singleton** (`src/service.py`)
   - Single instance shared by FastAPI + FastMCP
   - Avoids duplicated model/data load

5. **FastMCP Server** (`mcp_server_enhanced.py`)
   - MCP JSON-RPC 2.0 over streamable HTTP
   - Tools: `status`, `predict_antimicrobial_potential`
   - Resources: `resource://schema`, `resource://about`, `resource://strains`

6. **FastAPI Server** (`api_server.py`)
   - REST endpoints `/health` and `/predict`
   - Uses the same predictor singleton

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

**Endpoint**: `POST http://<host>:8001/mcp`

**Protocol**: JSON-RPC 2.0 (FastMCP streamable HTTP). Initialize once and include the `mcp-session-id` header in follow-up calls.

**Initialize**:
```json
{
  "jsonrpc": "2.0",
  "id": "init",
  "method": "initialize",
  "params": {
    "protocolVersion": "2024-11-05",
    "capabilities": {},
    "clientInfo": {"name": "client", "version": "0.1.0"}
  }
}
```

**Then notify initialized**:
```json
{
  "jsonrpc": "2.0",
  "method": "notifications/initialized",
  "params": {}
}
```

**Tool call format**:
```json
{
  "jsonrpc": "2.0",
  "id": "request-id",
  "method": "tools/call",
  "params": {
    "name": "predict_antimicrobial_potential",
    "arguments": {
      "smiles": "CCO",
      "aggregate_scores": true
    }
  }
}
```

**Response**: JSON-RPC result; tool outputs are in `result.structuredContent`.

### FastAPI (MODE=api)

- `GET http://<host>:8000/health`
- `POST http://<host>:8000/predict`

Request body is the `MoleculeInput` schema (see `README_API.md`).

## Usage Examples

### For LLM Agents (MCP JSON-RPC)

#### 1) Single Molecule - Detailed Prediction
```json
{
  "jsonrpc": "2.0",
  "id": "1",
  "method": "tools/call",
  "params": {
    "name": "predict_antimicrobial_potential",
    "arguments": {
      "smiles": "CCO",
      "aggregate_scores": false
    }
  }
}
```

#### 2) Single Molecule - Aggregated Summary
```json
{
  "jsonrpc": "2.0",
  "id": "2",
  "method": "tools/call",
  "params": {
    "name": "predict_antimicrobial_potential",
    "arguments": {
      "smiles": "CCO",
      "aggregate_scores": true
    }
  }
}
```

#### 3) Batch Processing with Custom IDs
```json
{
  "jsonrpc": "2.0",
  "id": "3",
  "method": "tools/call",
  "params": {
    "name": "predict_antimicrobial_potential",
    "arguments": {
      "molecules": [
        {"smiles": "CCO", "chem_id": "ethanol"},
        {"smiles": "CCN", "chem_id": "ethylamine"}
      ],
      "aggregate_scores": true
    }
  }
}
```

#### 4) Custom Thresholds
```json
{
  "jsonrpc": "2.0",
  "id": "4",
  "method": "tools/call",
  "params": {
    "name": "predict_antimicrobial_potential",
    "arguments": {
      "smiles": "CCO",
      "aggregate_scores": true,
      "app_threshold": 0.1,
      "min_nkill": 5
    }
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
   - MCP errors are returned via JSON-RPC `error`
   - Check response `error` before reading `result`

4. **Performance Optimization**
   - Batch multiple molecules in single requests
   - Use aggregated results for high-throughput screening
   - Cache molecular representations if repeating queries

## Deployment

### Docker Compose
```bash
docker compose build
docker compose up -d
```

### Environment Configuration
```yaml
environment:
  - MODE=mcp_enhanced      # api or mcp_enhanced
  - HOST=0.0.0.0
  - PORT_API=8000          # FastAPI
  - PORT_MCP=8001          # FastMCP
  - MCP_TRANSPORT=http
```

### Manual Setup
```bash
conda env create -f environment.yaml
conda activate mole
python mcp_server_enhanced.py
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

### MCP Session Errors
- Ensure you call `initialize` and then `notifications/initialized`
- Include `mcp-session-id` in subsequent requests

## References

- **Paper**: [Pre-trained molecular representations enable antimicrobial discovery](https://www.nature.com/articles/s41467-025-58804-4)
- **MolE Framework**: [GitHub Repository](https://github.com/rolayoalarcon/MolE)
- **Training Data**: Maier et al. (2018), Nature

## Contact

For technical support or API questions, refer to the project README or contact the repository maintainer.
