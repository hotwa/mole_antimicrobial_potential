# mcp_server.py
import argparse
from fastmcp import FastMCP
from predict_api import AntimicrobialPredictor, MoleculeInput

predictor = AntimicrobialPredictor()
mcp = FastMCP("Antimicrobial_MCP_Server")

@mcp.tool
async def predict_antimicrobial_activity(input_data: MoleculeInput) -> list:
    """
    Predict the antimicrobial potential of molecules.
    
    This function predicts the antibacterial activity for one or more small molecules.
    Supports automatic molecule ID matching, batch prediction, result aggregation, 
    and custom molecule IDs. It is suitable for LLM agents, bio-pharmaceutical R&D, 
    and automated high-throughput screening.

    This function returns prediction results immediately. 
    Batch inputs and custom identifiers are supported.

    Args:
        input_data (MoleculeInput): The input data specifying molecules and options.
            - molecules (list of dict, optional): Structured input for complex batches or when custom IDs are required.
                Each dictionary should include:
                    - smiles (str, required): The SMILES string representing the molecule.
                    - chem_id (str, optional): A custom identifier for the molecule (e.g., "ethanol"). 
                      If not provided, an automatic ID is generated.
            - smiles (str or list of str, optional): One or more SMILES strings. When using a list, 
              it should correspond to the chem_id list if provided. 
              If both molecules and smiles are provided, molecules takes precedence.
            - chem_id (str or list of str, optional): Custom molecule IDs, matching the order of SMILES if provided.
            - aggregate_scores (bool, optional): Whether to aggregate results. 
              True returns summary statistics for each molecule; 
              False returns predicted probabilities for each molecule-strain pair. (Default: False)
            - app_threshold (float, optional): Probability threshold for growth inhibition 
              (predicted probability >= threshold is considered inhibition; default 0.04374).
            - min_nkill (int, optional): Minimum number of strains inhibited to define broad-spectrum antibacterial 
              (default: 10).

    Example Usage:
        # Single molecule prediction (detailed by strain)
        predict(
            input_data={
                "smiles": "CCO"
            }
        )

        # Single molecule, aggregated result
        predict(
            input_data={
                "smiles": "CCO",
                "aggregate_scores": True
            }
        )

        # Multiple molecules, batch prediction
        predict(
            input_data={
                "smiles": ["CCO", "CCN"]
            }
        )

        # Multiple molecules, aggregation, and custom IDs
        predict(
            input_data={
                "smiles": ["CCO", "CCN"],
                "chem_id": ["ethanol", "ethylamine"],
                "aggregate_scores": True
            }
        )

        # Structured input (recommended for complex batch with custom IDs)
        predict(
            input_data={
                "molecules": [
                    {"smiles": "CCO", "chem_id": "ethanol"},
                    {"smiles": "CCN", "chem_id": "ethylamine"}
                ],
                "aggregate_scores": True
            }
        )

        # Custom inhibition threshold and broad-spectrum criteria
        predict(
            input_data={
                "smiles": "CCO",
                "aggregate_scores": True,
                "app_threshold": 0.1,
                "min_nkill": 5
            }
        )

    Note on JSON/Dictionary Inputs:
        When using agents or tools that require string-encoded JSON input 
        (such as certain LLM toolchains), ensure that dictionaries are 
        properly serialized and escaped. For example:
            input_data="{\"smiles\": [\"CCO\", \"CCN\"], \"aggregate_scores\": true}"

        If using structured data, make sure all nested structures are properly escaped.
        For REST APIs, direct JSON objects are supported.

    Returns:
        If aggregate_scores=False:
            List[dict]: Each element contains the prediction for a molecule-strain pair:
                [
                    {
                        "pred_id": "ethanol:Escherichia coli (NT5077)",
                        "antimicrobial_predictive_probability": 0.00005,
                        "growth_inhibition": 0
                    },
                    ...
                ]
        If aggregate_scores=True:
            List[dict]: Each element summarizes antibacterial scores and broad-spectrum properties for a molecule:
                [
                    {
                        "chem_id": "ethanol",
                        "apscore_total": -11.7,
                        "apscore_gnegative": -11.6,
                        "apscore_gpositive": -11.8,
                        "ginhib_total": 0,
                        "ginhib_gnegative": 0,
                        "ginhib_gpositive": 0,
                        "broad_spectrum": 0
                    },
                    ...
                ]

    Field Descriptions:
        - pred_id: Combination of molecule ID and strain (for non-aggregated results).
        - antimicrobial_predictive_probability: Predicted inhibition probability (0~1).
        - growth_inhibition: Indicates effective inhibition (1 = yes, 0 = no, controlled by app_threshold).
        - chem_id: Custom or auto-generated molecule ID.
        - apscore_total: Geometric mean score across all strains (lower = stronger antibacterial potential).
        - apscore_gnegative / apscore_gpositive: Scores for Gram-negative/positive strains.
        - ginhib_total: Total number of strains effectively inhibited.
        - ginhib_gnegative / ginhib_gpositive: Number of Gram-negative/positive strains inhibited.
        - broad_spectrum: Indicates broad-spectrum inhibition (1 = yes, 0 = no, per min_nkill criteria).

    Recommendations:
        - For high-throughput screening, use batch input (multiple SMILES or structured molecules).
        - Use aggregation for summary statistics, especially in drug discovery workflows.
        - Custom chem_id is recommended for database integration or complex queries.

    LLM/Agent Usage Notes:
        - Field types and sample inputs are auto-detected by compatible agents.
        - Supports OpenAPI, JSON Schema, Claude, and other major agent platforms.
        - For agent platforms that require stringified JSON input, ensure proper escaping as in the examples above.

    Documentation:
        - See project README, /mcp/schema endpoint, or llms.txt protocol for full API schema and field descriptions.
    """
    return await predictor.predict(input_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--transport", default="http")  # You can also use http, stdio
    args = parser.parse_args()
    mcp.run(host=args.host, port=args.port, transport=args.transport)
