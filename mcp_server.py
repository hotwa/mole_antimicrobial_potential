# mcp_server.py - 简化版本，仅保留MCP工具定义
import asyncio
from fastmcp import FastMCP
from predict_api import AntimicrobialPredictor, MoleculeInput

# Initialize predictor
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

    Args:
        input_data (MoleculeInput): The input data specifying molecules and options.
            - molecules (list of dict, optional): Structured input for complex batches
            - smiles (str or list of str, optional): One or more SMILES strings
            - chem_id (str or list of str, optional): Custom molecule IDs
            - aggregate_scores (bool, optional): Whether to aggregate results
            - app_threshold (float, optional): Probability threshold for growth inhibition
            - min_nkill (int, optional): Minimum number of strains inhibited for broad-spectrum

    Returns:
        List[dict]: Prediction results (detailed or aggregated based on aggregate_scores
    """
    return await predictor.predict(input_data)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--transport", default="http")
    args = parser.parse_args()
    
    print(f"[INFO] Starting MCP server on {args.host}:{args.port} via {args.transport}")
    mcp.run(host=args.host, port=args.port, transport=args.transport)
