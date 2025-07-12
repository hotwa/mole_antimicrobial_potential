import os
import re
import json
import pickle
import torch
import numpy as np
import pandas as pd
from typing import List, Union, Optional
from fastapi import FastAPI, Request, Response
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from scipy.stats.mstats import gmean
from sklearn.preprocessing import OneHotEncoder
from mole_representation import process_representation

class MoleculeInfo(BaseModel):
    smiles: str
    chem_id: Optional[str] = None

class MoleculeInput(BaseModel):
    molecules: Optional[List[MoleculeInfo]] = None
    smiles: Optional[Union[str, List[str]]] = None
    chem_id: Optional[Union[str, List[str]]] = None
    aggregate_scores: bool = Field(False, description="Whether to aggregate predictions")
    app_threshold: float = Field(0.04374140128493309, description="Threshold for growth inhibition")
    min_nkill: int = Field(10, description="Minimum strains for broad spectrum")

class AntimicrobialPredictor:
    """Service for antimicrobial activity prediction"""
    
    def __init__(self):
        self.device = "cuda:0" if torch.cuda.is_available() else "cpu"
        self.xgboost_model_path = "data/03.model_evaluation/MolE-XGBoost-08.03.2024_14.20.pkl"
        self.mole_model_path = "pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001"
        self.strain_categories = "data/01.prepare_training_data/maier_screening_results.tsv.gz"
        self.gram_information = "raw_data/maier_microbiome/strain_info_SF2.xlsx"
        
        # Load model once at initialization
        self.model = self._load_xgb_model()
        
        # Load strain data
        self.maier_screen = pd.read_csv(self.strain_categories, sep='\t', index_col=0)
        self.strain_ohe = self._prep_ohe(self.maier_screen.columns)
    
    def _load_xgb_model(self):
        """Load XGBoost model"""
        with open(self.xgboost_model_path, "rb") as file:
            return pickle.load(file)
    
    def _prep_ohe(self, categories):
        """Prepare one-hot encoding for strain variables"""
        ohe = OneHotEncoder(sparse=False)
        ohe.fit(pd.DataFrame(categories))
        return pd.DataFrame(
            ohe.transform(pd.DataFrame(categories)), 
            columns=categories, 
            index=categories
        )
    
    def _get_mole_representation(self, input_data: MoleculeInput):
        """Get MolE representation from SMILES, supports chem_id/molecules."""
        if input_data.molecules:
            smiles = [m.smiles for m in input_data.molecules]
            chem_ids = [m.chem_id if m.chem_id else f"mol{i+1}" for i, m in enumerate(input_data.molecules)]
        else:
            if isinstance(input_data.smiles, str):
                smiles = [input_data.smiles]
            else:
                smiles = input_data.smiles
            if input_data.chem_id:
                if isinstance(input_data.chem_id, str):
                    chem_ids = [input_data.chem_id]
                else:
                    chem_ids = input_data.chem_id
            else:
                chem_ids = [f"mol{i+1}" for i in range(len(smiles))]
        df = pd.DataFrame({"smiles": smiles, "chem_id": chem_ids})

        return process_representation(
            dataset_path=df,
            smile_column_str="smiles",
            id_column_str="chem_id",
            pretrained_dir=self.mole_model_path,
            device=self.device
        )

    
    def _add_strains(self, chemfeats_df):
        """Add strains to chemical features"""
        # Prepare chemical features
        chemfe = chemfeats_df.reset_index().rename(columns={"index": "chem_id"})
        chemfe["chem_id"] = chemfe["chem_id"].astype(str) 

        # Prepare OHE
        sohe = self.strain_ohe.reset_index().rename(columns={"index": "strain_name"})

        # Cartesian product merge
        xpred = chemfe.merge(sohe, how="cross")
        xpred["pred_id"] = xpred["chem_id"].str.cat(xpred["strain_name"], sep=":")

        xpred = xpred.set_index("pred_id")
        xpred = xpred.drop(columns=["chem_id", "strain_name"])
        
        return xpred
    
    def _gram_stain(self, label_df):
        """Add Gram stain information"""
        # Create copy of the label dataframe
        df_label = label_df.copy()
        
        # Read strain metadata
        maier_strains = pd.read_excel(
            self.gram_information,
            skiprows=[0,1, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54], 
            index_col="NT data base"
        )
        
        # Gather NT number
        df_label["nt_number"] = df_label["strain_name"].apply(
            lambda x: re.search(".*?\((NT\d+)\)", x).group(1)
        )

        # Create Gram strain dict
        gram_dict = maier_strains[["Gram stain"]].to_dict()["Gram stain"]

        # Add stain information
        df_label["gram_stain"] = df_label["nt_number"].apply(gram_dict.get)

        return df_label
    
    def _antimicrobial_potential(self, score_df):
        """Calculate antimicrobial potential scores"""
        # Separate chem_id from strain name
        score_df["chem_id"] = score_df["pred_id"].str.split(":", expand=True)[0]
        score_df["strain_name"] = score_df["pred_id"].str.split(":", expand=True)[1]

        # Add gram stain information
        pred_df = self._gram_stain(score_df)

        # Antimicrobial potential score
        apscore_total = pred_df.groupby("chem_id")["1"].apply(gmean).to_frame().rename(
            columns={"1": "apscore_total"}
        )
        apscore_total["apscore_total"] = np.log(apscore_total["apscore_total"])

        # Antimicrobial score by gram stain
        apscore_gram = pred_df.groupby(["chem_id", "gram_stain"])["1"].apply(gmean).unstack().rename(
            columns={"negative": "apscore_gnegative", "positive": "apscore_gpositive"}
        )
        apscore_gram["apscore_gnegative"] = np.log(apscore_gram["apscore_gnegative"])
        apscore_gram["apscore_gpositive"] = np.log(apscore_gram["apscore_gpositive"])
        
        # Number inhibited strains
        inhibted_total = pred_df.groupby("chem_id")["growth_inhibition"].sum().to_frame().rename(
            columns={"growth_inhibition": "ginhib_total"}
        )

        # Number inhibited strains per stain
        inhibted_gram = pred_df.groupby(["chem_id", "gram_stain"])["growth_inhibition"].sum().unstack().rename(
            columns={"negative": "ginhib_gnegative", "positive": "ginhib_gpositive"}
        )
        
        # Merge the results
        return apscore_total.join(apscore_gram).join(inhibted_total).join(inhibted_gram)
    
    async def predict(self, input_data: MoleculeInput):
        """
        预测小分子的抗菌潜力。（MCP工具标准接口，支持多种调用风格，详见下方用法示例）

        参数:
            input: MoleculeInput
                - smiles: str | List[str]    分子的SMILES字符串，或SMILES字符串列表
                - aggregate_scores: bool     是否聚合输出（True返回广谱/统计分数，False返回每个菌株概率，默认False）
                - chem_id: str | List[str]   可选，化合物自定义ID，和smiles一一对应
                - app_threshold: float       可选，抑制阈值，决定growth_inhibition
                - min_nkill: int             可选，定义广谱的“被抑制菌株数”下限
                - molecules: List[dict]      结构化输入，每个dict含smiles/chem_id（优先于上面参数）

        返回:
            PredictionResult
                - 若aggregate_scores=False: 返回每个分子的每个菌株预测概率/抑制情况
                - 若aggregate_scores=True:  返回每个分子的抗菌分数、抑制数、广谱判断

        用法示例:
            1. 单分子预测（默认菌株概率）
            >>> predict({"smiles": "CCO"})
            # 返回:
            [
                {
                    "pred_id": "mol1:Akkermansia muciniphila (NT5021)",
                    "antimicrobial_predictive_probability": 2.5e-6,
                    "growth_inhibition": 0
                },
                ...
            ]

            2. 单分子，聚合分数（只关心广谱结果）
            >>> predict({"smiles": "CCO", "aggregate_scores": True})
            # 返回:
            [
                {
                    "chem_id": "mol1",
                    "apscore_total": -11.7,
                    "apscore_gnegative": -11.6,
                    "apscore_gpositive": -11.8,
                    "ginhib_total": 0,
                    "ginhib_gnegative": 0,
                    "ginhib_gpositive": 0,
                    "broad_spectrum": 0
                }
            ]

            3. 多分子预测（SMILES数组）
            >>> predict({"smiles": ["CCO", "CCN"]})
            # 返回:
            [
                ...  # mol1/molecule1的每个菌株概率
                ...  # mol2/molecule2的每个菌株概率
            ]

            4. 多分子+聚合统计
            >>> predict({"smiles": ["CCO", "CCN"], "aggregate_scores": True})
            # 返回:
            [
                {"chem_id": "mol1", ...},
                {"chem_id": "mol2", ...}
            ]

            5. 自定义ID（chem_id与smiles一一对应）
            >>> predict({
                    "smiles": ["CCO", "CCN"],
                    "chem_id": ["ethanol", "ethylamine"],
                    "aggregate_scores": True
                })
            # 返回:
            [
                {"chem_id": "ethanol", ...},
                {"chem_id": "ethylamine", ...}
            ]

            6. 结构化molecules用法（推荐复杂场景）
            >>> predict({
                    "molecules": [
                        {"smiles": "CCO", "chem_id": "ethanol"},
                        {"smiles": "CCN", "chem_id": "ethylamine"}
                    ],
                    "aggregate_scores": True
                })
            # 返回:
            [
                {"chem_id": "ethanol", ...},
                {"chem_id": "ethylamine", ...}
            ]

            7. 自定义阈值和广谱下限
            >>> predict({
                    "smiles": "CCO",
                    "aggregate_scores": True,
                    "app_threshold": 0.1,
                    "min_nkill": 5
                })

        返回字段说明:
            - pred_id: 分子ID+菌株名称（如ethanol:Escherichia coli ...），仅aggregate_scores=False
            - antimicrobial_predictive_probability: 对每个菌株的抗菌预测概率
            - growth_inhibition: 0/1，是否被判定为抑制（概率高于阈值）
            - chem_id: 分子的自定义ID或mol1/mol2编号
            - apscore_total: 统计分数（越低越抗菌，需领域知识解读）
            - ginhib_total: 被抑制的菌株数
            - broad_spectrum: 0/1，是否“广谱抑制剂”

        注意:
            - 推荐聚合分数结果用于决策（aggregate_scores=True）。
            - 可混合用molecules结构体和smiles数组方式。
            - 不同分子/ID的输出顺序与输入一致。
            - MCP/LLM客户端可用工具发现和schema获取所有字段与默认值。
        """
        # Get MolE representation
        mole_representation = self._get_mole_representation(input_data)
        
        # Prepare strain-level predictions
        X_input = self._add_strains(mole_representation)
        
        # Make predictions
        y_pred = self.model.predict_proba(X_input)
        pred_df = pd.DataFrame(y_pred, columns=["0", "1"], index=X_input.index)
        
        # Binarize predictions using threshold
        pred_df["growth_inhibition"] = pred_df["1"].apply(
            lambda x: 1 if x >= input_data.app_threshold else 0
        )
        
        # Process results based on aggregation preference
        if input_data.aggregate_scores:
            pred_df = pred_df.reset_index()
            agg_df = self._antimicrobial_potential(pred_df)
            
            # Determine if chemical is broad spectrum
            agg_df["broad_spectrum"] = agg_df["ginhib_total"].apply(
                lambda x: 1 if x >= input_data.min_nkill else 0
            )
            
            # Convert to records for JSON serialization
            return agg_df.reset_index().to_dict(orient="records")
        else:
            # Report the antimicrobial predictive probabilities
            pred_df = pred_df.drop(columns=["0"])
            pred_df = pred_df.rename(columns={"1": "antimicrobial_predictive_probability"})
            
            # Convert to records for JSON serialization
            return pred_df.reset_index().to_dict(orient="records")

# Create FastAPI app
app = FastAPI(title="Antimicrobial Prediction MCP Tool")

# Initialize predictor service
predictor = AntimicrobialPredictor()

@app.post("/mcp")
async def mcp_endpoint(request: Request):
    """MCP endpoint that handles requests and streams responses using SSE"""
    # Parse request body
    try:
        body = await request.json()
    except Exception as e:
        return StreamingResponse(
            (f"event: error\ndata: {json.dumps({'error': str(e)})}\n\n",),
            media_type="text/event-stream"
        )
    # Extract MCP request data
    request_id = body.get("id", "unknown")
    method = body.get("method", "")
    params = body.get("params", {})
    
    # Prepare SSE response
    async def event_stream():
        print("DEBUG: method:", method, "type:", type(method))
        print("DEBUG: params:", params)
        try:
            if method == "predict":
                # Convert params to MoleculeInput model
                input_data = MoleculeInput(**params)
                print("DEBUG: input_data:", input_data)
                
                # Start event
                yield f"event: start\ndata: {json.dumps({'id': request_id})}\n\n"
                
                # Get prediction results
                results = await predictor.predict(input_data)
                
                # Send results
                yield f"event: data\ndata: {json.dumps({'id': request_id, 'result': results})}\n\n"
                
                # End event
                yield f"event: end\ndata: {json.dumps({'id': request_id})}\n\n"
            else:
                # Method not supported
                yield f"event: error\ndata: {json.dumps({'id': request_id, 'error': {'message': f'Method {method} not supported'}})}\n\n"
        except Exception as e:
            # Error event
            yield f"event: error\ndata: {json.dumps({'id': request_id, 'error': {'message': str(e)}})}\n\n"
    
    # Return streaming response
    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
        }
    )

@app.get("/mcp/schema")
async def mcp_schema():
    """Return the JSON schema for the MCP tool"""
    return {
        "openapi": {
            "info": {
                "title": "Antimicrobial Prediction Tool",
                "description": "Predicts antimicrobial activity from molecular SMILES",
                "version": "1.0.0"
            },
            "methods": {
                "predict": {
                    "description": "Predict antimicrobial activity for molecules",
                    "params": MoleculeInput.schema(),
                    "result": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "pred_id": {"type": "string"},
                                "antimicrobial_predictive_probability": {"type": "number"},
                                "growth_inhibition": {"type": "integer"}
                            }
                        }
                    }
                }
            }
        }
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)