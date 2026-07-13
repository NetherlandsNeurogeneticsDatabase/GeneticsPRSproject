#!/usr/bin/env python3
"""
================================================================================
PURPOSE:
    Parses an input file containing sentence data, checks for structural regex 
    matches, and routes hits through an LLM feature extraction routine multiple 
    times. It then calculates consensus votes/scores for target features and writes 
    the completed evaluations back to disk.


DEPENDENCIES (NON-PIP / LOCAL LIBRARY):
    * from scripts.utils import extract_features_from_text
      -> (Custom internal project function; must reside within your local workspace directory).

INPUT:
    * A JSON configuration configuration file containing:
        - "input_file": Path to data table (TSV/CSV).
        - "output_file": Path to save evaluated target columns.
        - "feature_cols": List of target feature column names to extract.
        - "id_cols": Unique key identifier headers to retain.
        - "patterns": List of regular expressions used to pre-screen row relevance.
        - "prompt_file": Path containing the LLM system prompt framework.
        - "n_runs", "token", "temp_val", "model_name" (Optional LLM tuning parameters).

PROCESS:
    1. Parse text-based parameters out of the targeted configuration JSON.
    2. Scan target sentences sequentially against regex rules to bypass empty elements.
    3. Run hit text structures through iterative generative model parsing blocks.
    4. Enforce consensus voting using a Collection Counter tracking matrix.

OUTPUT:
    * Writes out a structural duplicate of the target file containing identical 
      identifying ID columns accompanied by freshly extracted and scored features.
================================================================================
"""

import argparse
import json
import re
import pandas as pd
from collections import Counter
from tqdm import tqdm

# ==============================================================================
# NON-PIP / LOCAL PROJECT DEPENDENCY
# ==============================================================================
from utils import extract_features_from_text


def main(config_path):
    # --------------------------------------------------------------------------
    # 1. SETUP & CONFIGURATION LOADING
    # --------------------------------------------------------------------------
    with open(config_path, "r") as f:
        config = json.load(f)

    # Read data file based on structural definitions
    df = pd.read_csv(config["input_file"], sep="\t" if config.get("tsv", True) else ",")
    feature_cols = config["feature_cols"]
    id_cols = config["id_cols"]

    # Compile evaluation criteria patterns
    patterns = [re.compile(p, re.IGNORECASE | re.VERBOSE) for p in config["patterns"]]
    results = []

    # --------------------------------------------------------------------------
    # 2. SENTENCE PRE-SCREENING & FEATURE EXTRACTION LOOP
    # --------------------------------------------------------------------------
    try:
        for _, row in tqdm(df.iterrows(), total=len(df)):
            sentence = str(row.get("Sentence", "") or "")
            has_match = any(p.search(sentence) for p in patterns)

            # If no regex match is found, preserve structure but leave feature blank
            if not has_match:
                pred_dict = {col: None for col in feature_cols}
                for col in feature_cols:
                    pred_dict[f"{col}_score"] = None
            else:
                preds = []
                
                # UNCOMMENT TO MONITOR CURRENT PROCESSING SENTENCE:
                # print(f"{row['NBB_nr']} with {row['Sentence']}")
                
                # Execute evaluation iterations to establish distribution density
                for _ in range(config.get("n_runs", 10)):
                    p = extract_features_from_text(
                        sentence,
                        config["prompt_file"],
                        token=config.get("token", 20),
                        temp_val=config.get("temp_val", 0.1),
                        model_name=config.get("model_name", "meta-llama/Llama-3.1-8B-Instruct"),
                    )
                    preds.append(p)

                # --------------------------------------------------------------
                # 3. CONSENSUS AGGREGATION & SCORING
                # --------------------------------------------------------------
                pred_dict = {}
                for col in feature_cols:
                    # UNCOMMENT TO RUN TYPE VALIDATION CHECKS ON RAW LLM PARSER OUTPUTS:
                    # print("DEBUG extract_features_from_text output type:", type(p), col, "value:")
                    
                    values = [p.get(col) for p in preds if p and p.get(col) is not None]
                    
                    if values:
                        most_common, count = Counter(values).most_common(1)[0]
                        pred_dict[col] = most_common
                        pred_dict[f"{col}_score"] = round(count / len(preds), 1)
                    else:
                        pred_dict[col] = None
                        pred_dict[f"{col}_score"] = None

            # Map the unique index data keys side-by-side with calculated values
            results.append({**{col: row.get(col) for col in id_cols}, **pred_dict})

    except KeyboardInterrupt:
        print("\Interrupted by user — saving partial results...")

    # --------------------------------------------------------------------------
    # 4. EXPORT OUTPUT DATA
    # --------------------------------------------------------------------------
    finally:
        out_df = pd.DataFrame(results)
        out_df.to_csv(config["output_file"], index=False)
        print(f"Results saved to {config['output_file']} (partial if interrupted)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CAA/ABC pattern extraction")
    parser.add_argument("config", help="Path to JSON config file")
    args = parser.parse_args()
    main(args.config)



