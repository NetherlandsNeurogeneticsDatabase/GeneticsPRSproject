#!/usr/bin/env python3
import argparse
import json
import re
import pandas as pd
from collections import Counter
from tqdm import tqdm
from scripts.utils import extract_features_from_text


def main(config_path):
    # Load config
    with open(config_path, "r") as f:
        config = json.load(f)

    # Load data
    df = pd.read_csv(config["input_file"], sep="\t" if config.get("tsv", True) else ",")
    feature_cols = config["feature_cols"]
    id_cols = config["id_cols"]

    # Compile patterns
    patterns = [re.compile(p, re.IGNORECASE | re.VERBOSE) for p in config["patterns"]]

    results = []

    try:
        for _, row in tqdm(df.iterrows(), total=len(df)):
            sentence = str(row.get("Sentence", "") or "")
            has_match = any(p.search(sentence) for p in patterns)

            if not has_match:
                pred_dict = {col: None for col in feature_cols}
            else:
                preds = []
                #print(f"{row['NBB_nr']} with {row['Sentence']}")
                for _ in range(config.get("n_runs", 10)):
                    p = extract_features_from_text(
                        sentence,
                        config["prompt_file"],
                        token=config.get("token", 20),
                        temp_val=config.get("temp_val", 0.1),
                        model_name=config.get("model_name", "meta-llama/Llama-3.1-8B-Instruct"),
                    ) #or {col: None for col in feature_cols}
                    preds.append(p)

                # Aggregate
                pred_dict = {}
                for col in feature_cols:
                    #print("DEBUG extract_features_from_text output type:", type(p),col, "value:")
                    values = [p.get(col) for p in preds if p.get(col) is not None]
                    if values:
                        most_common, count = Counter(values).most_common(1)[0]
                        pred_dict[col] = most_common
                        pred_dict[f"{col}_score"] = round(count / len(preds), 1)
                    else:
                        pred_dict[col] = None
                        pred_dict[f"{col}_score"] = None

            results.append({**{col: row.get(col) for col in id_cols}, **pred_dict})

    except KeyboardInterrupt:
        print("\n⚠️ Interrupted by user — saving partial results...")

    finally:
        out_df = pd.DataFrame(results)
        out_df.to_csv(config["output_file"], index=False)
        print(f"Results saved to {config['output_file']} (partial if interrupted)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CAA/ABC pattern extraction")
    parser.add_argument("config", help="Path to JSON config file")
    args = parser.parse_args()
    main(args.config)
