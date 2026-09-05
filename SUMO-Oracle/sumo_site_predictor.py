#!/usr/bin/env python3
import argparse
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
import joblib
import json
import subprocess
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def parse_structure(structure_file,
                    fasta_file= None,
                    multifragment_file= None):

    stem = Path(structure_file).stem
    protein_name = stem.split("-")[1]
    filename = Path(structure_file).name

    if fasta_file and multifragment_file:
        with open(multifragment_file, 'r') as f:
            multifrag_files = {line.strip() for line in f}

        if filename in multifrag_files:
            print(f"{filename} is a multifragment protein, reading full sequence from {fasta_file}")

            with open(fasta_file, 'r') as f:
                header, seq = None, []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if header and protein_name in header:
                            return {protein_name: ''.join(seq)}
                        header = line[1:]
                        seq = []
                    else:
                        seq.append(line)
                if header and protein_name in header:
                    return {protein_name: ''.join(seq)}

            print(f"{protein_name} is not found in {fasta_file}")

    ext = Path(structure_file).suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    structure = parser.get_structure("s", structure_file)
    builder = PPBuilder()
    sequences = {}
    for model in structure:
        for chain in model:
            seq = "".join(str(pp.get_sequence()) for pp in builder.build_peptides(chain))
            if seq:
                sequences[protein_name] = seq
        break
    return sequences


def get_fragment_files(protein_name, input_dir,
                       multifrag_files):

    input_dir = Path(input_dir)
    ext = None

    for e in ['.cif', '.pdb']:
        if (input_dir / f"AF-{protein_name}-F1-model_v6{e}").exists():
            ext = e
            break

    if ext is None:
        return [], []

    expected = sorted(
        f for f in multifrag_files
        if f"AF-{protein_name}-" in f and f.endswith(ext)
    )

    found   = [f for f in expected if (input_dir / f).exists()]
    missing = [f for f in expected if not (input_dir / f).exists()]

    return found, missing


def find_lysine_positions(sequence):
    positions = []
    for i, residue in enumerate(sequence.upper()):
        if residue == 'K':
            positions.append(i + 1)
    return positions


def extract_kmer_window(sequence, position,
                        kmer_size= 5):
    k_idx = position - 1
    half_win = kmer_size // 2
    left  = max(0, k_idx - half_win)
    right = k_idx + half_win + 1
    window = sequence[left:right]
    if k_idx - half_win < 0:
        window = 'X' * (half_win - k_idx) + window
    if right > len(sequence):
        window = window + 'X' * (right - len(sequence))
    return window


def load_aaindex(aaindex_file):
    df = pd.read_csv(aaindex_file)
    df = df.set_index('AccNo')
    return df


def load_features(metadata_file):
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    return metadata["feature_columns"]


def extract_aaindex_features(kmer, aaindex_df,
                              features, kmer_size= 5):

    feature_values = {}
    half = kmer_size // 2
    kmer_pos = {f'pos{r}': half + r for r in range(-half, half + 1)}

    for feature in features:
        if feature.startswith('sum_'):
            index_name = feature[4:]
            vals = []
            has_x = False
            for aa in kmer:
                if aa == 'X':
                    has_x = True
                elif index_name in aaindex_df.index and aa in aaindex_df.columns:
                    vals.append(float(aaindex_df.loc[index_name, aa]))
            feature_values[feature] = sum(vals) if not has_x else None

        elif '_' in feature:
            parts = feature.split('_', 1)
            pos_prefix, index_name = parts[0], parts[1]
            if pos_prefix in kmer_pos:
                idx = kmer_pos[pos_prefix]
                aa = kmer[idx] if idx < len(kmer) else None
                if aa and aa != 'X' and index_name in aaindex_df.index and aa in aaindex_df.columns:
                    feature_values[feature] = float(aaindex_df.loc[index_name, aa])
                else:
                    feature_values[feature] = None

    return feature_values


def write_fasta(sequences, fasta_path):
    with open(fasta_path, 'w') as f:
        for protein_name, seq in sequences.items():
            f.write(f">{protein_name}\n{seq}\n")


def run_metapredict(fasta_file):
    cmd = ["metapredict-predict-disorder", fasta_file]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"metapredict failed:\n{result.stderr}")
    else:
        print(f"metapredict completed successfully")
    return result.returncode


def parse_metapredict_output(disorder_csv):
    data = []
    with open(disorder_csv, 'r') as f:
        for line in f:
            row = line.strip().split(",")
            data.append(row)

    max_columns = max(len(row) for row in data)
    for row in data:
        while len(row) < max_columns:
            row.append(None)

    df = pd.DataFrame(data)
    protein_col = df.iloc[:, 0].astype(str)
    score_cols = df.iloc[:, 2:]

    rows = []
    for protein, scores in zip(protein_col, score_cols.itertuples(index=False)):
        for pos, disorder in enumerate(scores, start=1):
            if pd.notnull(disorder):
                rows.append({"protein": protein, "pos": pos, "disorder": float(disorder)})

    return pd.DataFrame(rows)


def extract_all_helices_from_dssp_dir(dssp_dir, protein_name):

    all_helices = []
    dssp_dir = Path(dssp_dir)

    for dssp_cif in sorted(dssp_dir.glob(f"AF-{protein_name}-*_dssp.cif")):
        cif = MMCIF2Dict(str(dssp_cif))

        if "_struct_conf.conf_type_id" not in cif:
            continue

        def as_list(v):
            return [v] if isinstance(v, str) else v

        conf_type = as_list(cif["_struct_conf.conf_type_id"])
        beg_chain = as_list(cif["_struct_conf.beg_auth_asym_id"])
        beg_seqid = as_list(cif["_struct_conf.beg_auth_seq_id"])
        end_seqid = as_list(cif["_struct_conf.end_auth_seq_id"])
        conf_id   = as_list(cif["_struct_conf.id"])

        for i, t in enumerate(conf_type):
            if t.strip() == "HELX_RH_AL_P":
                all_helices.append({
                    "Protein":  protein_name,
                    "Chain":    beg_chain[i],
                    "Helix_ID": conf_id[i],
                    "Start":    int(beg_seqid[i]),
                    "End":      int(end_seqid[i])
                })

    result = pd.DataFrame(all_helices)
    if not result.empty:
        result = deduplicate_helices(result)
    return result


def deduplicate_helices(helix_df):

    if helix_df.empty:
        return helix_df

    helix_df = helix_df.copy()
    helix_df["length"] = helix_df["End"] - helix_df["Start"]
    helix_df = helix_df.sort_values("length", ascending=False).reset_index(drop=True)

    kept = []
    for _, row in helix_df.iterrows():
        overlaps = any(row["Start"] <= k["End"] and k["Start"] <= row["End"] for k in kept)
        if not overlaps:
            kept.append(row)

    return pd.DataFrame(kept).drop(columns="length").reset_index(drop=True)


def compute_dist_nearest(df, helix_df):

    if helix_df.empty or "Protein" not in helix_df.columns:
        print("No helix data, Dist_nearest will be None")
        df["dist_nearest_helix"] = None
        return df

    merged = df.merge(helix_df, on="Protein", how="left")
    merged["Dist_N"] = merged["K-pos"] - merged["Start"]
    merged["Dist_C"] = merged["K-pos"] - merged["End"]
    merged["In_helix"] = ((merged["Dist_N"] * merged["Dist_C"]) <= 0).map({True: "Yes", False: "No"})
    merged["abs_Dist_N"] = merged["Dist_N"].abs()
    merged["abs_Dist_C"] = merged["Dist_C"].abs()
    merged["dist_nearest_helix"] = merged[["abs_Dist_N", "abs_Dist_C"]].min(axis=1)

    in_ss     = merged[merged["In_helix"] == "Yes"]
    not_in_ss = merged[merged["In_helix"] == "No"]

    nearest_not_in = (
        not_in_ss.sort_values("dist_nearest_helix")
                 .groupby(["Protein", "K-pos"], as_index=False)
                 .first()
    )

    nearest = pd.concat([in_ss, nearest_not_in], ignore_index=True)
    nearest = nearest.sort_values("In_helix", ascending=False)
    nearest = nearest.drop_duplicates(subset=["Protein", "K-pos"], keep="first")

    return df.merge(nearest[["Protein", "K-pos", "dist_nearest_helix"]],
                    on=["Protein", "K-pos"], how="left")


def load_asa(asa_dir: str, protein_name: str, input_dir: str,
             multifrag_files: set, filename: str) -> pd.DataFrame:

    asa_dir = Path(asa_dir)
    input_dir = Path(input_dir)

    if filename in multifrag_files:
        cif_frags = sorted(input_dir.glob(f"AF-{protein_name}-*-model_v6.cif"))
        csv_files = []
        for cif in cif_frags:
            frag_id = cif.stem.replace("model_v6", "model_v4")
            csv_path = asa_dir / f"{frag_id}.csv"
            if csv_path.exists():
                csv_files.append(csv_path)
            else:
                print(f"ASA file not found: {csv_path.name}")
    else:
        stem = f"AF-{protein_name}-F1-model_v4"
        csv_path = asa_dir / f"{stem}.csv"
        csv_files = [csv_path] if csv_path.exists() else []

    if not csv_files:
        return pd.DataFrame()

    frames = []
    for f in csv_files:
        df = pd.read_csv(f)
        df = df[df["residue_name"].astype(str).str.strip().str.upper() == "LYS"].copy()
        df["accession"] = df["structure"].astype(str).str.split("-").str[1]
        df = df[df["accession"] == protein_name]
        if not df.empty:
            frames.append(df[["accession", "residue_number", "asa_abs"]])

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.drop_duplicates(subset=["accession", "residue_number"], keep="first")
    combined = combined.rename(columns={"residue_number": "K-pos",
                                        "accession": "Protein",
                                        "asa_abs": "absolute_asa"})
    return combined[["Protein", "K-pos", "absolute_asa"]]


def predict(df, model_file, metadata_file):

    model = joblib.load(model_file)

    with open(metadata_file, 'r') as f:
        metadata = json.load(f)

    feature_cols = metadata["feature_columns"]
    threshold_f1  = metadata["optimal_threshold_f1"]
    threshold_mcc = metadata["optimal_threshold_mcc"]

    X = df[feature_cols].to_numpy()
    df["sumoylation_probability"] = model.predict_proba(X)[:, 1]
    df["sumo_f1_prediction"]  = (df["sumoylation_probability"] >= threshold_f1).astype(int)
    df["sumo_mcc_prediction"] = (df["sumoylation_probability"] >= threshold_mcc).astype(int)

    print(f"Predictions complete: "
          f"{df['sumo_f1_prediction'].sum()} sites predicted (F1 threshold={threshold_f1}), "
          f"{df['sumo_mcc_prediction'].sum()} sites predicted (MCC threshold={threshold_mcc})")
    return df


def struct_to_kmer_dataframe(
    struct_file,
    input_dir,
    multifrag_files,
    aaindex_file= None,
    features_file= None,
    kmer_size= 5,
    fasta_file= None,
    multifragment_file= None,
    all_helices_file= None,
    dssp_dir= None,
    asa_dir= None,
    train_medians_file= None,
    model_file= None,
) -> pd.DataFrame:

    sequences = parse_structure(struct_file,
                                fasta_file=fasta_file,
                                multifragment_file=multifragment_file)

    if not sequences:
        print("No sequences found in structure file.")
        return pd.DataFrame()

    protein_name = list(sequences.keys())[0]
    filename = Path(struct_file).name

    if filename in multifrag_files:
        found, missing = get_fragment_files(protein_name, input_dir, multifrag_files)
        if missing:
            print(f"Missing fragment CIF files needed for DSSP:")
            for m in missing:
                print(f"       - {m}")
            print(f"DSSP will only use {len(found)} available fragments")
        else:
            print(f"All {len(found)} fragment CIF files present")

    fasta_extract_dir = Path(struct_file).parent.parent / "fasta_extract"
    fasta_extract_dir.mkdir(exist_ok=True)
    fasta_name = str(fasta_extract_dir / (Path(struct_file).stem + ".fasta"))
    write_fasta(sequences, fasta_name)
    run_metapredict(fasta_name)

    disorder_csv_path = Path("disorder_scores.csv")
    disorder_df = parse_metapredict_output(str(disorder_csv_path))
    disorder_csv_path.unlink(missing_ok=True)

    proteins, positions, kmers = [], [], []
    for protein_header, sequence in sequences.items():
        for pos in find_lysine_positions(sequence):
            kmer = extract_kmer_window(sequence, pos, kmer_size)
            proteins.append(protein_header)
            positions.append(pos)
            kmers.append(kmer)

    df = pd.DataFrame({'Protein': proteins, 'K-pos': positions, '5-mer': kmers})

    if aaindex_file and features_file:
        aaindex_df = load_aaindex(aaindex_file)
        features = load_features(features_file)

        feature_data = [extract_aaindex_features(kmer, aaindex_df, features, kmer_size)
                        for kmer in df['5-mer']]
        features_df = pd.DataFrame(feature_data)
        df = pd.concat([df, features_df], axis=1)

        df = df.merge(
            disorder_df.rename(columns={"protein": "Protein", "pos": "K-pos"}),
            on=["Protein", "K-pos"], how="left"
        )

        if all_helices_file and Path(all_helices_file).exists():
            helix_df = pd.read_csv(all_helices_file)
            helix_df = helix_df[helix_df['Protein'] == protein_name].copy()
            helix_df = deduplicate_helices(helix_df)
            df = compute_dist_nearest(df, helix_df)
        elif dssp_dir and Path(dssp_dir).exists():
            helix_df = extract_all_helices_from_dssp_dir(dssp_dir, protein_name)
            if not helix_df.empty:
                print(f"Extracted {len(helix_df)} helices from DSSP files")
                df = compute_dist_nearest(df, helix_df)
            else:
                print("No helices found in DSSP files")
                df["dist_nearest_helix"] = None
        else:
            print("No helix data provided")
            df["dist_nearest_helix"] = None

        if asa_dir and Path(asa_dir).exists():
            asa_df = load_asa(asa_dir, protein_name, input_dir, multifrag_files, filename)
            if not asa_df.empty:
                df = df.merge(asa_df, on=["Protein", "K-pos"], how="left")
                print(f"absolute_asa for {protein_name} ({asa_df['K-pos'].nunique()} lysine sites)")
            else:
                print(f"No ASA data found for {protein_name}")
                df["absolute_asa"] = None
        else:
            df["absolute_asa"] = None

        if train_medians_file and Path(train_medians_file).exists():
            train_medians = joblib.load(train_medians_file)
            for col, median in train_medians.items():
                if col in df.columns:
                    df[col] = df[col].fillna(median)
            print(f"Imputed missing values using train medians")
        else:
            print(f"train_medians.joblib not found")

        if model_file and Path(model_file).exists() and features_file:
            df = predict(df, model_file, features_file)
        else:
            print(f"Model file not found")

    return df


def process_folder(input_dir, args):

    input_dir = Path(input_dir)

    multifrag_files = set()
    if args.multifragment and Path(args.multifragment).exists():
        with open(args.multifragment, 'r') as f:
            multifrag_files = {line.strip() for line in f}

    cif_f1 = {f.stem.split("-")[1]: f for f in input_dir.glob("AF-*-F1-model_v6.cif")}
    pdb_f1 = {f.stem.split("-")[1]: f for f in input_dir.glob("AF-*-F1-model_v6.pdb")}
    merged = {**pdb_f1, **cif_f1}
    f1_files = sorted(merged.values(), key=lambda f: f.name)

    if not f1_files:
        print(f"No F1 structure files found in {input_dir}")
        return

    # print(f"Found {len(f1_files)} F1 file(s) to process")
    all_dfs = []

    for f1 in f1_files:
        filename = f1.name
        protein_name = filename.split("-")[1]
        print(f"\n{'='*60}")
        print(f"Processing: {protein_name}")

        if filename in multifrag_files:
            found, missing = get_fragment_files(protein_name, str(input_dir), multifrag_files)
            if missing:
                print(f"Missing fragment CIF files needed for DSSP:")
                for m in missing:
                    print(f"       - {m}")
                print(f"DSSP will only use {len(found)} available fragments")
            else:
                print(f"All {len(found)} fragment CIF files present")

        df = struct_to_kmer_dataframe(
            str(f1),
            input_dir=str(input_dir),
            multifrag_files=multifrag_files,
            aaindex_file=args.aaindex,
            features_file=args.metadata,
            kmer_size=5,
            fasta_file=args.fasta,
            multifragment_file=args.multifragment,
            all_helices_file=args.all_helices,
            dssp_dir=args.dssp_dir,
            asa_dir=args.asa_dir,
            train_medians_file=args.train_medians,
            model_file=args.model,
        )
        if not df.empty:
            all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(args.output, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-fa', '--fasta',
                        default='multifrag_prot_idmapping_2026_06_19.fasta')
    parser.add_argument('-mf', '--multifragment',
                        default='human_AF_multifragment.txt')
    parser.add_argument('-a', '--aaindex',
                        default='aaindex.csv')
    parser.add_argument('-m', '--metadata',
                        default='model_metadata.json')
    parser.add_argument('-ah', '--all_helices')
    parser.add_argument('-dd', '--dssp_dir',
                        default='dssp')
    parser.add_argument('-as', '--asa_dir',
                        default='asa')
    parser.add_argument('-tm', '--train_medians',
                        default='train_medians.joblib')
    parser.add_argument('-rf', '--model',
                        default='random_forest.joblib')

    args = parser.parse_args()

    process_folder(args.input_dir, args)


if __name__ == '__main__':
    main()