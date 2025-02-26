import os
import subprocess
from Bio import SeqIO
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# Step 1: Run Prokka for gene annotation
def run_prokka(input_fna, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    command = f"prokka --outdir {output_dir} --prefix genome {input_fna}"
    subprocess.run(command, shell=True, check=True)
    return os.path.join(output_dir, "genome.gff")

# Step 2: Run BProm for promoter prediction
def run_bprom(input_fna, output_gff):
    command = f"bprom -f {input_fna} -o {output_gff}"
    subprocess.run(command, shell=True, check=True)
    return output_gff

# Step 3: Run OperonMapper for operon prediction
def run_operon_mapper(input_fna, gff_file, output_gff):
    command = f"operon_mapper -i {input_fna} -g {gff_file} -o {output_gff}"
    subprocess.run(command, shell=True, check=True)
    return output_gff

# Step 4: Filter functional regions using BEDTools
def filter_regions(genome_bed, gff_file, output_bed):
    command = f"bedtools subtract -a {genome_bed} -b {gff_file} > {output_bed}"
    subprocess.run(command, shell=True, check=True)
    return output_bed

# Step 5: Extract features for neutral sites
def extract_features(neutral_sites_bed, genome_fasta):
    neutral_sites = pd.read_csv(neutral_sites_bed, sep="\t", header=None, names=["chrom", "start", "end"])
    neutral_sites["length"] = neutral_sites["end"] - neutral_sites["start"]
    
    # Example features: GC content, distance to nearest gene, etc.
    neutral_sites["gc_content"] = neutral_sites.apply(
        lambda row: calculate_gc_content(genome_fasta, row["chrom"], row["start"], row["end"]), axis=1
    )
    return neutral_sites

def calculate_gc_content(genome_fasta, chrom, start, end):
    for record in SeqIO.parse(genome_fasta, "fasta"):
        if record.id == chrom:
            sequence = str(record.seq[start:end])
            gc_count = sequence.upper().count("G") + sequence.upper().count("C")
            return gc_count / len(sequence) if len(sequence) > 0 else 0
    return 0

# Step 6: Train a machine learning model to score neutral sites
def train_model(features):
    # Example: Use GC content and length as features
    X = features[["length", "gc_content"]]
    y = [1] * len(features)  # Dummy labels (neutral sites are positive)
    
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Train a Random Forest classifier
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # Predict scores
    y_pred = model.predict_proba(X_test)[:, 1]
    print(f"Model accuracy: {accuracy_score(y_test, y_pred > 0.5)}")
    
    return model

# Step 7: Score neutral sites using the trained model
def score_neutral_sites(model, neutral_sites):
    X = neutral_sites[["length", "gc_content"]]
    neutral_sites["score"] = model.predict_proba(X)[:, 1]
    return neutral_sites

# Main pipeline
def main(input_fna):
    # Step 1: Run Prokka
    prokka_gff = run_prokka(input_fna, "prokka_output")
    
    # Step 2: Run BProm
    bprom_gff = run_bprom(input_fna, "bprom_output.gff")
    
    # Step 3: Run OperonMapper
    operon_gff = run_operon_mapper(input_fna, prokka_gff, "operon_output.gff")
    
    # Step 4: Filter functional regions
    neutral_sites_bed = filter_regions("genome.bed", prokka_gff, "non_coding_regions.bed")
    neutral_sites_bed = filter_regions(neutral_sites_bed, bprom_gff, "non_coding_no_promoters.bed")
    neutral_sites_bed = filter_regions(neutral_sites_bed, operon_gff, "neutral_sites_candidates.bed")
    
    # Step 5: Extract features
    neutral_sites = extract_features(neutral_sites_bed, input_fna)
    
    # Step 6: Train a machine learning model
    model = train_model(neutral_sites)
    
    # Step 7: Score neutral sites
    scored_sites = score_neutral_sites(model, neutral_sites)
    
    # Save final output
    scored_sites.to_csv("final_neutral_sites.csv", index=False)
    print("Pipeline completed. Neutral sites saved to final_neutral_sites.csv")

# Run the pipeline
if __name__ == "__main__":
    input_fna = "input.fna"  # Replace with your input file
    main(input_fna)