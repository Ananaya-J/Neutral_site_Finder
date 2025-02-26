import pandas as pd

# Load the CSV file
file_path = "/GCF_000005845.2_ASM584v2_genomic.csv"

# Read and clean column names
df = pd.read_csv(file_path)
df.columns = df.columns.str.strip().str.lower()  # Normalize column names

# Print available columns to check correctness
print("Available Columns:", df.columns)

# Ensure required columns exist
required_columns = {'start', 'end', 'strand'}
if not required_columns.issubset(df.columns):
    raise KeyError(f"Missing required columns! Found: {df.columns}")

# Check if a functional annotation column exists
function_col = None
for col in df.columns:
    if 'function' in col or 'product' in col or 'gene' in col:
        function_col = col
        break

# Convert start & end positions to integers
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)

# Remove duplicate genes that are already present in CDS
df = df.drop_duplicates(subset=['start', 'end'])

# Sort by start position
df = df.sort_values(by="start").reset_index(drop=True)

# Function to group ORFs into operons
def identify_operons(df, max_intergenic_distance=50, min_genes=3):
    """Groups genes into operons while ensuring at least 3 genes per operon."""
    operons = []
    current_operon = [df.iloc[0].to_dict()]  # Convert to dictionary

    for i in range(1, len(df)):
        prev_end = current_operon[-1]["end"]
        gene_data = df.iloc[i].to_dict()  # Convert each row to dictionary

        # Check if genes are close enough and on the same strand
        if gene_data["start"] - prev_end <= max_intergenic_distance and gene_data["strand"] == current_operon[-1]["strand"]:
            current_operon.append(gene_data)
        else:
            # Ensure only valid operons (≥3 genes) are considered
            if len(current_operon) >= min_genes:
                operons.append(current_operon)
            current_operon = [gene_data]

    # Add the last operon if it meets the criteria
    if len(current_operon) >= min_genes:
        operons.append(current_operon)

    return operons

# Identify operons
operons = identify_operons(df)

# Print all identified operons
for i, operon in enumerate(operons):
    print(f"\n🔹 Operon {i+1} (contains {len(operon)} genes):")
    for gene in operon:
        if function_col:
            print(f"  - Gene {gene['start']} to {gene['end']} ({gene['strand']}): {gene[function_col]}")
        else:
            print(f"  - Gene {gene['start']} to {gene['end']} ({gene['strand']})")

# Print total number of identified operons
print(f"\n✅ Total identified operons: {len(operons)}")