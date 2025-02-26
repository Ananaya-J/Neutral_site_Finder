import os
import gzip
import shutil
import re
import sqlite3
import pandas as pd
from urllib.request import urlretrieve
from Bio import SeqIO

def download_and_extract_fna(url, output_fna):
    """Download and extract the .fna.gz file."""
    compressed_file = output_fna + ".gz"
    print("Downloading genome...")
    urlretrieve(url, compressed_file)
    print("Extracting file...")
    with gzip.open(compressed_file, 'rb') as f_in, open(output_fna, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"File saved as {output_fna}")

def parse_fasta(fasta_file):
    """Load and return the genome sequence from a FASTA file."""
    with open(fasta_file, "r") as file:
        record = next(SeqIO.parse(file, "fasta"))
        return str(record.seq)

def find_promoters(sequence, patterns, flanking_length=15):
    """Identify promoters and extract sequences with flanking regions."""
    promoters = []
    for pattern in patterns:
        for match in re.finditer(pattern, sequence, re.IGNORECASE):
            start = max(0, match.start() - flanking_length)
            end = min(len(sequence), match.end() + flanking_length)
            extracted_seq = sequence[start:end]
            promoters.append((match.start(), match.group(), extracted_seq))
    return promoters

def save_to_csv(promoters, output_csv):
    """Save promoter sequences and locations to a CSV file."""
    df = pd.DataFrame(promoters, columns=["Position", "Promoter", "Sequence_with_flanks"])
    df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

def create_database(csv_file, db_file):
    """Create a database from the CSV file."""
    conn = sqlite3.connect(db_file)
    df = pd.read_csv(csv_file)
    df.to_sql("promoters", conn, if_exists="replace", index=False)
    conn.close()
    print(f"Database created: {db_file}")

if __name__ == "__main__":
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    fasta_file = "genome.fna"
    csv_file = "promoters.csv"
    db_file = "promoters.db"
    promoter_patterns = [r'TTGACA.{15,19}TATAAT', r'TATAAT', r'TTGACA']
    
    download_and_extract_fna(url, fasta_file)
    sequence = parse_fasta(fasta_file)
    promoters = find_promoters(sequence, promoter_patterns)
    save_to_csv(promoters, csv_file)
    create_database(csv_file, db_file)