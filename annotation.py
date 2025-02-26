import os
import gzip
import subprocess
import requests
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt

def check_requirements():
    """Ensure required Python packages are installed."""
    required_packages = ["requests", "pandas", "matplotlib", "biopython"]
    missing_packages = [pkg for pkg in required_packages if not is_package_installed(pkg)]
    
    if missing_packages:
        print(f"Installing missing packages: {', '.join(missing_packages)}")
        subprocess.run(["pip", "install"] + missing_packages, check=True)

def is_package_installed(package):
    """Check if a package is installed."""
    try:
        __import__(package)
        return True
    except ImportError:
        return False

def download_genome():
    """Download E. coli K-12 MG1655 genome from NCBI."""
    accession = "NC_000913.3"
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"
    genome_url = f"{base_url}/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gff_url = f"{base_url}/GCF_000005845.2_ASM584v2_genomic.gff.gz"

    os.makedirs("genome_data", exist_ok=True)

    genome_file = "genome_data/ecoli_genome.fna.gz"
    gff_file = "genome_data/ecoli_genome.gff.gz"

    try:
        if not os.path.exists(genome_file):
            print(f"Downloading genome file: {genome_url}")
            with open(genome_file, 'wb') as f:
                f.write(requests.get(genome_url).content)

        if not os.path.exists(gff_file):
            print(f"Downloading GFF annotation: {gff_url}")
            with open(gff_file, 'wb') as f:
                f.write(requests.get(gff_url).content)

        # Decompress files
        with gzip.open(genome_file, 'rb') as f_in, open("genome_data/ecoli_genome.fna", 'wb') as f_out:
            f_out.write(f_in.read())

        with gzip.open(gff_file, 'rb') as f_in, open("genome_data/ecoli_genome.gff", 'wb') as f_out:
            f_out.write(f_in.read())

        print("Genome and annotation downloaded and extracted successfully.")
        return "genome_data/ecoli_genome.fna", "genome_data/ecoli_genome.gff"
    
    except Exception as e:
        print(f"Error downloading genome: {e}")
        return None, None

def run_prokka_annotation(genome_file, output_dir="genome_data/prokka_results"):
    """Run Prokka annotation."""
    os.makedirs(output_dir, exist_ok=True)
    prefix = "ecoli"

    # Check if Prokka is installed
    if subprocess.run(["which", "prokka"], capture_output=True).returncode != 0:
        print("Prokka is not installed. Install it using: conda install -c bioconda prokka")
        return None

    try:
        print(f"Running Prokka annotation on {genome_file}...")
        prokka_cmd = [
            "prokka", "--outdir", output_dir, "--prefix", prefix, "--kingdom", "Bacteria",
            "--genus", "Escherichia", "--species", "coli", "--strain", "K-12", "--locustag", "ECOLI",
            "--compliant", "--cpus", "4", "--force", genome_file  # Added --force to overwrite existing results
        ]
        subprocess.run(prokka_cmd, check=True)

        prokka_gff = os.path.join(output_dir, f"{prefix}.gff")
        prokka_csv = os.path.join(output_dir, f"{prefix}_genes.csv")

        print(f"Prokka annotation completed. Results in {output_dir}")
        convert_gff_to_csv(prokka_gff, prokka_csv)
        return prokka_gff, prokka_csv

    except subprocess.CalledProcessError as e:
        print(f"Error running Prokka: {e}")
        return None, None  # Ensure function returns a tuple, preventing unpacking errors

def convert_gff_to_csv(gff_file, output_file="genome_data/ecoli_genes.csv"):
    """Convert GFF annotation to CSV format."""
    genes = []
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'gene':
                    continue

                gene_info = {
                    'chromosome': parts[0], 'source': parts[1], 'type': parts[2],
                    'start': int(parts[3]), 'end': int(parts[4]), 'strand': parts[6]
                }
                attrs = {kv.split('=')[0]: kv.split('=')[1] for kv in parts[8].split(';') if '=' in kv}
                gene_info.update(attrs)

                genes.append(gene_info)

        df = pd.DataFrame(genes)
        df.to_csv(output_file, index=False)
        print(f"Converted GFF to CSV: {output_file}")
    
    except Exception as e:
        print(f"Error converting GFF to CSV: {e}")

def compare_annotations(ncbi_csv, prokka_csv, output_file="genome_data/annotation_comparison.csv"):
    """Compare NCBI and Prokka gene annotations."""
    try:
        ncbi_genes = pd.read_csv(ncbi_csv)
        prokka_genes = pd.read_csv(prokka_csv)

        overlaps = prokka_genes.merge(ncbi_genes, on=['start', 'end', 'strand'], suffixes=('_prokka', '_ncbi'))
        overlaps.to_csv(output_file, index=False)

        print(f"Annotation comparison saved to {output_file}")
    
    except Exception as e:
        print(f"Error comparing annotations: {e}")

def plot_comparison(ncbi_csv, prokka_csv):
    """Plot comparison of gene counts."""
    try:
        ncbi_count = len(pd.read_csv(ncbi_csv))
        prokka_count = len(pd.read_csv(prokka_csv))

        plt.figure(figsize=(6, 4))
        plt.bar(['NCBI', 'Prokka'], [ncbi_count, prokka_count], color=['blue', 'red'])
        plt.ylabel('Number of genes')
        plt.title('Comparison of NCBI and Prokka annotations')
        plt.savefig("genome_data/annotation_comparison.png")
        print("Comparison plot saved as genome_data/annotation_comparison.png")

    except Exception as e:
        print(f"Error generating comparison plot: {e}")

def main():
    """Main function to run genome annotation pipeline."""
    check_requirements()
    
    genome_file, gff_file = download_genome()
    if genome_file and gff_file:
        prokka_gff, prokka_csv = run_prokka_annotation(genome_file)
        
        if prokka_gff and prokka_csv:
            ncbi_csv = convert_gff_to_csv(gff_file)
            if ncbi_csv:
                compare_annotations(ncbi_csv, prokka_csv)
                plot_comparison(ncbi_csv, prokka_csv)
    
    print("Pipeline execution completed.")

if __name__ == "__main__":
    main()

