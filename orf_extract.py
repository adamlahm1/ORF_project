#Set path to ORF folder itinerary through cd 'path' and create docker image
'''
first run the first part of the code in git bash:
> sed -n '31,113p' orf_extract.py > part1.py
> python part1.py
then create a Dockerfile in file path dockerfiles/blast/v2.16.0 within the file :
FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y wget libgomp1

RUN pip install --upgrade pip && \
    pip install matplotlib pandas biopython

WORKDIR /app

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    mv ncbi-blast-2.16.0+ /usr/local/ && \
    rm ncbi-blast-2.16.0+-x64-linux.tar.gz

ENV PATH="/usr/local/ncbi-blast-2.16.0+/bin:$PATH"

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/db/v5/swissprot.tar.gz && \
    mkdir -p /db/swissprot && mv swissprot* /db/swissprot && \
    tar -xzvf /db/swissprot/swissprot.tar.gz -C /db/swissprot

COPY . /app

'''
### 1st part
from Bio import SeqIO
import matplotlib.pyplot as plt


def read_fasta_sequences(filepath):
    """Reads a FASTA file and returns a list of (header, sequence) tuples."""
    try:
        return [(record.id, str(record.seq)) for record in SeqIO.parse(filepath, "fasta")]
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        return []


def extract_orfs_from_sequence(sequence, seq_id):
    """Extract non-overlapping ORFs from a single sequence."""
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    i = 0
    while i <= len(sequence) - 3:
        if sequence[i:i+3] == start_codon:
            start = i
            i += 3
            while i <= len(sequence) - 3:
                codon = sequence[i:i+3]
                if codon in stop_codons:
                    end = i + 3
                    reading_frame = start % 3
                    orfs.append({
                        "id": f"{seq_id}_ORF{len(orfs)+1}",
                        "seq": sequence[start:end],
                        "length": end - start,
                        "source": seq_id,
                        "start": start,
                        "end": end,
                        "frame": reading_frame
                    })
                    break
                i += 3
        else:
            i += 1
    return orfs


def write_gff(orfs, gff_output):
    """Write ORF data to a GFF file."""
    with open(gff_output, "w") as f:
        for orf in orfs:
            f.write(
                f"{orf['source']}\tORF_Pipeline\tCDS\t{orf['start']+1}\t{orf['end']}\t.\t+\t.\t"
                f"ID={orf['id']};length={orf['length']};frame={orf['frame']}\n"
            )


def write_fasta(orfs, fasta_output):
    """Write ORF sequences to a multi-FASTA file."""
    with open(fasta_output, "w") as f:
        for orf in orfs:
            f.write(f">{orf['id']}\n{orf['seq']}\n")


def run_orf_pipeline(fasta_file, gff_file, fasta_out_file):
    """Main function to coordinate ORF extraction and output."""
    sequences = read_fasta_sequences(fasta_file)
    all_orfs = []

    for index, (seq_id, sequence) in enumerate(sequences):
        orfs = extract_orfs_from_sequence(sequence, f"Seq{index + 1}")
        all_orfs.extend(orfs)

    write_gff(all_orfs, gff_file)
    write_fasta(all_orfs, fasta_out_file)
    print(f"Extracted {len(all_orfs)} ORFs from {len(sequences)} sequences.")


# to Run 1st part
run_orf_pipeline(
    fasta_file="Homo_sapiens_cdna_assembled.fasta",
    gff_file="orfs_output.gff",
    fasta_out_file="orfs_output.fasta"
)

#Bash
'''
Open docker desktop on computer and in git bash after setting path to the project file (ORF): docker build -t blast:v2.16.0 .      #to build the image
After image has been mounted enter the docker environment and mount volume with : docker run -it --rm -v "C:\Users\Adam\OneDrive\Documents\M2 DM Bio\S2\Data mining in biosci\ORF":/ORF blast:v2.16.0 bash
Once in the environment set path : cd /ORF , then export BLASTDB=/db/swissprot/swissprot so that blast can look for taxid in the right place
then run :blastx -db /db/swissprot/swissprot -query /ORF/orfs_output.fasta \
       -taxids 9606 \
       -outfmt 7 \
       -out /ORF/blast_results.tsv \
       -num_threads 4

finally run the 2nd part of the code through:
> sed -n '130,232p' orf_extract.py > part2.py
> python part2.py
'''

### 2nd part
import pandas as pd
import matplotlib.pyplot as plt

def parse_blast(blast_file, evalue_thresh=1e-5):
    valid, invalid = [], []
    seen = set()

    with open(blast_file) as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            qid = parts[0]
            evalue = float(parts[10])

            if qid not in seen:
                seen.add(qid)
                if evalue < evalue_thresh:
                    valid.append((qid, evalue))
                else:
                    invalid.append((qid, evalue))
    return valid, invalid


def append_gff(gff_path, valid_hits):
    # Read GFF safely
    df = pd.read_csv(gff_path, sep='\t', header=None, comment='#')

    if df.shape[1] != 9:
        raise ValueError(f"GFF file should have 9 columns, found {df.shape[1]}. Please check the format.")

    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    # Add CDS annotations
    new_rows = []
    for qid, evalue in valid_hits:
        match = df[df['attributes'].str.contains(qid)]
        if not match.empty:
            row = match.iloc[0]
            new_row = row.copy()
            new_row['type'] = 'CDS'
            new_row['attributes'] += f";validated_by=BLAST;evalue={evalue}"
            new_rows.append(new_row)

    df_all = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    df_all.to_csv("annotated_orfs.gff", sep='\t', header=False, index=False)


def compute_fpr(valid, invalid):
    total = len(valid) + len(invalid)
    if total == 0:
        return 0.0
    return len(invalid) / total


def plot_orf_distributions(gff_path, valid_ids, invalid_ids):
    df = pd.read_csv(gff_path, sep='\t', header=None, comment='#')
    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    df_orfs = df[df['type'] == 'CDS'].copy()
    df_orfs['length'] = df_orfs['end'] - df_orfs['start'] + 1
    df_orfs['id'] = df_orfs['attributes'].str.extract(r'ID=([^;]+)')

    df_orfs['label'] = df_orfs['id'].apply(
        lambda x: 'valid' if x in set(valid_ids) else 'invalid' if x in set(invalid_ids) else 'unknown'
    )

    plt.figure()
    plt.hist(df_orfs['length'], bins=40, edgecolor='black')
    plt.title("ORF Length Distribution (All)")
    plt.xlabel("Length (bp)")
    plt.ylabel("Frequency")
    plt.savefig("orf_lengths_all.png")

    plt.figure()
    for label in ['valid', 'invalid']:
        subset = df_orfs[df_orfs['label'] == label]
        plt.hist(subset['length'], bins=30, alpha=0.6, label=label)
    plt.title("ORF Length Distribution: Valid vs Invalid")
    plt.xlabel("Length (bp)")
    plt.ylabel("Count")
    plt.legend()
    plt.savefig("orf_lengths_valid_vs_invalid.png")
    plt.show()

#to Run 2nd part

blast_file = "blast_results.tsv"
gff_file = "orfs_output.gff"

valid_hits, invalid_hits = parse_blast(blast_file)
append_gff(gff_file, valid_hits)

fpr = compute_fpr(valid_hits, invalid_hits)
print(f"False Positive Rate: {fpr:.2f}")

plot_orf_distributions(gff_file,
                       valid_ids=[x[0] for x in valid_hits],
                       invalid_ids=[x[0] for x in invalid_hits])
