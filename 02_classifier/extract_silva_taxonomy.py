import gzip
from Bio import SeqIO

input_fasta = "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
output_taxonomy = "silva-138-2-taxonomy.txt"

with gzip.open(input_fasta, "rt") as fasta_in, open(output_taxonomy, "w") as out:
    for record in SeqIO.parse(fasta_in, "fasta"):
        seq_id = record.id
        try:
            taxonomy = record.description.split(None, 1)[1]
        except IndexError:
            taxonomy = "Unclassified"
        out.write(f"{seq_id}\t{taxonomy}\n")
