from Bio import SeqIO

input_fasta = "SILVA_138.2_SSURef_NR99_tax_silva.fasta"
output_fasta = "SILVA_138.2_SSURef_NR99_tax_silva_DNA.fasta"

with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(input_fasta, "fasta"):
        record.seq = record.seq.transcribe().back_transcribe()  # U → T
        SeqIO.write(record, out_f, "fasta")
