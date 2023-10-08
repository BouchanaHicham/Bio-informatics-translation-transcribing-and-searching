from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez

# [*] Means Stop => Check the Standard_Genetic_Code.jpg I Downloaded
Standard_Genetic_Code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Task 1: Read a genetic sequence from a FASTA file
fasta_file = "sequence.fasta"
sequence_record = SeqIO.read(fasta_file, "fasta")

print("sequence_record",sequence_record)
print(" ------------ [Sequence] ------------ ")
print(sequence_record.seq)

# Task 2: Transcribe the genetic sequence into mRNA
mRNA_sequence = sequence_record.seq.transcribe()
print(" ------------ [mRNA Sequence] ------------ ")
print(mRNA_sequence)

# Task 3: Translate mRNA into a protein sequence
# Amino Acid Sequence = Protein Sequences
protein_sequence = mRNA_sequence.translate()
print(" ------------ [Protein Sequence] ------------ ")
print(protein_sequence)

print(" ------------ [Protein Sequence Using My Func] ------------ ")

def translate_hardcoded(mRNA_sequence):
    # Transcribe and translate the mRNA sequence
    protein_sequence = ""
    for i in range(0, len(mRNA_sequence), 3):
        codon = mRNA_sequence[i: i + 3]  # We get our codon [3 Letters]
        amino_acid = Standard_Genetic_Code.get(codon, "X")  # Use "X" for unknown codons
        protein_sequence += amino_acid

    print(protein_sequence)

translate_hardcoded(mRNA_sequence)

user_input = int(input("[1]: Search for homologous sequences using BLAST\n[2]: Search Gene Id In NCBI Gene Bank\n"))
if(user_input == 1):
    # Task 4: Search for homologous sequences using BLAST,
    # BLAST: Basic Local Alignment Search Tool
    print("\nSearching for homologous sequences using BLAST...")
    result_handle = NCBIWWW.qblast("blastn", "nt", str(sequence_record.seq))
    #It uses the 'blastn' algorithm to search against the nt (nucleotide) database,
    # and it uses the DNA sequence as the query.
    blast_records = NCBIXML.parse(result_handle)
    #This allows you to iterate through the results.
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            print("Alignment Title:", alignment.title)
            print("Alignment Length:", alignment.length)
            print("Alignment E-value:", alignment.hsps[0].expect)
            # E-value (a measure of the statistical significance of the alignment)
            print("Alignment Sequence:", alignment.hsps[0].sbjct)
            print()
elif(user_input == 2):

    # Task 5: Explore the NCBI GenBank database
    Entrez.email = "xxmasterxx249@gmail.com"
    gene_id = int(input("Enter Your Gene Id:\n"))
    #gene_id = "3043"  # Here is the NCBI Gene ID for the Beta-Globin (HBB) gene in humans: Gene ID: 3043

    db_choice = int(input("[1]: Gene\n[2]: Nucleotide\n"))
    if(db_choice == 1):
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
        # Read and print the gene information
        gene_record = handle.read()

    elif(db_choice == 2):
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
        gene_record = SeqIO.read(handle, "genbank")
        handle.close()

    print("\nGenBank Record for Gene ID:", gene_id)
    print(gene_record)

# Question For Prof (Making Sure):
# Transcription: DNA => mRNA and Vice versa
# Translating: mRNA or DNA => Protein/amino acid
# Alignment Sequence is somewhere in the DNA sequence ? like a portion ?
# And why doesn't the ID search come up the same as the one I got 3043 HBB, and why is the Seq('') diff
