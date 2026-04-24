from Bio import SeqIO

proteins=[]
output_file=open("results/analysis.txt","w")
print("Reading FASTA file...\n")
print("\n")
output_file.write("Reading FASTA file...\n")
output_file.write("\n")
for record in SeqIO.parse("data/multi_seq.fasta","fasta"):

    print("Sequence ID:", record.id)
    output_file.write("Sequence ID:"+record.id+"\n")

    seq = record.seq

    print("DNA Sequence:", seq)
    output_file.write("DNA Sequence:"+str(seq)+"\n")

    print("Length:", len(seq))
    output_file.write("Length:"+str(len(seq))+"\n")

     # GC Content
    g = seq.count("G")
    c = seq.count("C")

    gc_content = (g + c) / len(seq) * 100

    print("GC Content:", round(gc_content, 2), "%")
    output_file.write("GC Content:"+ str(round(gc_content, 2))+ "%"+"\n")

    # Calculate AT Content
    a = seq.count("A")
    t = seq.count("T")

    at_content = (a + t) / len(seq) * 100

    print("AT Content:", round(at_content, 2), "%")
    output_file.write("AT Content:"+ str(round(at_content, 2))+ "%"+"\n")

    # Base Counts
    print("A count:", a)
    output_file.write("A count:"+str(a)+"\n")
    print("T count:", t)
    output_file.write("A count:"+str(t)+"\n")
    print("G count:", g)
    output_file.write("A count:"+str(g)+"\n")
    print("C count:", c)
    output_file.write("A count:"+str(c)+"\n")

    # DNA to RNA Conversion
    rna = seq.transcribe()

    print("RNA Sequence:", rna)
    output_file.write("RNA Sequence:"+ str(rna)+"\n")

    # RNA to Protein Translation
    protein = seq.translate()

    print("Protein Sequence:", protein)
    output_file.write("Protein Sequence:"+ str(protein)+"\n")

    # store protein 
    proteins.append((record.id, protein))


    print("********************************")
    print("********************************")
    output_file.write("********************************"+"\n")
    output_file.write("********************************"+"\n")

# Find longest protein sequence
longest_id = ""
longest_protein = ""

for pid, prot in proteins:
    if len(prot) > len(longest_protein):
        longest_protein = prot
        longest_id = pid

print("\n===== LONGEST PROTEIN =====")
output_file.write("\n===== LONGEST PROTEIN ====="+"\n")
print("Sequence ID:", longest_id)
output_file.write("Sequence ID:"+ longest_id+"\n")
print("Length:", len(longest_protein))
output_file.write("Length:"+str(len(longest_protein))+"\n")
print("Protein:", longest_protein)
output_file.write("Protein:"+str(longest_protein)+"\n")

output_file.close()