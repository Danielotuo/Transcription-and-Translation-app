import Bio
from Bio import SeqIO
from Bio.Seq import Seq, transcribe, translate
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from Bio.Data import CodonTable
from collections import Counter
import streamlit as st
from Bio.SeqUtils import seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

st.title("Transcription and Translation")

menu = ["Menu", "Protein Synthesis", "Upload FASTA File"]
choice = st.sidebar.selectbox("Select Activity", menu)

if choice == "Menu":
    st.write("Transcription and translation take the information"
             " in DNA and use it to produce proteins."
             " Translation is the process where the information carried in "
             "mRNA molecules is used to create proteins. The specific "
             "sequence of nucleotides in the mRNA molecule provide the code for"
             "the production of a protein with a specific sequence of amino acids.")
    st.image("dnatr.gif", use_column_width=True)
    st.subheader("Use the sidebar to select activity")

elif choice == "Protein Synthesis":
    st.subheader("Convert DNA to mRNA or to Protein sequence")
    dna_seq = st.text_input("Enter DNA sequence")

    if dna_seq:
        mRNA = transcribe(dna_seq)
        st.write("Transcribed Sequence is: ", mRNA)

        protein_seq = translate(dna_seq)
        st.write("Translated Sequence is: ", protein_seq)

        st.image("dna.jpg", use_column_width=True)

elif choice == "Upload FASTA File":
    st.subheader("Protein sequence")
    sequence_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])

    if sequence_file:
        dna = SeqIO.read(sequence_file, "fasta")

        if st.checkbox("Details of sequence"):
            st.write(dna)

        if st.checkbox("Length of sequence"):
            st.write("length of dna sequence: ", len(dna))

        dna_seq = dna.seq
        translated_dna = dna_seq.translate()

        if st.checkbox("Translation: Each row representing an amino acid sequence"):
            # split amino acids before stop codon. stop codon terminates translation
            AA = translated_dna.split('*')
            protein_sequence = [str(i) for i in AA]
            st.write(protein_sequence)

        # change to 3 letter amino acids instead of 1 letter AA
        if st.checkbox("View 3 letter Amino Acids"):
            st.write(seq3(translated_dna))

        # Amino acid count
        AA_analysed = ProteinAnalysis(str(translated_dna))
        AA_freq = AA_analysed.count_amino_acids()

        if st.checkbox("Amino acid Count"):
            st.write(AA_analysed.get_amino_acids_percent())

        # Visualize the amino acid count
        if st.checkbox("Visualize Amino Acid count"):
            plt.bar(AA_freq.keys(), AA_freq.values())
        st.pyplot()

