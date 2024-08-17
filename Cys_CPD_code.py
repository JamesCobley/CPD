import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import requests
from io import BytesIO
from itertools import combinations  # Import combinations

def fetch_protein_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta = response.text
        sequence = ''.join(fasta.split('\n')[1:]).replace(' ', '')
        return sequence
    else:
        st.error(f"Failed to fetch data for UniProt ID {uniprot_id}")
        return None

def get_cysteine_positions(sequence):
    """Get positions of cysteine residues (C) in the protein sequence"""
    cysteine_positions = [i for i, res in enumerate(sequence) if res == 'C']
    return cysteine_positions

def generate_heatmap(num_cysteines, cysteine_positions):
    """Generate a heatmap based on the number of cysteines"""
    proteoforms = []
    for i in range(num_cysteines + 1):  # From 0 to num_cysteines cysteines oxidized
        for comb in combinations(cysteine_positions, i):
            proteoform = np.zeros(num_cysteines, dtype=int)
            proteoform[list(comb)] = 1
            proteoforms.append(proteoform)
    
    data = np.array(proteoforms)
    
    # Create a heatmap
    fig, ax = plt.subplots(figsize=(12, 60))
    sns.heatmap(data, cmap="Greys", cbar=False, linewidths=0.5, linecolor='gray', ax=ax)
    ax.set_facecolor('black')  # Set the background color of the cells that have no line
    plt.title("Cysteine Redox Proteoforms", fontsize=20)
    plt.xlabel("Cysteine Sites", fontsize=15)
    plt.ylabel("Proteoforms", fontsize=15)
    plt.xticks(ticks=np.arange(0.5, num_cysteines + 0.5), labels=np.arange(1, num_cysteines + 1), fontsize=12)
    plt.yticks(fontsize=10)
    
    # Save figure to a BytesIO object
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)  # Rewind the buffer to the beginning
    return buf

# Streamlit app
st.title('Cysteine Redox Proteoforms Heatmap')

uniprot_id = st.text_input("Enter UniProt Accession Number:", "P12345")

if uniprot_id:
    sequence = fetch_protein_sequence(uniprot_id)
    if sequence:
        cysteine_positions = get_cysteine_positions(sequence)
        num_cysteines = len(cysteine_positions)

        st.write(f"Protein Sequence Length: {len(sequence)}")
        st.write(f"Number of Cysteines: {num_cysteines}")
        st.write("Cysteine Positions:", cysteine_positions)

        if num_cysteines > 0:
            buf = generate_heatmap(num_cysteines, cysteine_positions)
            
            # Display the plot
            st.image(buf, use_column_width=True, caption='Cysteine Redox Proteoforms Heatmap')

            # Download button
            buf.seek(0)  # Rewind buffer again before providing download link
            st.download_button(
                label="Download Heatmap",
                data=buf,
                file_name="Cysteine_Redox_Proteoforms.png",
                mime="image/png"
            )
