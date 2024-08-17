import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import requests
from io import BytesIO
from itertools import combinations
from scipy.special import comb  # Import comb for binomial coefficients

def fetch_protein_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    try:
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)
        response.raise_for_status()  # Raises an HTTPError for bad responses
        fasta = response.text
        sequence = ''.join(fasta.split('\n')[1:]).replace(' ', '')
        return sequence
    except requests.RequestException as e:
        st.error(f"Failed to fetch data for UniProt ID {uniprot_id}: {e}")
        return None

def get_cysteine_positions(sequence):
    """Get positions of cysteine residues (C) in the protein sequence"""
    return [i + 1 for i, res in enumerate(sequence) if res == 'C']  # Adjust to start from 1

def generate_heatmap(cysteine_positions):
    """Generate a heatmap based on the cysteine positions"""
    num_cysteines = len(cysteine_positions)
    proteoforms = []
    
    # Generate all unique proteoforms
    for i in range(num_cysteines + 1):
        for combi in combinations(range(num_cysteines), i):
            proteoform = np.zeros(num_cysteines, dtype=int)
            proteoform[list(combi)] = 1
            proteoforms.append(proteoform)
    
    data = np.array(proteoforms)
    
    # Calculate Pascal's Triangle row
    pascal_row = [int(comb(num_cysteines, k, exact=True)) for k in range(num_cysteines + 1)]
    
    # Calculate percentage redox grades
    redox_grades = [(100 * (num_cysteines - k) / num_cysteines) for k in range(num_cysteines + 1)]
    
    # Combine Pascal row with redox percentages
    pascal_with_percentages = [f"{value} ({percentage:.1f}%)" for value, percentage in zip(pascal_row, redox_grades)]
    
    # Number of redox grades
    num_redox_grades = len(pascal_row)
    
    st.write(f"Number of % redox grades: {num_redox_grades}")
    st.write(f"i space: {len(proteoforms)}")
    st.write(f"Pascal's Triangle Row for {num_cysteines} Cysteines with Redox Percentages: {pascal_with_percentages}")
    
    # Create a heatmap
    fig, ax = plt.subplots(figsize=(12, 60))
    sns.heatmap(data, cmap="Greys", cbar=False, linewidths=0.5, linecolor='gray', ax=ax)
    ax.set_facecolor('black')
    plt.title("Cysteine Redox Proteoforms", fontsize=20)
    plt.xlabel("Cysteine Sites", fontsize=15)
    plt.ylabel("Proteoforms", fontsize=15)
    plt.xticks(ticks=np.arange(0.5, num_cysteines + 0.5), labels=cysteine_positions, fontsize=12)
    plt.yticks(fontsize=10)
    
    # Save figure to a BytesIO object
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)  # Rewind the buffer to the beginning
    return buf

# Streamlit app
st.title('Cysteine Redox Proteoforms Diagram')

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
            buf = generate_heatmap(cysteine_positions)
            
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
