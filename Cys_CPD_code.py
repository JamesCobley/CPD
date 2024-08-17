import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from io import BytesIO

# Initialize parameters
num_cysteines = 10  # Total number of cysteines

# Generate all unique proteoforms
proteoforms = []
for i in range(num_cysteines + 1):  # From 0 to 10 cysteines oxidized
    for comb in combinations(range(num_cysteines), i):
        proteoform = np.zeros(num_cysteines, dtype=int)
        proteoform[list(comb)] = 1
        proteoforms.append(proteoform)

# Convert list of arrays to a 2D NumPy array
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

# Streamlit app
st.title('Cysteine Redox Proteoforms Heatmap')
st.write('This heatmap visualizes the various proteoforms based on the number of cysteines oxidized.')

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
