# Exported from Promoter.ipynb cell 14 (Colab, 28 Feb 2026)
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd
import numpy as np
import plotly.graph_objects as go

# 1. Load the data
sub_df = pd.read_csv(os.path.join(DATA_DIR, 'Reduced_Subsystem_Intersection.csv'))
sys_df = pd.read_csv(os.path.join(DATA_DIR, 'Reduced_System_Intersectioncsv.csv'))

# 2. Filter out requested subsystems
EXCLUDE_SUBSYSTEMS = ["Pre-Golgi", "Post-Golgi"]
sub_df = sub_df[~sub_df['Subsystem'].isin(EXCLUDE_SUBSYSTEMS)]

# 3. Merge to link genes to their System and Subsystem
merged_df = pd.merge(
    sub_df[['CHO GENE SYMBOL', 'Subsystem', 'p-val']],
    sys_df[['CHO GENE SYMBOL', 'System']],
    on='CHO GENE SYMBOL'
)

# 4. Aggregate Subsystems (Outer Ring)
# We use the unique gene count and the best p-value for each subsystem
sub_agg = merged_df.groupby(['System', 'Subsystem']).agg(
    p_val      = ('p-val', 'min'),
    Gene_Count = ('CHO GENE SYMBOL', 'nunique')
).reset_index()

# Handle p=0 to avoid infinity in log scale
sub_agg['p_val'] = sub_agg['p_val'].replace(0, 1e-10).fillna(1.0)
sub_agg['neg_log10_p'] = -np.log10(sub_agg['p_val'])

# 5. Aggregate Systems (Inner Ring)
# IMPORTANT: Gene_Count here must be the SUM of its children to avoid white gaps
sys_agg = sub_agg.groupby('System').agg(
    p_val      = ('p_val', 'min'),
    Gene_Count = ('Gene_Count', 'sum')
).reset_index()
sys_agg['neg_log10_p'] = -np.log10(sys_agg['p_val'].replace(0, 1e-10))

# 6. Build Plotly Hierarchy Lists
ids, labels, parents, values, colors = [], [], [], [], []

# Level 1: System Level (Inner Ring)
for _, row in sys_agg.iterrows():
    ids.append(row['System'])
    labels.append(row['System'])
    parents.append("")  # These connect to the center
    values.append(row['Gene_Count'])
    colors.append(row['neg_log10_p'])

# Level 2: Subsystem Level (Outer Ring)
for _, row in sub_agg.iterrows():
    # ID is System|Subsystem to ensure uniqueness
    ids.append(f"{row['System']}|{row['Subsystem']}")
    labels.append(row['Subsystem'])
    parents.append(row['System'])
    values.append(row['Gene_Count'])
    colors.append(row['neg_log10_p'])

# 7. Define the Smooth Blue Gradient
# From light Ice Blue to Deep Midnight Blue
custom_blues = [
    [0.0, '#E1F5FE'], # Lightest Ice Blue (Low significance)
    [0.2, '#81D4FA'], # Sky Blue
    [0.5, '#0288D1'], # Vivid Blue
    [0.8, '#01579B'], # Deep Blue
    [1.0, '#002244']  # Midnight Blue (High significance)
]

# 8. Create the Plot
fig = go.Figure(go.Sunburst(
    ids=ids,
    labels=labels,
    parents=parents,
    values=values,
    branchvalues="total", # This is what fills the circle completely
    marker=dict(
        colors=colors,
        colorscale=custom_blues,
        colorbar=dict(title="-log₁₀(p)"),
        line=dict(color='white', width=1) # Adds a crisp border to hide gaps
    ),
    hovertemplate='<b>%{label}</b><br>Gene Count: %{value}<br>-log10(p): %{color:.4f}<extra></extra>'
))

fig.update_layout(
    title='System and Subsystem Enrichment Mapping',
    width=1000,
    height=1000,
    margin=dict(t=50, l=10, r=10, b=10)
)

# Radial orientation makes text easier to read on outer rings
fig.update_traces(insidetextorientation='radial', textfont_size=12)

# Show/Save the plot
fig.show()
