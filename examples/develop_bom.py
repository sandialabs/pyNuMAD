"""Python analogue of develop_BOM.m

Reads a blade YAML, defines blade segments (spanwise and/or chord-wise),
and exports a segmented Bill-of-Materials to Excel for TEA cost analysis.
"""
import os
import numpy as np
import pandas as pd
from pynumad.objects.blade import Blade
from pynumad.graphics.graphics import plot_blade_geometry_pv, plot_segment_mass

yaml_file = os.path.join(
    os.path.dirname(__file__), "example_data", "IEA-22-280-RWT.yaml"
)

# Choose a segmentation case: 'one-piece', 'two-piece', 'chord-wise', 'three-piece'
segment_case = "chord-wise"

if segment_case == "one-piece":
    segments = [
        {"span_nd": [0.0, 1.0], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
    ]
elif segment_case == "two-piece":
    segments = [
        {"span_nd": [0.0, 0.5], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
        {"span_nd": [0.5, 1.0], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
    ]
elif segment_case == "chord-wise":
    segments = [
        {"span_nd": [0.0, 0.5], "hp_extents": ["le", "b"], "lp_extents": ["le", "b"]},
        {"span_nd": [0.0, 0.5], "hp_extents": ["b", "te"], "lp_extents": ["b", "te"]},
        {"span_nd": [0.5, 1.0], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
    ]
elif segment_case == "three-piece":
    segments = [
        {"span_nd": [0.00, 0.33], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
        {"span_nd": [0.33, 0.67], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
        {"span_nd": [0.67, 1.00], "hp_extents": ["le", "te"], "lp_extents": ["le", "te"]},
    ]

blade = Blade(yaml_file)
blade.update_blade()

output_file = f"bom_TEA_{segment_case}.xlsx"
sheet_name = "IEA-15 BOM"

bom_tea = blade.export_bom_tea(output_file, sheet_name, segments)

csv_file = output_file.replace(".xlsx", ".csv")
bom_tea.to_csv(csv_file, index=False)

print(f"BOM written to {output_file} and {csv_file} ({len(bom_tea)} rows)")
print(bom_tea.head(20).to_string(index=False))

bom = blade.bill_of_materials
g_to_kg = 0.001
print("\nMass by segment:")
for df, label in [(bom.hp, "HP"), (bom.lp, "LP"), (bom.sw, "SW")]:
    if df.empty:
        continue
    for seg_id, grp in df.groupby("segment_id"):
        mass_kg = grp["weight"].sum() * g_to_kg
        print(f"  {label} {seg_id}: {mass_kg:.2f} kg")
print(f"Total dry weight: {bom.dryweight:.2f} kg")

g_to_kg = 0.001
ispan_vals = blade.ispan

def _distribute(df, ispan):
    density = np.zeros(len(ispan))
    for _, row in df.iterrows():
        span_len = row["endsta"] - row["beginsta"]
        if span_len <= 0:
            continue
        mask = (ispan >= row["beginsta"]) & (ispan <= row["endsta"])
        density[mask] += row["weight"] * g_to_kg / span_len
    return density

hp_d = _distribute(bom.hp, ispan_vals) if not bom.hp.empty else np.zeros(len(ispan_vals))
lp_d = _distribute(bom.lp, ispan_vals) if not bom.lp.empty else np.zeros(len(ispan_vals))
sw_d = _distribute(bom.sw, ispan_vals) if not bom.sw.empty else np.zeros(len(ispan_vals))
total_d = hp_d + lp_d + sw_d

mass_dist = pd.DataFrame({
    "span_m":     ispan_vals,
    "hp_kg_m":    hp_d,
    "lp_kg_m":    lp_d,
    "sw_kg_m":    sw_d,
    "total_kg_m": total_d,
})
print("\nSpanwise linear mass density (kg/m):")
print(mass_dist.to_string(index=False, float_format="%.3f"))

plot_segment_mass(bom, segments)
import matplotlib.pyplot as plt
plt.show()

plot_blade_geometry_pv(blade)
