"""Python analogue of develop_BOM.m

Reads a blade YAML, defines blade segments (spanwise and/or chord-wise),
and exports a segmented Bill-of-Materials to Excel for TEA cost analysis.
"""
import os
from pynumad.objects.blade import Blade

yaml_file = os.path.join(
    os.path.dirname(__file__), "example_data", "IEA-15MW_modified.yaml"
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
