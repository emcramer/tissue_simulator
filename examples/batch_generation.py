"""
Batch generation example: Generate multiple tissue sections for analysis.
"""

from tissue_simulator import TissueSection
import numpy as np
import matplotlib.pyplot as plt

# Parameters for batch generation
num_sections = 5
tissue_params = {
    'height': 400,
    'width': 400,
    'thickness': 80,
    'cell_radii': {
        'type_a': (5, 10),
        'type_b': (7, 14)
    }
}

# Generate multiple tissue sections
tissues = []
packing_fractions = []
cell_counts = []

print("Generating batch of tissue sections...")
for i in range(num_sections):
    print(f"\nSection {i+1}/{num_sections}")
    
    tissue = TissueSection(**tissue_params)
    num_cells = tissue.generate_cells(max_attempts=1500)
    
    tissues.append(tissue)
    cell_counts.append(num_cells)
    
    stats = tissue.get_cell_statistics()
    packing_fractions.append(stats['packing_fraction'])
    
    print(f"  Cells: {num_cells}")
    print(f"  Packing fraction: {stats['packing_fraction']:.3f}")
    
    # Export each section
    tissue.export_to_csv(f"tissue_section_{i+1}.csv")

# Statistical analysis
print("\n=== Batch Statistics ===")
print(f"Mean cell count: {np.mean(cell_counts):.1f} ± {np.std(cell_counts):.1f}")
print(f"Mean packing fraction: {np.mean(packing_fractions):.3f} ± {np.std(packing_fractions):.3f}")

# Create comparison plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Cell count distribution
axes[0].bar(range(1, num_sections + 1), cell_counts, color='steelblue', alpha=0.7)
axes[0].axhline(np.mean(cell_counts), color='red', linestyle='--', label='Mean')
axes[0].set_xlabel('Section Number')
axes[0].set_ylabel('Cell Count')
axes[0].set_title('Cell Count per Section')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Packing fraction distribution
axes[1].bar(range(1, num_sections + 1), packing_fractions, color='coral', alpha=0.7)
axes[1].axhline(np.mean(packing_fractions), color='red', linestyle='--', label='Mean')
axes[1].set_xlabel('Section Number')
axes[1].set_ylabel('Packing Fraction')
axes[1].set_title('Packing Fraction per Section')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('batch_analysis.png', dpi=150)
print("\nComparison plot saved as 'batch_analysis.png'")
plt.show()
