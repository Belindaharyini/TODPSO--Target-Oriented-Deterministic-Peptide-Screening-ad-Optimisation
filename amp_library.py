# -*- coding: utf-8 -*-
import random
import pandas as pd

# Parameters
batch_size = 10_000_000
top_per_batch = 1000
total_batches = 10
seed = 56

# Amino acid groups
negatively_charged = ['D', 'E']
positively_charged = ['R', 'K']
hydrophobic = ['A', 'G','P', 'V', 'L', 'I', 'M', 'F', 'W', 'Y']
polar = ['S', 'T', 'N', 'Q']
other = ['C', 'H']
all_amino_acids = negatively_charged + positively_charged + hydrophobic + polar + other

# Core Motif
core_motif = ['D', 'X', 'E', 'X', 'X', 'X', 'X', 'R', 'X', 'X', 'X', 'X', 'K', 'X', 'X', 'X', 'X', 'X', 'X', 'X']

# Scoring Function
def calculate_properties(peptide):
    hydrophobicity = sum(1 for aa in peptide if aa in hydrophobic) / len(peptide)
    net_charge = sum(1 for aa in peptide if aa in positively_charged) - sum(1 for aa in peptide if aa in negatively_charged)
    protease_resistance = peptide.endswith(('G', 'P', 'K'))
    aggregation_risk = sum(1 for i in range(len(peptide) - 2) if peptide[i] in hydrophobic and peptide[i+1] in hydrophobic)

    score = (
        (0.3 <= hydrophobicity <= 0.5) * 40 +
        (abs(net_charge) <= 3) * 30 +
        protease_resistance * 20 +
        (aggregation_risk <= 2) * 10
    )
    return hydrophobicity, net_charge, protease_resistance, aggregation_risk, score

# Deterministic Peptide Generator
def generate_batch(batch_id):
    random.seed(seed + batch_id)
    peptide_library = []
    for i in range(batch_size):
        peptide = list(core_motif)
        for pos, aa in enumerate(peptide):
            if aa == 'X':
                peptide[pos] = random.choice(all_amino_acids)
        peptide = ''.join(peptide)

        # Calculate properties
        hydrophobicity, net_charge, protease_resistance, aggregation_risk, score = calculate_properties(peptide)

        peptide_library.append({
            "ID": f"B{batch_id:02}_PEP_{i+1:07}",
            "Peptide": peptide,
            "Hydrophobicity": hydrophobicity,
            "Net_Charge": net_charge,
            "Protease_Resistance": protease_resistance,
            "Aggregation_Risk": aggregation_risk,
            "Score": score
        })

    # Filter Top Peptides per Batch
    batch_df = pd.DataFrame(peptide_library)
    top_batch_df = batch_df.sort_values(by='Score', ascending=False).head(top_per_batch)
    top_batch_df.to_csv(f'top_{top_per_batch}_peptides_batch_{batch_id}.csv', index=False)
    print(f"Batch {batch_id} processed and saved.")
    return top_batch_df

all_top_peptides = []
for batch_id in range(1, total_batches + 1):
    print(f"Processing Batch {batch_id}...")
    top_batch_df = generate_batch(batch_id)
    all_top_peptides.append(top_batch_df)

# Filtered 1000 peptides
final_df = pd.concat(all_top_peptides)
final_top_1000 = final_df.sort_values(by='Score', ascending=False).head(1000)
final_top_1000.to_csv('final_top_1000_peptides.csv', index=False)

print("Final Top 1000 peptides saved as 'final_top_1000_peptides.csv'")
