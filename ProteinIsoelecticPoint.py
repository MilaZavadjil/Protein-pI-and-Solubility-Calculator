import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

#Import protein sequence
from pyfaidx import Fasta
sequences = Fasta('uniprot_sprot.fasta')
protein_sequence = sequences['sp|P10636|TAU_HUMAN']
print(protein_sequence)

#Isoelectric point calculator
protein_sequence = str(protein_sequence)
protein = IP(protein_sequence)
result_pI = protein.pi()
print(result_pI)

# Percent Solubility Calculator
def nonlinear_percent_solubility(pH, pI, a, b):
    if pH < pI:
        # Protein has a net positive charge
        return 100 * np.exp(a * (pI - pH)) / (1 + np.exp(a * (pI - pH)))
    elif pH > pI:
        # Protein has a net negative charge
        return 100 * np.exp(b * (pH - pI)) / (1 + np.exp(b * (pH - pI)))
    else:
        # pH is at the pI, minimal net charge
        return 50  # 50% solubility at the pI

# Example values
pI = result_pI
a = 1.2
b = 1.2
pH_values = np.linspace(-1, 14, 100)

# Calculate estimated percent solubility at each pH using a nonlinear function
percent_solubility_values = [nonlinear_percent_solubility(pH, pI, a, b) for pH in pH_values]

# Plot percent solubility as a function of pH
plt.plot(pH_values, percent_solubility_values)
plt.xlabel('pH')
plt.ylabel('Percent Solubility')
plt.title('Nonlinear Percent Solubility of Protein as a Function of pH')
plt.grid(True)
plt.show()