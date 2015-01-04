__author__ = 'Fule Liu'


methods_all = ['Kmer', 'RevKmer', 'PseKNC',
               'PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC', 'AC', 'CC', 'ACC']

k_2_DNA_methods = ['PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General']
k_3_DNA_methods = ['PC-PseTNC-General', 'SC-PseTNC-General']

theta_1_methods = ['PseDNC', 'PC-PseDNC-General', 'PC-PseTNC-General', 'PC-PseAAC', 'PC-PseAAC-General']
theta_2_methods = ['SC-PseDNC-General', 'SC-PseTNC-General', 'SC-PseAAC', 'SC-PseAAC-General']

di_inds_6_DNA = ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist']
tri_inds_DNA = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

di_inds_RNA = ['Twist (RNA)', 'Tilt (RNA)', 'Roll (RNA)', 'Rise (RNA)', 'Slide (RNA)', 'Shift (RNA)',
                'Stacking energy (RNA)', 'Enthalpy (RNA)1', 'Entropy (RNA)', 'Free energy (RNA)', 'Hydrophilicity (RNA)']
inds_3_Protein = ['Hydrophobicity', 'Hydrophilicity', 'Mass']
