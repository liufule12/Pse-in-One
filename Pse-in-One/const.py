__author__ = 'Fule Liu'


METHODS_ALL = ['Kmer', 'RevKmer', 'PseKNC',
               'PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC', 'AC', 'CC', 'ACC']
METHODS_DNA = ['Kmer', 'RevKmer', 'PseKNC', 'PseDNC',
               'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC']
METHODS_RNA = ['Kmer', 'PC-PseDNC-General', 'SC-PseDNC-General', 'DAC', 'DCC', 'DACC']
METHODS_PROTEIN = ['Kmer', 'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General', 'AC', 'CC', 'ACC']

METHODS_AC = ['DAC', 'TAC', 'AC']
METHODS_CC = ['DCC', 'TCC', 'CC']
METHODS_ACC = ['DACC', 'TACC', 'ACC']

K_2_DNA_METHODS = ['PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'DAC', 'DCC', 'DACC']
K_3_DNA_METHODS = ['PC-PseTNC-General', 'SC-PseTNC-General', 'TAC', 'TCC', 'TACC']

THETA_1_METHODS = ['PseDNC', 'PC-PseDNC-General', 'PC-PseTNC-General', 'PC-PseAAC', 'PC-PseAAC-General']
THETA_2_METHODS = ['SC-PseDNC-General', 'SC-PseTNC-General', 'SC-PseAAC', 'SC-PseAAC-General']

DI_INDS_6_DNA = ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist']
TRI_INDS_DNA = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

DI_INDS_RNA = ['Twist (RNA)', 'Tilt (RNA)', 'Roll (RNA)', 'Rise (RNA)', 'Slide (RNA)', 'Shift (RNA)',
               'Stacking energy (RNA)', 'Enthalpy (RNA)1', 'Entropy (RNA)', 'Free energy (RNA)', 'Hydrophilicity (RNA)']
INDS_3_PROTEIN = ['Hydrophobicity', 'Hydrophilicity', 'Mass']
