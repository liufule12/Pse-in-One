__author__ = 'Fule Liu'


METHODS_ALL = ['Kmer', 'RevKmer', 'PseKNC',
               'PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC', 'AC', 'CC', 'ACC']

METHODS_DNA = ['Kmer', 'RevKmer', 'PseKNC', 'PseDNC',
               'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC']
METHODS_DNA_KMER = ['Kmer', 'RevKmer']
METHODS_DNA_ACC = ['DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC']
METHODS_DNA_PSE = ['PseDNC', 'PseKNC',
                   'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',]

METHODS_RNA = ['Kmer', 'PC-PseDNC-General', 'SC-PseDNC-General', 'DAC', 'DCC', 'DACC']
METHODS_RNA_KMER = ['Kmer']
METHODS_RNA_ACC = ['DAC', 'DCC', 'DACC']
METHODS_RNA_PSE = ['PC-PseDNC-General', 'SC-PseDNC-General']

METHODS_PROTEIN = ['Kmer', 'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General', 'AC', 'CC', 'ACC']
METHODS_PROTEIN_KMER = ['Kmer']
METHODS_PROTEIN_ACC = ['AC', 'CC', 'ACC']
METHODS_PROTEIN_PSE = ['PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General']

KMER_FILENAME = 'kmer.py'
ACC_FILENAME = 'acc.py'
PSE_FILENAME = 'pse.py'

METHODS_AC = ['DAC', 'TAC', 'AC']
METHODS_CC = ['DCC', 'TCC', 'CC']
METHODS_ACC = ['DACC', 'TACC', 'ACC']

K_2_DNA_METHODS = ['PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'DAC', 'DCC', 'DACC']
K_3_DNA_METHODS = ['PC-PseTNC-General', 'SC-PseTNC-General', 'TAC', 'TCC', 'TACC']

THETA_1_METHODS = ['PseDNC', 'PC-PseDNC-General', 'PC-PseTNC-General', 'PC-PseAAC', 'PC-PseAAC-General']
THETA_2_METHODS = ['SC-PseDNC-General', 'SC-PseTNC-General', 'SC-PseAAC', 'SC-PseAAC-General']

DI_INDS_6_DNA = ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist']
TRI_INDS_DNA = ['Dnase I', 'Bendability (DNAse)']

DI_INDS_RNA = ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']
INDS_3_PROTEIN = ['Hydrophobicity', 'Hydrophilicity', 'Mass']
