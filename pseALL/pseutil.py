__author__ = 'Fule Liu'

import sys
import os
import pickle
from math import pow

from pseALL.util import frequency
from pseALL.util import get_data
from pseALL.kmer import make_kmer_list
from pseALL.data import index_list


"""Prepare for PseKNC."""


class AAIndex():
    def __init__(self, head, index_dict):
        self.head = head
        self.index_dict = index_dict

    def __str__(self):
        return "%s\n%s" % (self.head, self.index_dict)


def pseknc(input_data, k, w, lamada, phyche_list, alphabet, extra_index_file=None, all_prop=False, theta_type=1):
    """This is a complete process in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    phyche_list = get_phyche_list(k, phyche_list, extra_index_file=extra_index_file, all_prop=all_prop)

    # Get phyche_vals.
    if alphabet == index_list.DNA or alphabet == index_list.RNA:
        if extra_index_file is not None:
            extra_phyche_index = get_extra_index(extra_index_file)
            from pseALL.util import normalize_index

            phyche_vals = get_phyche_value(k, phyche_list, alphabet,
                                           normalize_index(extra_phyche_index, alphabet, is_convert_dict=True))
        else:
            phyche_vals = get_phyche_value(k, phyche_list, alphabet)
    elif alphabet == index_list.PROTEIN:
        phyche_vals = get_aaindex(phyche_list)
        if extra_index_file is not None:
            phyche_vals.extend(extend_aaindex(extra_index_file))

    seq_list = get_data(input_data, alphabet)

    return make_pseknc_vector(seq_list, phyche_vals, k, w, lamada, alphabet, theta_type)


def ipseknc(input_data, k, w, lamada, phyche_list, alphabet, extra_index_file=None, all_prop=False):
    """This is a complete process in iPseKNC, k is kmer, but the index is just for dinucleotide.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    phyche_list = get_phyche_list(k=2, phyche_list=phyche_list, extra_index_file=extra_index_file, all_prop=all_prop)

    # Get phyche_vals.
    if extra_index_file is not None:
        extra_phyche_index = get_extra_index(extra_index_file)
        from pseALL.util import normalize_index

        phyche_vals = get_phyche_value(k=2, phyche_list=phyche_list, alphabet=alphabet,
                                       extra_phyche_index=normalize_index(extra_phyche_index, alphabet, is_convert_dict=True))
    else:
        phyche_vals = get_phyche_value(k=2, phyche_list=phyche_list, alphabet=alphabet)

    seq_list = get_data(input_data, alphabet)

    return make_pseknc_vector(seq_list, phyche_vals, k, w, lamada, alphabet, theta_type=3)


def get_phyche_list(k, phyche_list, extra_index_file, all_prop=False):
    """Get phyche_list and check it.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    if phyche_list is None or len(phyche_list) == 0:
        if extra_index_file is None and all_prop is False:
            error_info = 'Error, The phyche_list, extra_index_file and all_prop can\'t be all False.'
            raise ValueError(error_info)

    from pseALL.data import index_list

    # Set all_prop_list.
    all_prop_list = []
    try:
        if alphabet == index_list.DNA:
            if k == 2:
                all_prop_list = index_list.didna_list
            elif k == 3:
                all_prop_list = index_list.tridna_list
            else:
                error_info = 'Error, the k value must be 2 or 3.'
                raise ValueError(error_info)
        elif alphabet == index_list.RNA:
            if k == 2:
                all_prop_list = index_list.dirna_list
            else:
                error_info = 'Error, the k or alphabet error.'
                raise ValueError(error_info)
        elif alphabet == index_list.PROTEIN:
            all_prop_list = index_list.pro_list
        else:
            error_info = "Error, the alphabet must be dna, rna or protein."
            raise ValueError(error_info)
    except:
        raise

    # Set and check physicochemical properties.
    try:
        # Set all properties.
        if all_prop is True:
            phyche_list = all_prop_list
        # Check phyche properties.
        else:
            for e in phyche_list:
                if e not in all_prop_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    raise NameError(error_info)
    except:
        raise

    return phyche_list


def get_extra_index(filename):
    extra_index_vals = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip().split('\t')
            vals = [float(val) for val in line[1:]]
            extra_index_vals.append(vals)

    return extra_index_vals


def get_aaindex(index_list):
    """Get the aaindex from data/aaindex.data.

    :param index_list: the index we want to get.
    :return: a list of AAIndex obj.
    """
    new_aaindex = []
    with open('data/aaindex.data', 'rb') as f:
        aaindex = pickle.load(f)
        for index_vals in aaindex:
            if index_vals.head in index_list:
                new_aaindex.append(index_vals)

    return new_aaindex


def extend_aaindex(filename):
    """Extend the user-defined AAIndex from user's file.
    :return: a list of AAIndex obj.
    """
    from pseALL.scrip.extract_aaindex import extra_aaindex, norm_index_vals

    aaindex = extra_aaindex(filename)
    for ind, e in enumerate(aaindex):
        aaindex[ind] = AAIndex(e.head, norm_index_vals(e.index_dict))

    return aaindex


def get_phyche_value(k, phyche_list, alphabet, extra_phyche_index=None):
    """Generate DNA or RNA phyche_value.

    :param k: int, the value of k-tuple.
    :param phyche_list: physicochemical properties list.
    :param extra_phyche_index: dict, the key is the olinucleotide (string),
                                     the value is its physicochemical property value (list).
                               It means the user-defined physicochemical indices.
    """
    if extra_phyche_index is None:
        extra_phyche_index = {}

    phyche_value = extend_phyche_index(get_phyche_index(k, phyche_list, alphabet), extra_phyche_index)

    return phyche_value


def extend_phyche_index(original_index, extend_index):
    """Extend DNA or RNA {phyche:[value, ... ]}"""
    if extend_index is None or len(extend_index) == 0:
        return original_index
    for key in list(original_index.keys()):
        original_index[key].extend(extend_index[key])
    return original_index


def get_phyche_factor_dic(k, alphabet):
    """Get all DNA or RNA {nucleotide: [(phyche, value), ...]} dict."""
    full_path = os.path.realpath(__file__)
    if 2 == k and alphabet == index_list.DNA:
        file_path = "%s/data/didna.data" % os.path.dirname(full_path)
    elif 2 == k and alphabet == index_list.RNA:
        file_path = "%s/data/dirna.data" % os.path.dirname(full_path)
    elif 3 == k:
        file_path = "%s/data/mmc4.data" % os.path.dirname(full_path)
    else:
        sys.stderr.write("The k can just be 2 or 3.")
        sys.exit(0)

    try:
        with open(file_path, 'rb') as f:
            phyche_factor_dic = pickle.load(f)
    except:
        with open(file_path, 'r') as f:
            phyche_factor_dic = pickle.load(f)

    return phyche_factor_dic


def get_phyche_index(k, phyche_list, alphabet):
    """get phyche_value according phyche_list."""
    phyche_value = {}
    if 0 == len(phyche_list):
        for nucleotide in make_kmer_list(k, alphabet):
            phyche_value[nucleotide] = []
        return phyche_value

    nucleotide_phyche_value = get_phyche_factor_dic(k, alphabet)
    for nucleotide in make_kmer_list(k, alphabet):
        if nucleotide not in phyche_value:
            phyche_value[nucleotide] = []
        for e in nucleotide_phyche_value[nucleotide]:
            if e[0] in phyche_list:
                phyche_value[nucleotide].append(e[1])

    return phyche_value


"""Calculate PseKNC."""


def parallel_cor_function(nucleotide1, nucleotide2, phyche_index):
    """Get the cFactor.(Type1)"""
    temp_sum = 0.0
    phyche_index_values = list(phyche_index.values())
    len_phyche_index = len(phyche_index_values[0])
    for u in range(len_phyche_index):
        temp_sum += pow(float(phyche_index[nucleotide1][u]) - float(phyche_index[nucleotide2][u]), 2)

    return temp_sum / len_phyche_index


def series_cor_function(nucleotide1, nucleotide2, big_lamada, phyche_value):
    """Get the series correlation Factor(Type 2)."""
    return float(phyche_value[nucleotide1][big_lamada]) * float(phyche_value[nucleotide2][big_lamada])


def pro_cor_fun1(ri, rj, aaindex_list):
    _sum = 0.0
    len_index = len(aaindex_list)
    for aaindex in aaindex_list:
        _sum += pow(aaindex.index_dict[ri] - aaindex.index_dict[rj], 2)
    return _sum / len_index


def pro_cor_fun2(ri, rj, aaindex):
    return aaindex.index_dict[ri] * aaindex.index_dict[rj]


def get_parallel_factor(k, lamada, sequence, phyche_value, alphabet):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)

    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - k - i + 1):
            nucleotide1 = sequence[j: j + k]
            nucleotide2 = sequence[j + i: j + i + k]
            if alphabet == index_list.DNA or alphabet == index_list.RNA:
                temp_sum += parallel_cor_function(nucleotide1, nucleotide2, phyche_value)
            elif alphabet == index_list.PROTEIN:
                temp_sum += pro_cor_fun1(nucleotide1, nucleotide2, phyche_value)

        theta.append(temp_sum / (l - k - i + 1))

    return theta


def get_series_factor(k, lamada, sequence, phyche_value, alphabet):
    """Get the corresponding series factor theta list."""
    theta = []
    l_seq = len(sequence)
    if alphabet == index_list.DNA or alphabet == index_list.RNA:
        temp_values = list(phyche_value.values())
        max_big_lamada = len(temp_values[0])
    elif alphabet == index_list.PROTEIN:
        max_big_lamada = len(phyche_value)

    for small_lamada in range(1, lamada + 1):
        for big_lamada in range(max_big_lamada):
            temp_sum = 0.0
            for i in range(0, l_seq - k - small_lamada + 1):
                nucleotide1 = sequence[i: i + k]
                nucleotide2 = sequence[i + small_lamada: i + small_lamada + k]
                if alphabet == index_list.DNA or alphabet == index_list.RNA:
                    temp_sum += series_cor_function(nucleotide1, nucleotide2, big_lamada, phyche_value)
                elif alphabet == index_list.PROTEIN:
                    temp_sum += pro_cor_fun2(nucleotide1, nucleotide2, phyche_value[big_lamada])

            theta.append(temp_sum / (l_seq - k - small_lamada + 1))

    return theta


def make_pseknc_vector(sequence_list, phyche_value, k=2, w=0.05, lamada=1, alphabet=index_list.DNA, theta_type=1):
    """Generate the pseknc vector."""
    kmer = make_kmer_list(k, alphabet)
    vector = []

    for sequence in sequence_list:
        if len(sequence) < k or lamada + k > len(sequence):
            error_info = "Sorry, the sequence length must be larger than " + str(lamada + k)
            sys.stderr.write(error_info)
            sys.exit(0)

        # Get the nucleotide frequency in the DNA sequence.
        fre_list = [frequency(sequence, str(key)) for key in kmer]
        fre_sum = float(sum(fre_list))

        # Get the normalized occurrence frequency of nucleotide in the DNA sequence.
        fre_list = [e / fre_sum for e in fre_list]

        # Get the theta_list.
        if 1 == theta_type:
            theta_list = get_parallel_factor(k, lamada, sequence, phyche_value, alphabet)
        elif 2 == theta_type:
            theta_list = get_series_factor(k, lamada, sequence, phyche_value, alphabet)
        elif 3 == theta_type:
            theta_list = get_parallel_factor(k=2, lamada=lamada, sequence=sequence,
                                             phyche_value=phyche_value, alphabet=alphabet)
        theta_sum = sum(theta_list)

        # Generate the vector according the Equation 9.
        denominator = 1 + w * theta_sum

        temp_vec = [round(f / denominator, 3) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 3))

        vector.append(temp_vec)

    return vector


if __name__ == '__main__':


    # Test dna type1.
    print("Test di_dna, type1.")
    alphabet = index_list.DNA
    res = pseknc(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], k=2, w=0.5, lamada=10,
                 phyche_list=['Tilt', 'Roll', 'Rise', 'Slide', 'Shift'],
                 extra_index_file="data/test_ext_dna.txt", alphabet=alphabet)

    for e in res:
        print(e)

    print("Test tri_dna, type1.")
    alphabet = index_list.DNA
    res = pseknc(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], k=3, w=0.5, lamada=10,
                 phyche_list=['Dnase I'],
                 extra_index_file="data/test_ext_tridna.txt", alphabet=alphabet)

    for e in res:
        print(e)

    # Test dna type2
    print("Test di_dna, type2.")
    alphabet = index_list.DNA
    res = pseknc(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], k=2, w=0.5, lamada=2,
                 phyche_list=['Tilt', 'Roll', 'Twist'],
                 extra_index_file=None, alphabet=alphabet, theta_type=2)

    for e in res:
        print(len(e), e)

    print("Test tri_dna, type2.")
    alphabet = index_list.DNA
    res = pseknc(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], k=3, w=0.5, lamada=2,
                 phyche_list=['Dnase I'],
                 extra_index_file="data/test_ext_tridna.txt", alphabet=alphabet, theta_type=2)

    for e in res:
        print(e)

    # Test iPseKNC.
    print("Test iPseKNC.")
    alphabet = index_list.DNA
    res = ipseknc(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], k=3, w=0.5, lamada=10,
                  phyche_list=['Tilt', 'Roll', 'Rise', 'Slide', 'Shift'],
                  extra_index_file="data/test_ext_dna.txt", alphabet=alphabet)

    for e in res:
        print(len(e), e)

    # Test rna.
    default_indexs = ['Twist (RNA)', 'Tilt (RNA)', 'Roll (RNA)', 'Rise (RNA)', 'Slide (RNA)', 'Shift (RNA)',
                      'Stacking energy (RNA)', 'Enthalpy (RNA)1', 'Entropy (RNA)', 'Free energy (RNA)',
                      'Hydrophilicity (RNA)']

    _default_indexs = ['Twist (RNA)', 'Tilt (RNA)', 'Roll (RNA)', 'Rise (RNA)', 'Slide (RNA)', 'Shift (RNA)',
                       'Stacking energy (RNA)', 'Enthalpy (RNA)1', 'Entropy (RNA)', 'Free energy (RNA)']

    print("Test rna, type1.")
    alphabet = index_list.RNA
    res = pseknc(input_data=['GACUGAACUGCACUUUGGUUUCAUAUUAUUUGCUC'], k=2, w=0.05, lamada=3,
                 phyche_list=_default_indexs, extra_index_file="data/test_ext_rna.txt",
                 alphabet=alphabet, theta_type=1)
    print(res)

    print("Test rna, type2.")
    alphabet = index_list.RNA
    res = pseknc(input_data=['GACUGAACUGCACUUUGGUUUCAUAUUAUUUGCUC'], k=2, w=0.05, lamada=3,
                 phyche_list=['Slide (RNA)', 'Adenine content', 'Hydrophilicity (RNA)'], extra_index_file=None,
                 alphabet=alphabet, theta_type=2)
    print(res)

    with open('data/test.fasta') as f:
        res = get_data(f, alphabet='ACGT', desc=True)
        for e in res:
            print(e)

    # Test protein.
    default_pro = ['Hydrophobicity', 'Hydrophilicity', 'Mass']
    alphabet = index_list.PROTEIN
    res = pseknc(input_data=open('data/test.fasta'), k=1, w=0.05, lamada=2,
                 phyche_list=['Hydrophobicity', 'Hydrophilicity'], extra_index_file="data/test_ext_pro.txt",
                 alphabet=alphabet, theta_type=1)

    for e in res:
        print(len(e), e)