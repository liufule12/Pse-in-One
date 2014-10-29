__author__ = 'Fule Liu'

import sys
import os
import pickle
from math import pow

from pseALL.util import frequency
from pseALL.util import get_data
from pseALL.kmerutil import make_kmer_list
from pseALL.data import index_list


"""Prepare for PseKNC."""


class AAIndex():
    def __init__(self, head, index_dict):
        self.head = head
        self.index_dict = index_dict

    def __str__(self):
        return "%s\n%s" % (self.head, self.index_dict)


def pseknc(input_data, k, w, lamada, phyche_list, extra_phyche_index, alphabet, all_prop=False, theta_type=1):
    """This is a complete process in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_phyche_index: for DNA and RNA, it is a list. For protein, it is a filename.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    phyche_list = get_phyche_list(k, phyche_list, all_prop=all_prop)
    # Get phyche_vals.
    if alphabet == index_list.DNA or alphabet == index_list.RNA:
        phyche_vals = get_phyche_value(k, phyche_list, extra_phyche_index, alphabet)
    elif alphabet == index_list.PROTEIN:
        phyche_vals = get_aaindex(phyche_list)
        if extra_phyche_index != None:
            phyche_vals.extend(extend_aaindex(extra_phyche_index))
        for e in phyche_vals:
            print(e)

    seq_list = get_data(input_data, alphabet)

    return make_pseknc_vector(seq_list, phyche_vals, k, w, lamada, alphabet, theta_type)


def get_phyche_list(k, phyche_list, all_prop=False):
    """Get phyche_list and check it.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    try:
        if (phyche_list is None or len(phyche_list) == 0) and all_prop is False:
            error_infor = 'Error, The phyche_list and all_prop can\'t be both False'
            raise ValueError(error_infor)
    except:
        raise

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

    # print(all_prop_list)
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


def get_phyche_value(k, phyche_list, extra_phyche_index, alphabet):
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
        theta_sum = sum(theta_list)

        # Generate the vector according the Equation 9.
        denominator = 1 + w * theta_sum

        temp_vec = [round(f / denominator, 3) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 3))

        vector.append(temp_vec)

    return vector


if __name__ == '__main__':
    extra_phyche_index = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1],
                          'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
                          'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
                          'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17, 1],
                          'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
                          'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
                          'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39, 1],
                          'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
                          'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
                          'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59, 1],
                          'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
                          'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
                          'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39, 1],
                          'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
                          'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
                          'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1]}

    # Test dna.
    print("Test dna.")
    alphabet = index_list.DNA
    res = pseknc(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], k=2, w=0.05, lamada=3,
                 phyche_list=['Twist', 'Tilt', 'Roll', 'Rise', 'Slide', 'Shift'], extra_phyche_index=None,
                 alphabet=alphabet)

    for e in res:
        print(e)

    # Test rna.
    print("Test rna.")
    alphabet = index_list.RNA
    res = pseknc(input_data=['GACUGAACUGCACUUUGGUUUCAUAUUAUUUGCUC'], k=2, w=0.05, lamada=3,
                 phyche_list=['Slide (RNA)', 'Adenine content', 'Hydrophilicity (RNA)'], extra_phyche_index=None,
                 alphabet=alphabet, theta_type=2)
    print(res)

    with open('data/test.fasta') as f:
        res = get_data(f, alphabet='ACGT', desc=True)
        for e in res:
            print(e)

    # Test protein.
    alphabet = index_list.PROTEIN
    res = pseknc(input_data=open('data/test.fasta'), k=1, w=0.05, lamada=1,
                 phyche_list=['Hydrophobicity'], extra_phyche_index=None,
                 alphabet=alphabet, theta_type=2)

    for e in res:
        print(len(e), e)

        # For test.
        # from repDNA.psenac import PCPseDNC
        #
        # pcpsednc = PCPseDNC(lamada=3, w=0.05)
        # res = pcpsednc.make_pcpsednc_vec(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'],
        #                                  phyche_index=['Twist', 'Tilt', 'Roll', 'Rise', 'Slide', 'Shift'],
        #                                  extra_phyche_index=extra_phyche_index)
        # print("Test ")
        # for e in res:
        #     print(e)
