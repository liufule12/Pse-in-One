# coding=utf-8
__author__ = 'Fule Liu'


import sys
import os
import pickle

from pseALL.kmer import make_kmer_list


def convert_phyche_index_to_dict(phyche_index, alphabet):
    """Convert phyche index from list to dict."""
    # for e in phyche_index:
    #     print e
    len_index_value = len(phyche_index[0])
    k = 0
    for i in range(1, 10):
        if len_index_value < 4**i:
            error_infor = 'Sorry, the number of each index value is must be 4^k.'
            sys.stdout.write(error_infor)
            sys.exit(0)
        if len_index_value == 4**i:
            k = i
            break
    from repDNA.nacutil import make_kmer_list
    kmer_list = make_kmer_list(k, alphabet)
    # print kmer_list
    len_kmer = len(kmer_list)
    phyche_index_dict = {}
    for kmer in kmer_list:
        phyche_index_dict[kmer] = []
    # print phyche_index_dict
    phyche_index = list(zip(*phyche_index))
    for i in range(len_kmer):
        phyche_index_dict[kmer_list[i]] = list(phyche_index[i])

    return phyche_index_dict


def standard_deviation(value_list):
    """Return standard deviation."""
    from math import sqrt
    from math import pow
    n = len(value_list)
    average_value = sum(value_list) * 1.0 / n
    return sqrt(sum([pow(e - average_value, 2) for e in value_list]) * 1.0 / (n - 1))


def normalize_index(phyche_index, is_convert_dict=False):
    """Normalize the physicochemical index."""
    normalize_phyche_value = []
    for phyche_value in phyche_index:
        average_phyche_value = sum(phyche_value) * 1.0 / len(phyche_value)
        sd_phyche = standard_deviation(phyche_value)
        normalize_phyche_value.append([round((e - average_phyche_value) / sd_phyche, 2) for e in phyche_value])

    if is_convert_dict is True:
        return convert_phyche_index_to_dict(normalize_phyche_value)

    return normalize_phyche_value


def add_property_id(property_name, property_dict, property_value):
    """This function is for function read_index_file.
    """
    for i in range(1, 100):
        temp_property = property_name + str(i)
        if temp_property not in property_dict:
            property_dict[temp_property] = property_value
            return property_dict


def read_index_file(filename):
    """Read DiProDB file, extra the DNA and RNA index_vals dict."""
    with open(filename) as f:
        lines = f.readlines()
        dna_dict, rna_dict = {}, {}
        dna_nucleic_acid = ['B-DNA', 'DNA', 'DNA/RNA']
        rna_nucleic_acid = ['A-RNA', 'RNA', 'DNA/RNA']

        for line in lines[1:]:
            line = line.rstrip().split('\t')
            nucleic_acid = line[-1]
            property_name = line[1]
            property_value = [float(e) for e in line[2:-1]]

            # Add a property index in DNA.
            if nucleic_acid in dna_nucleic_acid:
                if property_name in dna_dict:
                    dna_dict = add_property_id(property_name, dna_dict, property_value)
                else:
                    dna_dict[property_name] = property_value

            # Add a property index in RNA.
            if nucleic_acid in rna_nucleic_acid:
                if property_name in rna_dict:
                    rna_dict = add_property_id(property_name, rna_dict, property_value)
                else:
                    rna_dict[property_name] = property_value

    dna_res = dna_dict.items()
    rna_res = rna_dict.items()
    # for e in dna_res:
    #     print(e)
    # print(len(dna_res))
    # print(len(rna_res))

    return dna_dict, rna_dict


def combine_dna_dict(dna_dict, alphabet, write_file):
    tri_dna_file = os.path.abspath('..') + '/data/mmc3.data'
    with open(tri_dna_file, 'rb') as f:
        dinucle_property_val = pickle.load(f)
        kmer_list = make_kmer_list(k=2, alphabet=alphabet)

        for _property, vals in dna_dict.items():
            for ind, val in enumerate(vals):
                for e in dinucle_property_val[kmer_list[ind]]:
                    if e[0] == _property:
                        print(kmer_list[ind], e, _property)
                        break
                dinucle_property_val[kmer_list[ind]].append((_property, val))

        with open(write_file, 'wb') as f_write:
            pickle.dump(dinucle_property_val, f_write, protocol=2)

    return dinucle_property_val


def write_rna(rna_dict, filename):
    kmer_list = make_kmer_list(k=2, alphabet='ACGU')

    dinucle_prop_vals = {}
    for _property, vals in rna_dict.items():
        for ind, val in enumerate(vals):
            dinucle = kmer_list[ind]
            if dinucle not in dinucle_prop_vals:
                dinucle_prop_vals[dinucle] = []
            dinucle_prop_vals[dinucle].append((_property, val))

    print(len(dinucle_prop_vals['AA']), dinucle_prop_vals)

    with open(filename, 'wb') as f_write:
        pickle.dump(dinucle_prop_vals, f_write, protocol=2)


if __name__ == '__main__':
    # 读取DiProGB数据库DNA、RNA数据，存储为字典。
    father_path = os.path.abspath('..')

    dna_dict, rna_dict = read_index_file(father_path + '/data/diindex_name.txt')

    # 归一化，以字典形式存储与didna.data文件。
    dinucle_name, dinucle_vals = dna_dict.keys(), dna_dict.values()
    print(dinucle_name)
    print(len(dinucle_name))
    print(rna_dict.keys())
    print(len(rna_dict.keys()))
    dinucle_vals = normalize_index(dinucle_vals)
    for ind, name in enumerate(dinucle_name):
        dna_dict[name] = dinucle_vals[ind]

    # 合并原mmc3.data与didna.data文件。
    combine_dna_dict(dna_dict, 'ACGT', father_path + '/data/didna.data')

    # 写2核酸rna索引文件
    write_rna(rna_dict, father_path + '/data/dirna.data')

    with open(father_path + '/data/didna.data', 'rb') as f:
        res = pickle.load(f)
        for e in res.items():
            print(len(e[1]), e)

    with open(father_path + '/data/dirna.data', 'rb') as f:
        res = pickle.load(f)
        for e in res.items():
            print(len(e[1]), e)
