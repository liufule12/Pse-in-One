__author__ = 'Fule Liu'

import time

import const
from util import get_data, check_args, read_k
from pse import get_phyche_list, get_extra_index, get_phyche_value, get_aaindex, extend_aaindex, AAIndex
from data import index_list


def acc(input_data, k, lag, phyche_list, alphabet, extra_index_file=None, all_prop=False, theta_type=1):
    """This is a complete acc in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    :param theta_type: the value 1, 2 and 3 for ac, cc or acc.
    """
    phyche_list = get_phyche_list(k, phyche_list,
                                  extra_index_file=extra_index_file, alphabet=alphabet, all_prop=all_prop)
    # print(phyche_list)
    # Get phyche_vals.
    if alphabet == index_list.DNA or alphabet == index_list.RNA:
        if extra_index_file is not None:
            extra_phyche_index = get_extra_index(extra_index_file)
            from util import normalize_index
            phyche_vals = get_phyche_value(k, phyche_list, alphabet,
                                           normalize_index(extra_phyche_index, alphabet, is_convert_dict=True))
        else:
            phyche_vals = get_phyche_value(k, phyche_list, alphabet)
    elif alphabet == index_list.PROTEIN:
        phyche_vals = get_aaindex(phyche_list)
        # print(phyche_vals)
        if extra_index_file is not None:
            phyche_vals.extend(extend_aaindex(extra_index_file))

    seqs = get_data(input_data, alphabet)
    if alphabet == index_list.PROTEIN:
        # Transform the data format to dict {acid: [phyche_vals]}.
        phyche_keys = phyche_vals[0].index_dict.keys()
        phyche_vals = [e.index_dict.values() for e in phyche_vals]
        new_phyche_vals = zip(*[e for e in phyche_vals])
        phyche_vals = {key: list(val) for key, val in zip(phyche_keys, new_phyche_vals)}

    if theta_type == 1:
        return make_ac_vec(seqs, lag, phyche_vals, k)
    elif theta_type == 2:
        return make_cc_vec(seqs, lag, phyche_vals, k)
    elif theta_type == 3:
        return make_acc_vec(seqs, lag, phyche_vals, k)


def make_ac_vec(sequence_list, lag, phyche_value, k):

    # Get the length of phyche_vals.
    phyche_values = list(phyche_value.values())
    len_phyche_value = len(phyche_values[0])

    vec_ac = []
    for sequence in sequence_list:
        len_seq = len(sequence)
        each_vec = []

        for temp_lag in range(1, lag + 1):
            for j in range(len_phyche_value):

                # Calculate average phyche_value for a nucleotide.
                ave_phyche_value = 0.0
                for i in range(len_seq - temp_lag - k + 1):
                    nucleotide = sequence[i: i + k]
                    ave_phyche_value += float(phyche_value[nucleotide][j])
                ave_phyche_value /= len_seq

                # Calculate the vector.
                temp_sum = 0.0
                for i in range(len_seq - temp_lag - k + 1):
                    nucleotide1 = sequence[i: i + k]
                    nucleotide2 = sequence[i + temp_lag: i + temp_lag + k]
                    temp_sum += (float(phyche_value[nucleotide1][j]) - ave_phyche_value) * (
                        float(phyche_value[nucleotide2][j]))

                each_vec.append(round(temp_sum / (len_seq - temp_lag - k + 1), 8))
        vec_ac.append(each_vec)

    return vec_ac


def make_cc_vec(sequence_list, lag, phyche_value, k):
    phyche_values = list(phyche_value.values())
    len_phyche_value = len(phyche_values[0])

    vec_cc = []
    for sequence in sequence_list:
        len_seq = len(sequence)
        each_vec = []

        for temp_lag in range(1, lag + 1):
            for i1 in range(len_phyche_value):
                for i2 in range(len_phyche_value):
                    if i1 != i2:
                        # Calculate average phyche_value for a nucleotide.
                        ave_phyche_value1 = 0.0
                        ave_phyche_value2 = 0.0
                        for j in range(len_seq - temp_lag - k + 1):
                            nucleotide = sequence[j: j + k]
                            ave_phyche_value1 += float(phyche_value[nucleotide][i1])
                            ave_phyche_value2 += float(phyche_value[nucleotide][i2])
                        ave_phyche_value1 /= len_seq
                        ave_phyche_value2 /= len_seq

                        # Calculate the vector.
                        temp_sum = 0.0
                        for j in range(len_seq - temp_lag - k + 1):
                            nucleotide1 = sequence[j: j + k]
                            nucleotide2 = sequence[j + temp_lag: j + temp_lag + k]
                            temp_sum += (float(phyche_value[nucleotide1][i1]) - ave_phyche_value1) * \
                                        (float(phyche_value[nucleotide2][i2]) - ave_phyche_value2)
                        each_vec.append(round(temp_sum / (len_seq - temp_lag - k + 1), 8))

        vec_cc.append(each_vec)

    return vec_cc


def make_acc_vec(seqs, lag, phyche_values, k):
    from functools import reduce
    zipped = list(zip(make_ac_vec(seqs, lag, phyche_values, k), make_cc_vec(seqs, lag, phyche_values, k)))
    return [reduce(lambda x, y: x + y, e) for e in zipped]


def main(args):
    with open(args.inputfile) as f:
        k = read_k(args.alphabet, args.method, 0)

        # Get index_list.
        if args.i is not None:
            from pse import read_index
            ind_list = read_index(args.i)
        else:
            ind_list = []

        default_e = []
        # Set Pse default index_list.
        if args.alphabet == 'DNA':
            args.alphabet = index_list.DNA
            if k == 2:
                default_e = const.DI_INDS_6_DNA
            elif k == 3:
                default_e = const.TRI_INDS_DNA
        elif args.alphabet == 'RNA':
            args.alphabet = index_list.RNA
            default_e = const.DI_INDS_RNA
        elif args.alphabet == 'Protein':
            args.alphabet = index_list.PROTEIN
            default_e = const.INDS_3_PROTEIN

        theta_type = 1
        if args.method in const.METHODS_AC:
            theta_type = 1
        elif args.method in const.METHODS_CC:
            theta_type = 2
        elif args.method in const.METHODS_ACC:
            theta_type = 3
        else:
            print("Method error!")

        # ACC.
        if args.e is None and len(ind_list) == 0 and args.a is False:
            # Default Pse.
            res = acc(f, k, args.lag, default_e, args.alphabet,
                      extra_index_file=args.e, all_prop=args.a, theta_type=theta_type)
        else:
            res = acc(f, k, args.lag, ind_list, args.alphabet,
                      extra_index_file=args.e, all_prop=args.a, theta_type=theta_type)

    # Write correspond res file.
    if args.f == 'tab':
        from util import write_tab
        write_tab(res, args.outputfile)
    elif args.f == 'svm':
        from util import write_libsvm
        write_libsvm(res, [args.l] * len(res), args.outputfile)
    elif args.f == 'csv':
        from util import write_csv
        write_csv(res, args.outputfile)

    # print(len(res[0]), res[0])


if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is acc module for generate acc vector.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfile',
                       help="The input file, in valid FASTA format.")
    parse.add_argument('outputfile',
                       help="The outputfile stored results.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
                       help="The alphabet of sequences.")
    parse.add_argument('method', type=str,
                       help="The method name of autocorrelation.")

    parse.add_argument('-lag', type=int, default=2,
                       help="The value of lag.")
    parse.add_argument('-i',
                       help="The indices file user choose.\n"
                            "Default indices:\n"
                            "DNA dinucleotide: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                            "DNA trinucleotide: Dnase I, Bendability (DNAse).\n"
                            "RNA: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                            "Protein: Hydrophobicity, Hydrophilicity, Mass.")
    parse.add_argument('-e',
                       help="The user-defined indices file.")
    parse.add_argument('-all_index', dest='a', action='store_true', help="Choose all physicochemical indices")
    parse.add_argument('-no_all_index', dest='a', action='store_false',
                       help="Do not choose all physicochemical indices, default.")
    parse.set_defaults(a=False)
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")
    parse.add_argument('-l', default='+1', choices=['+1', '-1'],
                       help="The libSVM output file label.")

    args = parse.parse_args()
    # print(args)

    if check_args(args, 'acc.py'):
        print("Calculating...")
        start_time = time.time()
        main(args)
        print("Done.")
        print("Used time: %ss" % (time.time() - start_time))

    # # Test ACC for DNA.
    # print("Test ACC for DNA.")
    # print(acc(open('data/test_dna.fasta'), k=2, lag=2, theta_type=3,
    #           phyche_list=['Tilt'], alphabet=index_list.DNA, extra_index_file='data/test_ext_dna.txt'))
    #
    # from repDNA.ac import DACC
    # dacc = DACC(lag=2)
    # print(dacc.make_dacc_vec(open('data/test_dna.fasta'), phyche_index=['Tilt', 'Twist']))
    #
    # print(acc(open('data/test_dna.fasta'), k=3, lag=2, theta_type=3,
    #           phyche_list=['Dnase I'], alphabet=index_list.DNA, extra_index_file='data/test_ext_tridna.txt'))
    #
    # from repDNA.ac import TACC
    # tacc = TACC(lag=2)
    # print(tacc.make_tacc_vec(open('data/test_dna.fasta'), phyche_index=['Dnase I', 'Nucleosome']))
    #
    # # Test ACC for RNA.
    # print("Test ACC for RNA")
    # res = acc(open('data/test_rna.fasta'), k=2, lag=3, theta_type=3,
    #           phyche_list=['Twist (RNA)'], alphabet=index_list.RNA, extra_index_file='data/test_ext_rna.txt')
    # print(len(res[0]), res)
    #
    # # Test ACC for PROTEIN.
    # print("Test ACC for PROTEIN.")
    # res = acc(open('data/test_pro.fasta'), k=1, lag=2, theta_type=3,
    #           phyche_list=['Hydrophobicity', 'Hydrophilicity', 'Mass'], alphabet=index_list.PROTEIN)
    # print(len(res[0]), res)