__author__ = 'Fule Liu'

import sys
import re

from util import frequency
from util import get_data
from data import index_list

sys.setrecursionlimit(99999999)


def make_kmer_list(k, alphabet):
    # Base case.
    if k == 1:
        return alphabet

    # Handle k=0 from user.
    if k == 0:
        return []

    # Error case.
    if k < 1:
        sys.stderr.write("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = [kmer + alphabet[i_letter] for kmer in make_kmer_list(k - 1, alphabet)
                    for i_letter in range(alphabet_length)]

    return return_value


def find_revcomp(sequence, revcomp_dictionary):
    # Save time by storing reverse complements in a hash.
    if sequence in revcomp_dictionary:
        return revcomp_dictionary[sequence]

    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)

    return_value = ""
    for letter in rev_sequence:
        if letter == "A":
            return_value += "T"
        elif letter == "C":
            return_value += "G"
        elif letter == "G":
            return_value += "C"
        elif letter == "T":
            return_value += "A"
        elif letter == "N":
            return_value += "N"
        else:
            error_info = ("Unknown DNA character (%s)\n" % letter)
            sys.exit(error_info)

    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value

    return return_value


def _cmp(a, b):
    return (a > b) - (a < b)


def make_revcomp_kmer_list(kmer_list):
    revcomp_dictionary = {}
    new_kmer_list = [kmer for kmer in kmer_list if _cmp(kmer, find_revcomp(kmer, revcomp_dictionary)) <= 0]
    return new_kmer_list


def make_kmer_vector(k, alphabet, filename, revcomp=False):
    """Generate kmer vector."""
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)

        if revcomp and re.search(r'[^acgtACGT]', ''.join(alphabet)) is not None:
            sys.exit("Error, Only DNA sequence can be reverse compliment.")

        vector = []
        kmer_list = make_kmer_list(k, alphabet)
        for seq in seq_list:
            count_sum = 0

            # Generate the kmer frequency dict.
            kmer_count = {}
            for kmer in kmer_list:
                temp_count = frequency(seq, kmer)
                if not revcomp:
                    if kmer not in kmer_count:
                        kmer_count[kmer] = 0
                    kmer_count[kmer] += temp_count
                else:
                    rev_kmer = find_revcomp(kmer, {})
                    if kmer <= rev_kmer:
                        if kmer not in kmer_count:
                            kmer_count[kmer] = 0
                        kmer_count[kmer] += temp_count
                    else:
                        if rev_kmer not in kmer_count:
                            kmer_count[rev_kmer] = 0
                        kmer_count[rev_kmer] += temp_count

                count_sum += temp_count

            # Normalize.
            if not revcomp:
                count_vec = [kmer_count[kmer] for kmer in kmer_list]
            else:
                revc_kmer_list = make_revcomp_kmer_list(kmer_list)
                count_vec = [kmer_count[kmer] for kmer in revc_kmer_list]
            count_vec = [round(float(e)/count_sum, 3) for e in count_vec]

            vector.append(count_vec)

    return vector


if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is a kmer module for generate kmer vector.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfile',
                       help="The input file, in valid FASTA format.")
    parse.add_argument('outputfile',
                       help="The outputfile stored results.")
    parse.add_argument('k', type=int, choices=range(1, 7),
                       help="The k value of kmer.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'PROTEIN'],
                       help="The alphabet of sequences.")
    parse.add_argument('-r', default=0, type=int, choices=[1, 0],
                       help="Whether need to reverse complement.\n"
                            "1 means True, 0 means False. (default = 0)")
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")

    args = parse.parse_args()

    # Set revcomp parameter.
    if args.r != 1:
        args.r = False
    elif args.r == 1 and args.alphabet != 'DNA':
        print("Error, the -r parameter can only be used in DNA.")
    elif args.r == 1 and args.alphabet == 'DNA':
        args.r = True

    # Set alphabet parameter.
    if args.alphabet == 'DNA':
        args.alphabet = index_list.DNA
    elif args.alphabet == 'RNA':
        args.alphabet = index_list.RNA
    elif args.alphabet == 'PROTEIN':
        args.alphabet = index_list.PROTEIN

    res = make_kmer_vector(k=args.k, alphabet=args.alphabet, filename=args.inputfile, revcomp=args.r)

    # Write correspond res file.
    if args.f == 'svm':
        from util import write_libsvm
        write_libsvm(res, [0] * len(res), args.outputfile)
    elif args.f == 'tab':
        from util import write_tab
        write_tab(res, args.outputfile)

    print("Done.")