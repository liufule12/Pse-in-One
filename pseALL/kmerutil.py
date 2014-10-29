__author__ = 'Fule Liu'

import sys
import re

from pseALL.util import frequency
from pseALL.util import get_data
from pseALL.data import index_list

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
            print("Unknown DNA character (%s)\n" % letter)
            sys.exit(0)

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
            print("Error, Only DNA sequence can be reverse compliment.")
            sys.exit(0)

        kmer_list = make_kmer_list(k, alphabet)

        count_sum = 0
        vector = []
        for seq in seq_list:
            # Generate the kmer frequency dict.
            kmer_count = {}
            for kmer in kmer_list:
                temp_count = frequency(seq, kmer)
                if revcomp:
                    rev_kmer = find_revcomp(kmer, {})
                    if kmer <= rev_kmer:
                        if kmer not in kmer_count:
                            kmer_count[kmer] = 0
                        kmer_count[kmer] += temp_count
                    else:
                        if rev_kmer not in kmer_count:
                            kmer_count[rev_kmer] = 0
                        kmer_count[rev_kmer] += temp_count
                else:
                    if kmer not in kmer_count:
                        kmer_count[kmer] = 0
                    kmer_count[kmer] += temp_count
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
    read_file = "data/test.fasta"

    res = make_kmer_vector(2, index_list.RNA, read_file)
    print(len(res[0]), res)

    res = make_kmer_vector(2, index_list.PROTEIN, read_file)
    print(len(res[0]), res)