__author__ = 'Fule Liu'

import sys
import itertools


"""Used for process original data."""


class Seq:
    def __init__(self, name, seq, no):
        self.name = name
        self.seq = seq.upper()
        self.no = no
        self.length = len(seq)

    def __str__(self):
        """Output seq when 'print' method is called."""
        return "%s\tNo:%s\tlength:%s\n%s" % (self.name, str(self.no), str(self.length), self.seq)


def is_under_alphabet(s, alphabet):
    """Judge the string is within the scope of the alphabet or not.

    :param s: The string.
    :param alphabet: alphabet.

    Return True or the error character.
    """
    for e in s:
        if e not in alphabet:
            return e

    return True


def is_fasta(seq):
    """Judge the Seq object is in FASTA format.
    Two situation:
    1. No seq name.
    2. Seq name is illegal.
    3. No sequence.

    :param seq: Seq object.
    """
    if not seq.name:
        error_info = 'Error, sequence ' + str(seq.no) + ' has no sequence name.'
        print(seq)
        sys.stderr.write(error_info)
        return False
    if -1 != seq.name.find('>'):
        error_info = 'Error, sequence ' + str(seq.no) + ' name has > character.'
        sys.stderr.write(error_info)
        return False
    if 0 == seq.length:
        error_info = 'Error, sequence ' + str(seq.no) + ' is null.'
        sys.stderr.write(error_info)
        return False

    return True


def read_fasta(f):
    """Read a fasta file.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return Seq obj list.
    """
    name, seq = '', ''
    count = 0
    seq_list = []
    lines = f.readlines()
    for line in lines:
        if not line:
            break

        if '>' == line[0]:
            if 0 != count or (0 == count and seq != ''):
                if is_fasta(Seq(name, seq, count)):
                    seq_list.append(Seq(name, seq, count))
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    count += 1
    if is_fasta(Seq(name, seq, count)):
        seq_list.append(Seq(name, seq, count))
    else:
        sys.exit(0)

    return seq_list


def read_fasta_yield(f):
    """Yields a Seq object.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    """
    name, seq = '', ''
    count = 0
    while True:
        line = f.readline()
        if not line:
            break

        if '>' == line[0]:
            if 0 != count or (0 == count and seq != ''):
                if is_fasta(Seq(name, seq, count)):
                    yield Seq(name, seq, count)
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    if is_fasta(Seq(name, seq, count)):
        yield Seq(name, seq, count)
    else:
        sys.exit(0)


def read_fasta_check_dna(f, alphabet):
    """Read the fasta file, and check its legality.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return the seq list.
    """
    seq_list = []
    for e in read_fasta_yield(f):
        res = is_under_alphabet(e.seq, alphabet)
        if res:
            seq_list.append(e)
        else:
            error_info = 'Sorry, sequence ' + str(e.no) \
                         + ' has character ' + str(res) + '.(The character must be ' + alphabet + ').'
            sys.exit(error_info)

    return seq_list


def get_sequence_check_dna(f, alphabet):
    """Read the fasta file.

    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return the sequence list.
    """
    sequence_list = []
    for e in read_fasta_yield(f):
        res = is_under_alphabet(e.seq, alphabet)
        if res is not True:
            error_info = 'Error, sequence ' + str(e.no) \
                         + ' has character ' + str(res) + '.(The character must be ' + alphabet + ').'
            sys.exit(error_info)
        else:
            sequence_list.append(e.seq)

    return sequence_list


def is_sequence_list(sequence_list, alphabet):
    """Judge the sequence list is within the scope of alphabet and change the lowercase to capital."""
    count = 0
    new_sequence_list = []

    for e in sequence_list:
        e = e.upper()
        count += 1
        res = is_under_alphabet(e, alphabet)
        if res is not True:
            error_info = 'Sorry, sequence ' + str(count) \
                         + ' has illegal character ' + str(res) + '.(The character must be A, C, G or T)'
            sys.stderr.write(error_info)
            return False
        else:
            new_sequence_list.append(e)

    return new_sequence_list


def get_data(input_data, alphabet, desc=False):
    """Get sequence data from file or list with check.

    :param input_data: type file or list
    :param desc: with this option, the return value will be a Seq object list(it only works in file object).
    :return: sequence data or shutdown.
    """
    if hasattr(input_data, 'read'):
        if desc is False:
            return get_sequence_check_dna(input_data, alphabet)
        else:
            return read_fasta_check_dna(input_data, alphabet)
    elif isinstance(input_data, list):
        input_data = is_sequence_list(input_data, alphabet)
        if input_data is not False:
            return input_data
        else:
            sys.exit(0)
    else:
        error_info = 'Sorry, the parameter in get_data method must be list or file type.'
        sys.exit(error_info)


"""Some basic function for generate feature vector."""


def frequency(tol_str, tar_str):
    """Generate the frequency of tar_str in tol_str.

    :param tol_str: mother string.
    :param tar_str: substring.
    """
    i, j, tar_count = 0, 0, 0
    len_tol_str = len(tol_str)
    len_tar_str = len(tar_str)
    while i < len_tol_str and j < len_tar_str:
        if tol_str[i] == tar_str[j]:
            i += 1
            j += 1
            if j >= len_tar_str:
                tar_count += 1
                i = i - j + 1
                j = 0
        else:
            i = i - j + 1
            j = 0

    return tar_count


def write_libsvm(vector_list, label_list, write_file):
    """Write the vectors into disk in livSVM format."""
    len_vector_list = len(vector_list)
    len_label_list = len(label_list)
    if len_vector_list == 0:
        sys.exit("The vector is none.")
    if len_label_list == 0:
        sys.exit("The label is none.")
    if len_vector_list != len_label_list:
        sys.exit("The length of vector and label is different.")

    with open(write_file, 'w') as f:
        for ind1, vec in enumerate(vector_list):
            temp_write = str(label_list[ind1])
            for ind2, val in enumerate(vec):
                temp_write += ' ' + str(ind2+1) + ':' + str(vec[ind2])
            f.write(temp_write)
            f.write('\n')


def write_tab(_vecs, write_file):
    """Write the vectors into disk in tab format."""
    with open(write_file, 'w') as f:
        for vec in _vecs:
            f.write(str(vec[0]))
            for val in vec[1:]:
                f.write('\t' + str(val))
            f.write('\n')


def write_csv(_vecs, write_file):
    """Write the vectors into disk in csv format."""
    import csv
    with open(write_file, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for vec in _vecs:
            spamwriter.writerow(vec)


def convert_phyche_index_to_dict(phyche_index, alphabet):
    """Convert phyche index from list to dict."""
    # for e in phyche_index:
    #     print e
    len_index_value = len(phyche_index[0])
    k = 0
    for i in range(1, 10):
        if len_index_value < 4**i:
            sys.exit("Sorry, the number of each index value is must be 4^k.")
        if len_index_value == 4**i:
            k = i
            break
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


def make_kmer_list(k, alphabet):
    try:
        return ["".join(e) for e in itertools.product(alphabet, repeat=k)]
    except TypeError:
        print("TypeError: k must be an inter and larger than 0, alphabet must be a string.")
        raise TypeError
    except ValueError:
        print("TypeError: k must be an inter and larger than 0")
        raise ValueError


def standard_deviation(value_list):
    """Return standard deviation."""
    from math import sqrt
    from math import pow
    n = len(value_list)
    average_value = sum(value_list) * 1.0 / n
    return sqrt(sum([pow(e - average_value, 2) for e in value_list]) * 1.0 / (n - 1))


def normalize_index(phyche_index, alphabet, is_convert_dict=False):
    """Normalize the physicochemical index."""
    normalize_phyche_value = []
    for phyche_value in phyche_index:
        average_phyche_value = sum(phyche_value) * 1.0 / len(phyche_value)
        sd_phyche = standard_deviation(phyche_value)
        normalize_phyche_value.append([round((e - average_phyche_value) / sd_phyche, 2) for e in phyche_value])

    if is_convert_dict is True:
        return convert_phyche_index_to_dict(normalize_phyche_value, alphabet)

    print(normalize_phyche_value)
    return normalize_phyche_value


def read_k(alphabet, _method, k):
    import const
    if alphabet == 'Protein':
        return 1
    elif alphabet == 'RNA':
        return 2

    if _method in const.K_2_DNA_METHODS:
        return 2
    elif _method in const.K_3_DNA_METHODS:
        return 3
    elif _method == 'PseKNC':
        return k
    else:
        print("Error in read_k.")


def check_args(args, filename):
    """Check pse and acc method args."""
    import const
    if 'w' in args:
        if args.w < 0 or args.w > 1:
            print("Error: The value of w must be no less than 0 and no larger than 1.")
            return False
    if 'method' in args:
        if args.alphabet == 'DNA' and args.method not in const.METHODS_DNA:
            if filename == const.ACC_FILENAME:
                print("Error: the DNA method parameter can only be " + str(const.METHODS_DNA_ACC))
            if filename == const.PSE_FILENAME:
                print("Error: the DNA method parameter can only be " + str(const.METHODS_DNA_PSE))
            else:
                print("Error: the DNA method parameter error.")
            return False
        elif args.alphabet == 'RNA' and args.method not in const.METHODS_RNA:
            if filename == const.ACC_FILENAME:
                print("Error: the RNA method parameter can only be " + str(const.METHODS_RNA_ACC))
            if filename == const.PSE_FILENAME:
                print("Error: the RNA method parameter can only be " + str(const.METHODS_RNA_PSE))
            else:
                print("Error: the RNA method parameter error.")
            return False
        elif args.alphabet == 'Protein' and args.method not in const.METHODS_PROTEIN:
            if filename == const.ACC_FILENAME:
                print("Error: the protein method parameter can only be " + str(const.METHODS_PROTEIN_ACC))
            if filename == const.PSE_FILENAME:
                print("Error: the protein method parameter can only be " + str(const.METHODS_PROTEIN_PSE))
            else:
                print("Error: the protein method parameter error.")
            return False
    if 'k' in args:
        if args.k <= 0:
            print("Error: the value of k must be an inter and larger than 0.")
            return False
    return True


if __name__ == '__main__':
    # phyche_index = \
    #     [[0.026, 0.036, 0.031, 0.033, 0.016, 0.026, 0.014, 0.031, 0.025, 0.025, 0.026, 0.036, 0.017, 0.025, 0.016, 0.026],
    #      [0.038, 0.038, 0.037, 0.036, 0.025, 0.042, 0.026, 0.037, 0.038, 0.036, 0.042, 0.038, 0.018, 0.038, 0.025, 0.038],
    #      [0.020, 0.023, 0.019, 0.022, 0.017, 0.019, 0.016, 0.019, 0.020, 0.026, 0.019, 0.023, 0.016, 0.020, 0.017, 0.020],
    #      [1.69, 1.32, 1.46, 1.03, 1.07, 1.43, 1.08, 1.46, 1.32, 1.20, 1.43, 1.32, 0.72, 1.32, 1.07, 1.69],
    #      [2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
    #      [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]

    phyche_index = \
        [[0.026, 0.036, 0.031, 0.033, 0.016, 0.026, 0.014, 0.031, 0.025, 0.025, 0.026, 0.036, 0.017, 0.025, 0.016, 0.026],
         [0.038, 0.038, 0.037, 0.036, 0.025, 0.042, 0.026, 0.037, 0.038, 0.036, 0.042, 0.038, 0.018, 0.038, 0.025, 0.038]]

    for e in (normalize_index(phyche_index, alphabet='ACGT')):
        print(e)

    phyche_index_dict = normalize_index(phyche_index, alphabet='ACGT', is_convert_dict=True)
    print(phyche_index_dict)

    # print(normalize_index(phyche_index, True))
