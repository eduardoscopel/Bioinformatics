import sys
import getopt

def main(argv):
    alignment_file = ''
    ref = ''

    try:
        # Define command-line options and arguments
        short_options = 'a:r:h:'
        long_options = ['alignment_file=', 'ref=', 'help']

        # Get command-line arguments and options
        options, args = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)

    if not options or '-h' in argv:
        print_usage()
        sys.exit(2)

    print('\nID subs in protein alignment\n')

    for arg, val in options:
        if arg in ('-a', '--alignment_file'):
            alignment_file = val
        elif arg in ('-r', '--ref'):
            ref = ">" + val

    # Open the fasta file
    with open(alignment_file, 'r') as f1:
        f_handler = f1.read().splitlines()

    # Extract sequences and strains from the file
    seq = f_handler[1::2]
    strains = f_handler[0::2]

    # If no reference is specified, use the first sequence as the reference
    if not ref:
        ref = strains[0]

    ref_index = strains.index(ref)
    ref_seq = seq[ref_index]

    final_list = [["strain", "RefPosMut"]]

    seq_length = len(ref_seq)

    for i in range(len(seq)):
        tmp = []
        if seq[i] != ref_seq:
            tmp.append(strains[i][1:])
            # Compare each position in the sequence
            for j in range(seq_length):
                if seq[i][j] != "X" and seq[i][j] != ref_seq[j]:
                    tmp.append(ref_seq[j] + str(j + 2) + seq[i][j])
            if len(tmp) > 1:
                final_list.append(tmp)

    for item in final_list:
        # Write results to separate files
        with open(str(item[0]) + '.txt', 'w') as f2:
            for line in item:
                f2.write(line + '\n')

def print_usage():
    print('This program identifies point substitutions in protein alignment files.\n')
    print('Usage:')
    print('-h, --help')
    print('-a, --alignment_file, multi-sequence alignment file in .fasta format [required]')
    print('-r, --ref, reference strain name [first sequence]')

if __name__ == "__main__":
    main(sys.argv[1:])
