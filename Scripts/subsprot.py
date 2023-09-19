import sys, getopt, math

def main(argv):

    alignment_file = ''
    ref = ''

    arg_list = sys.argv[1:]
    short_options = 'a:r:h:'
    long_options = ['alignment_file=', 'ref=', 'help']
    usage = 'subsprot.py\n\nExample:\n subsprot.py -a input.fasta -r "ref name"\
 \n\nUsage:\n -h, --help \n -a, --alignment_file, multi sequence alignment file in .fasta format [required]\
 \n -r, --ref, reference strain name [first sequence] \n\
 \n'

    try:
        options, args = getopt.getopt(arg_list,short_options,long_options)
    except getopt.GetoptError:
        print usage
        sys.exit(2)

    if len(arg_list)<1 or not options:
        print usage
        sys.exit(2)

    print '\n', 'ID subs in protein alignment', '\n'

    for arg, val in options:
        if arg in ('-a', '--alignment_file'):
            alignment_file = val
        elif arg in ('-r', '--ref'):
            ref = ">" + val
        elif arg in ('-h', '--help'):
            print 'This program identifies point substitutions in protein alignment files.\n'
            print usage
            sys.exit()
        else:
            print usage
            sys.exit()

    # Open fasta file
    with open(alignment_file,'r') as f1:
        f_handler = f1.read().splitlines()

    # assign sequences and strains lists
    seq = f_handler[1::2]
    strains = f_handler[0::2]
    print("ref is" + str(ref))

    # find reference sequence in strains list
    if ref not in strains:
        ref = strains[0]

    print("ref is" + str(ref))
    ref_index = strains.index(ref)
    ref_seq = seq[ref_index]

    final_list = [["strain","RefPosMut"]]
    tmp_list = list()
    seq_length = len(ref_seq)
    c=0
    for i in range(len(seq)):
        tmp = list()
        if seq[i] != ref_seq:
            tmp.append(strains[i][1:])
            for j in range(seq_length):
                if seq[i][j] != "X" and seq[i][j] != ref_seq[j]:
                    tmp.append(ref_seq[j]+str(j+2)+seq[i][j])
            if tmp == [strains[i][1:]]:
                tmp = list()
            else:
                final_list.append(tmp)

    for item in final_list:
        with open(str(item[0])+'.txt','w') as f2:
            for line in item:
                print >>f2, line
            f2.close


        #f2.write("%s\n" % item[1:])

#print('\n'.join(' '.join(map(str,sl)) for sl in final_list))
if __name__ == "__main__":
    main(sys.argv[1:])
