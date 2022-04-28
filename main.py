from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys, getopt


# function to calculate gc precentage in DNA
def gc(DNA):
    """This command takes a seq as a string and returns the gc percentage of it."""
    for i in DNA:
        if i not in 'AGCTN':
            return 'Invalid Seq'
    DNA = DNA.upper()
    nBases = DNA.count('N')
    gcBases = DNA.count('G') + DNA.count('C')
    precantage = (gcBases / (len(DNA) - nBases)) * 100
    return precantage


# function to find complement of DNA
def complement(DNA):
    DNA = DNA.upper()
    baseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    letters = list(DNA)
    letters = [baseComplement[base] for base in letters]
    return ''.join(letters)


# function to find reverse of DNA
def reverse(DNA):
    return DNA[::-1]


# funciton to find reverse complement of DNA
def reverse_complement(DNA):
    """This command takes a seq as a string and returns its reverse complement."""
    for i in DNA:
        if i not in 'AGCT':
            return 'Invalid Seq'
    seq = complement(DNA)
    seq = reverse(seq)
    return seq


def transcribe(DNA):
    """This command takes a seq as a string and returns its transcription."""
    DNA = DNA.upper()
    for i in DNA:
        if i not in 'AGCT':
            return 'Invalid Seq'
    return DNA.replace('T', 'U')


def calc_nbases(DNA):
    """This command takes a seq and calculates its nbases."""
    DNA = DNA.upper()
    for i in DNA:
        if i not in 'AGCTN':
            return 'Invalid Seq'
    return DNA.count('N')


# function to check if there amino acid in protein sequence
def check_Protein(protein):
    protein = protein.upper()
    for i in protein:
        if i not in 'ABCDEFGHIKLMNPQRSTVWXYZ':
            return False
    return True


def check_DNA(DNA):
    DNA = DNA.upper()
    for i in DNA:
        if i not in 'AGCT':
            return False
    return True


def check_RNA(RNA):
    RNA = RNA.upper()
    for i in RNA:
        if i not in 'AGCU':
            return False
    return True


def is_valid(Seq, type):
    """This command takes a seq and a type (protein, dna, rna) and returns a Boolean value of whether it’s a valid
    type or not """
    type = type.lower()
    if type == 'protein':
        if check_Protein(Seq):
            return "it's a valid protein"
        else:
            return "it's not valid protein"
    elif type == 'dna':
        if check_DNA(Seq):
            return "it's a valid Dna"
        else:
            return "it's not valid Dna"
    elif type == 'rna':
        if check_RNA(Seq):
            return "it's a valid Rna"
        else:
            return "it's not valid Rna"
    else:
        return 'Invalid Type or invalid Sequance'


def filter_nbases(Seq):
    """This command takes a seq and returns the Seq after removing n bases."""
    Seq = Seq.upper()
    for i in Seq:
        if i not in 'AGCTN':
            return 'Invalid Seq'
    Seq = Seq.replace("N", "")
    return Seq


def output_alignment(alignments, output):
    f = open(output, 'w')
    for alignment in alignments:
        nonFormattedAlignment = str(alignment)
        f.write(nonFormattedAlignment)
        f.write('\n')
        formattedAlignment = str(format_alignment(*alignment))
        f.write(formattedAlignment)
        f.write('\n')
    f.close()


def seq_alignment(seq1, seq2, output=""):
    """This command takes 2 sequences as input and get all their alignments along with the score. The -o is an
    optional parameter if we need the output to be written on a file instead of the screen. """
    for i in seq1:
        if i not in 'AGCT':
            print('Seq1 Invalid')
            return
    for j in seq2:
        if j not in 'AGCT':
            print('Seq2 Invalid')
            return
    alignments = pairwise2.align.globalxx(seq1, seq2)  # global alignment
    if output == '':
        for alignment in alignments:
            print(alignment)
            print(format_alignment(*alignment))
    else:
        output_alignment(alignments, output)
        print('Alignmnet Done to File ', output)


def seq_alignment_files(file1, file2, outputfile=""):
    """This command takes 2 fasta files as input, each file contains a single sequence. It reads the 2 sequences from
    files and get all their alignments along with the score. The -o is an optional parameter if we need the output to
    be written on a file instead of the screen. """
    try:
        seq1 = SeqIO.read(file1, 'fasta')
        seq2 = SeqIO.read(file2, 'fasta')
    except OSError as Error:
        print(Error)
        return 'Please Enter a valid File name'
    alignments = pairwise2.align.globalxx(seq1, seq2)  # global alignment
    if outputfile == '':
        for alignment in alignments:
            print(alignment)
            print(format_alignment(*alignment))
    else:
        output_alignment(alignments, outputfile)
        print('Alignmnet Done to File ', outputfile)


def write_blast(blast_record, output):
    f = open(output, 'w')
    f.write("*******************Description*******************\n")
    for description in blast_record.descriptions:
        f.write("Title " + str(description.title) + '\n')
        f.write("Score " + str(description.score) + '\n')
        f.write("e " + str(description.e) + '\n')
        f.write("Number of Alignments " + str(description.num_alignments) + '\n')

    f.write("*******************Alignments*******************\n")
    for alignment in blast_record.alignments:
        f.write("Title " + str(alignment.title) + '\n')
        f.write("Length " + str(alignment.length) + '\n')
        f.write("*****hsp*****\n")
        for hsp in alignment.hsps:
            f.write("score " + str(hsp.score) + '\n')
            f.write("bits " + str(hsp.bits) + '\n')
            f.write("expect " + str(hsp.expect) + '\n')
            f.write("alignments " + str(hsp.num_alignments) + '\n')
            f.write("identities " + str(hsp.identities) + '\n')
            f.write("positives " + str(hsp.positives) + '\n')
            f.write("gaps " + str(hsp.gaps) + '\n')
            f.write("strand " + str(hsp.strand) + '\n')
            f.write("frame " + str(hsp.frame) + '\n')
            f.write("query " + str(hsp.query) + '\n')
            f.write("query_start " + str(hsp.query_start) + '\n')
            f.write("match " + str(hsp.match) + '\n')
            f.write("sbjct " + str(hsp.sbjct) + '\n')
            f.write("sbjct_start " + str(hsp.sbjct_start) + '\n')

    f.write("multiple_alignment " + str(blast_record.multiple_alignment) + '\n')


def online_alignment(seq, output=''):
    """This command takes a sequence and uses BLAST to search the internet for its alignments. The output should be
    all the information in the resultant BLAST record. The -o is an optional parameter if we need the output to be
    written on a file instead of the screen. """
    try:
        print("Working...")
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        if output == '':
            print("*******************Description*******************")
            for description in blast_record.descriptions:
                print("Title ", description.title)
                print("Score ", description.score)
                print("e ", description.e)
                print("Number of Alignments ", description.num_alignments)

            print("*******************Alignments*******************")
            for alignment in blast_record.alignments:
                print("Title ", alignment.title)
                print("Length ", alignment.length)
                print("*****hsp*****")
                for hsp in alignment.hsps:
                    print("score ", hsp.score)
                    print("bits ", hsp.bits)
                    print("expect ", hsp.expect)
                    print("alignments ", hsp.num_alignments)
                    print("identities ", hsp.identities)
                    print("positives ", hsp.positives)
                    print("gaps ", hsp.gaps)
                    print("strand ", hsp.strand)
                    print("frame ", hsp.frame)
                    print("query ", hsp.query)
                    print("query_start ", hsp.query_start)
                    print("match ", hsp.match)
                    print("sbjct ", hsp.sbjct)
                    print("sbjct_start ", hsp.sbjct_start)

            print("multiple_alignment ", blast_record.multiple_alignment)
        else:
            write_blast(blast_record, output)
    except OSError as error:
        print("Can not Connect ")
        print(error)
        return
    print("Done.")


def merge_fasta(file1, file2, *files, output=''):
    """This command takes any number of fasta files (at least two) and merge their contents into one fasta output
    file. There is an option to write the merge result in a file using -o option, otherwise the merge result will be
    displayed on the console. """

    try:
        file1 = SeqIO.parse(file1, 'fasta')
        file2 = SeqIO.parse(file2, 'fasta')
    except OSError as Error:
        print(Error)
        return
    FilesList = []
    FilesList.append(file1)
    FilesList.append(file2)
    for file in files:
        try:
            print(file)
            file = SeqIO.parse(file, 'fasta')
        except OSError as Error:
            print(Error)
            return
        FilesList.append(file)

    if output == '':
        for file in FilesList:
            for record in file:
                print(record.id)
                print(record.description)
                print(record.seq)
    else:
        with open('output.fasta', 'w')as outputFile:
            for file in FilesList:
                SeqIO.write(file, outputFile, 'fasta')


def convert_to_fasta(file):
    """This command converts the input genbank file with multiple records onto a fasta formatted file. The output is
    to be written in a different output fasta file. """
    if file[-3:] == 'gbk':
        output = file[:-3] + 'fasta'
        try:
            with open(file)as input:
                sequance = SeqIO.parse(input, 'genbank')
                SeqIO.write(sequance, output, 'fasta')
                print('Conversion Done to %s' % output)
        except OSError as Error:
            print(Error)
            return
    else:
        print('File must be genbank\n', convert_to_fasta.__doc__)


def usage():
    print('gc: \n', gc.__doc__, '\n')
    print('transcribe:\n', transcribe.__doc__, '\n')
    print('reverse_complement:\n', reverse_complement.__doc__, '\n')
    print('calc_nbases:\n', calc_nbases.__doc__, '\n')
    print('is_valid:\n', is_valid.__doc__, '\n')
    print('filter_nbases:\n', filter_nbases.__doc__, '\n')
    print('seq_alignment:\n', seq_alignment.__doc__, '\n')
    print('seq_alignment_files:\n', seq_alignment_files.__doc__, '\n')
    print('online_alignment:\n', online_alignment.__doc__, '\n')
    print('merge_fasta:\n', merge_fasta.__doc__, '\n')
    print('convert_to_fasta:\n', convert_to_fasta.__doc__, '\n')


try:
    oplist, args = getopt.gnu_getopt(sys.argv[1:], 'o:h', 'help')
except getopt.GetoptError as Error:
    print(Error)
    sys.exit()
if args != []:
    if args[0] == 'gc':
        if not oplist:
            if len(args) == 2:
                print(gc(args[1]))
            else:
                print('Less or more number of prammeters:\n', gc.__doc__)
        else:
            print('This command dose not have option %s\n' % oplist[0][0], gc.__doc__)
    elif args[0] == 'transcribe':
        if not oplist:
            if len(args) == 2:
                print(transcribe(args[1]))
            else:
                print('Less or more number of prammeters:\n', transcribe.__doc__)
        else:
            print('This command dose not have option %s\n' % oplist[0][0], transcribe.__doc__)
    elif args[0] == 'reverse_complement':
        if not oplist:
            if len(args) == 2:
                print(reverse_complement(args[1]))
            else:
                print('Less or more number of prammeters:\n', reverse_complement.__doc__)
        else:
            print('This command dose not have option %s\n' % oplist[0][0], reverse_complement.__doc__)
    elif args[0] == 'calc_nbases':
        if not oplist:
            if len(args) == 2:
                print(calc_nbases(args[1]))
            else:
                print('Less or more number of prammeters:\n', calc_nbases.__doc__)
        else:
            print('This command dose not have option %s\n' % oplist[0][0], calc_nbases.__doc__)
    elif args[0] == 'is_valid':
        if not oplist:
            if len(args) == 3:
                print(is_valid(args[1], args[2]))
            else:
                print('Less or more number of prammeters:\n', is_valid.__doc__)
        else:
            print('This command dose not have option %s\n' % oplist[0][0], is_valid.__doc__)

    elif args[0] == 'filter_nbases':
        if not oplist:
            if len(args) == 2:
                print(filter_nbases(args[1]))
            else:
                print('Less or more number of prammeters:\n', filter_nbases.__doc__)
        print('This command dose not have option %s\n' % oplist[0][0], filter_nbases.__doc__)

    elif args[0] == 'seq_alignment':
        if not oplist:
            if len(args) == 3:
                seq_alignment(args[1], args[2])
            else:
                print('Less or more number of prammeters:\n', seq_alignment.__doc__)
        else:
            if '-o' in oplist[0]:
                if len(args) == 3:
                    seq_alignment(args[1], args[2], oplist[0][1])
                else:
                    print('Less or more number of prammeters:\n', seq_alignment.__doc__)
            else:
                print('Wrong option %s not recognized' % oplist[0][0])
    elif args[0] == 'seq_alignment_files':
        if not oplist:
            if len(args) == 3:
                seq_alignment_files(args[1], args[2])
            else:
                print('Less or more number of prammeters:\n', seq_alignment_files.__doc__)
        else:
            if '-o' in oplist[0]:
                if len(args) == 3:
                    seq_alignment_files(args[1], args[2], oplist[0][1])
                else:
                    print('Less or more number of prammeters:\n', seq_alignment_files.__doc__)
            else:
                print('Wrong option %s not recognized' % oplist[0])
    elif args[0] == 'online_alignment':
        if not oplist:
            if len(args) == 2:
                online_alignment(args[1])
            else:
                print('Less or more number of prammeters:\n', online_alignment.__doc__)
        else:
            if '-o' in oplist[0]:
                if len(args) == 2:
                    online_alignment(args[1], oplist[0][1])
                    print('Alignmnet Done to File ', oplist[0][1])
                else:
                    print('Less or more number of prammeters:\n', online_alignment.__doc__)
            else:
                print('Wrong option %s not recognized' % oplist[0])
    elif args[0] == 'merge_fasta':
        if not oplist:
            if len(args) == 3:
                merge_fasta(args[1], args[2])
            elif len(args) > 3:
                files = []
                for i in range(3, len(args)):
                    files.append(args[i])
                files = tuple(files)
                merge_fasta(args[1], args[2], files)
            else:
                print('less argumnents for merge_fasta\n', merge_fasta.__doc__)
        else:
            if len(args) == 3 and '-o' in oplist[0]:
                merge_fasta(args[1], args[2], oplist[0][1])
            elif len(args) > 3 and '-o' in oplist[0]:
                files = []
                for i in range(3, len(args)):
                    files.append(args[i])
                files = tuple(files)
                merge_fasta(args[1], args[2], files, oplist[0][1])
            else:
                print('less argumnents for merge_fasta or Wrong option\n', merge_fasta.__doc__)
    elif args[0] == 'convert_to_fasta':
        if not oplist:
            if len(args) == 2:
                convert_to_fasta(args[1])
            else:
                print('Less or more number of prammeters:\n', convert_to_fasta.__doc__)
        else:
            print('This command dose not have option %s\n' % oplist[0][0], convert_to_fasta.__doc__)

    else:
        print('Wrong Command please See usage')
        usage()

elif oplist:
    if '-h' in oplist[0] or '--help' in oplist[0]:
        usage()
