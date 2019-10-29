"""
Python 3 script containing functions for writing a linear alignment to a .SAM file
"""

def create_sam_header(filename, RNAME, RLEN, version=1.4, sorting_order="unsorted"):

    # creating a new SAM file
    sam_file = open(filename, 'w')

    # adding headers
    sam_file.write("@HD\tVN:%.1f\tSO:%s\n" % (version, sorting_order))
    sam_file.write("@SQ\tSN:%s\tLN:%d\n" % (RNAME, RLEN))
    sam_file.close()


def append_sam_alignment(filename, reverse=False, secondary=False, QNAME="*", RNAME="*", POS=0, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*"):

    # calculating the flag
    FLAG = 0
    if reverse: FLAG += 16
    if secondary: FLAG += 256

    # open existing SAM file append one linear alignment as a line to the alignment body
    sam_file = open(filename, "a")
    sam_file.write("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL))
    sam_file.close()


def main():
    filename = "test.sam"
    ref_name = "ref"
    create_sam_header(filename, RNAME=ref_name, RLEN=100)
    append_sam_alignment(filename, reverse=True, secondary=True, QNAME="read1", RNAME=ref_name, POS=2, CIGAR="12M")


if __name__ == '__main__':
    main()