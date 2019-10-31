"""
Python 3 script containing functions for writing a linear alignment to a .SAM file
"""

import sys, os

def create_sam_header(filepath, RNAME, RLEN, version=1.4, sorting_order="unsorted"):

    # creating a new SAM file
    sam_file = open(filepath, 'w')

    # adding headers
    sam_file.write("@HD\tVN:%.1f\tSO:%s\n" % (version, sorting_order))
    sam_file.write("@SQ\tSN:%s\tLN:%d\n" % (RNAME, RLEN))
    sam_file.close()


def append_sam_alignment(filepath, reverse=False, secondary=False, QNAME="*", RNAME="*", POS=0, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*"):

    # calculating the flag
    FLAG = 0
    if reverse: FLAG += 16
    if secondary: FLAG += 256

    # open existing SAM file append one linear alignment as a line to the alignment body
    sam_file = open(filepath, "a")
    sam_file.write("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL))
    sam_file.close()

def compare_sam__files(filepath1, filepath2, fields):

    # Check if file path exist
    if not os.path.isfile(filepath1):
        print("File path {} does not exist.".format(filepath1))
        sys.exit()
    if not os.path.isfile(filepath2):
        print("File path {} does not exist.".format(filepath2))
        sys.exit()

    # open files
    file1 = open(filepath1, "r")
    file2 = open(filepath2, "r")
    identical = True

    # remove all lines starting with @
    lines1 = list(enumerate(file1))
    i = 0
    while i < len(lines1):
        if(lines1[i][1][0]) == "@":
            del lines1[i]
        else:
            i = i+1

    lines2 = list(enumerate(file2))
    i = 0
    while i < len(lines2):
        if (lines2[i][1][0]) == "@":
            del lines2[i]
        else:
            i = i + 1

    file1.close()
    file2.close()

    # Exit if remaining lines are not equal in number
    if len(lines1) != len(lines2):
        print("Alignment parts of files do not have equal lengths !")
        sys.exit()

    # Compare relevant elements of alignments
    for i in range(0, len(lines1)):
        split1 = lines1[i][1].split("\t")
        split2 = lines2[i][1].split("\t")
        for f in fields:
            if split1[f-1] != split2[f-1]:
                identical = False
                print("Field %d not identical in line %d of file 1 and line %d of file 2" % (f, lines1[i][0]+1, lines2[i][0]+1))

    if identical:
        print("The files are identical")


def main():
    filepath1 = "test1.sam"
    filepath2 = "test2.sam"
    ref_name = "ref"
    create_sam_header(filepath1, RNAME=ref_name, RLEN=100)
    append_sam_alignment(filepath1, reverse=True, secondary=True, QNAME="read11", RNAME=ref_name, POS=2, CIGAR="12M")
    append_sam_alignment(filepath1, reverse=False, secondary=True, QNAME="read12", RNAME=ref_name, POS=2, CIGAR="12M")
    append_sam_alignment(filepath1, reverse=True, secondary=False, QNAME="read13", RNAME=ref_name, POS=2, CIGAR="12M")
    create_sam_header(filepath2, RNAME=ref_name, RLEN=100)
    append_sam_alignment(filepath2, reverse=True, secondary=True, QNAME="read21", RNAME=ref_name, POS=2, CIGAR="12M")
    append_sam_alignment(filepath2, reverse=False, secondary=True, QNAME="read22", RNAME=ref_name, POS=2, CIGAR="12M")
    append_sam_alignment(filepath2, reverse=True, secondary=False, QNAME="read23", RNAME=ref_name, POS=2, CIGAR="12M")

    fields = [1, 2, 3, 4, 6]
    compare_sam__files("test1.sam", "test2.sam", fields)

if __name__ == '__main__':
    main()