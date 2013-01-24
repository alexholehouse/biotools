#!/usr/bin/python
# Alex Holehouse, Jan 15th 2013
# alex.holehouse@gmail.com
#
# Simple script or method which extracts the accession values or names from
# a SwissProt FASTA file
#
#

def __internal_handler(filename, request):
    
    if request == "acc":
        r = 1
    if request == "name":
        r = 2
    if request == "both":
        r = 3
    
   
    with open(args.filename[0]) as f:
        content = f.readlines()

    acc = []

    # search for sp| so we only search for swissprot records
    for line in content:
        if line.find("sp|") > -1 or line.find("tr|") > -1:

            if r == 3:
                temp = line.split("|")[2]
                nextVal = temp.split(" OS")[0]
                nextVal = line.split("|")[1] + " " + nextVal

            elif r == 2:
                temp = line.split("|")[r]
                nextVal = temp.split(" OS")[0]

            else:
                nextVal = line.split("|")[r]

            acc.append(nextVal)
            

    # use set to remove duplicates
    return acc


def get_acc_from_file(filename):
    """ Input:  Name of file as string
    
        Output: Returns a list of accessions ordered as they
                were in the file
    """

    return(__internal_handler(filename, "acc"))


def get_name_from_file(filename):
    """ Input:  Name of file as string
    
        Output: Returns a list of names ordered as they
                were in the file
    """

    return(__internal_handler(filename, "name"))


def get_acc_and_name_from_file(filename):
    """ Input:  Name of file as string
    
        Output: Returns a list of acessions and names ordered
                as they were in the file
    """

    return(__internal_handler(filename, "both"))



            

if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Extract accession values from a SwissProt FASTA file and print to STDOUT')
    parser.add_argument('filename', metavar='filename',  nargs='+',
                        help='FASTA file which contains the SwissProt sequences of interest ')


    parser.add_argument('--name', dest='nameR', action='store_const',
                        const=True, default=False,
                        help='get names')

    parser.add_argument('--acc', dest='accR', action='store_const',
                        const=True, default=False,
                        help='get accessions')

    parser.add_argument('--both', dest='bothR', action='store_const',
                        const=True, default=False,
                        help='get accessions and name together')

    parser.add_argument('--rows', dest='vert', action='store_const',
                        const=True, default=False,
                        help='print with one item per line')
    
    parser.add_argument('--list', dest='hoz', action='store_const',
                        const=True, default=False,
                        help='print as a comma seperated list')


    parser.add_argument('--num', dest='num', action='store_const',
                        const=True, default=False,
                        help='number of accessions in file')


    args = parser.parse_args()


    if args.num:
        out = get_acc_from_file(args.filename[0])
        print len(out)
        exit(0)


    if not args.nameR and not args.accR and not args.bothR:
        print "[ERROR] - Please select one from --name, --accession or --both"
        exit(1)

    if (args.nameR and args.accR) or (args.nameR and args.bothR) or (args.accR and args.bothR):
        print "[ERROR] - Please select one from --name, --accession or --both"
        exit(1)

    
    if args.nameR:
        out = get_name_from_file(args.filename[0])
    
    if args.accR:
        out = get_acc_from_file(args.filename[0])

    if args.bothR:
        out = get_acc_and_name_from_file(args.filename[0])

    ## Formatting ouput

    

    elif args.vert and args.hoz:
        print "[ERROR]: Please select one of --vertical or --horizontal"
        exit(1)
        
    if args.vert:
        for i in out:
            print i

    elif args.hoz:
        counter = 0
        size = len(out)
        for i in out:
            sys.stdout.write(i)
            counter = counter+1
            if counter < size: 
                sys.stdout.write(", ")
        print ""

    else:
        print out

    
            





    





