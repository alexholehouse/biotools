#!/usr/bin/python
# translate_codes
#
# Alex Holehouse
# Febuary 2013
# v 0.1
#
# Simple script which takes a FASTA style sequence and converts it into a CAMPARI-ready sequence file. Note this DOES NOT 
# switch out His for CAMPARI readable histadine 3 amino acid codes (Hie and Hid)
# 



# ---------------------------------------
# Public function which converts a single letter
# amino acid code to a three letter code. Returns
# none on failure
def convert_AA_1_to_3(one_AACode):
    """
        Convert a one aminod acid letter code to a three
        amino acid letter code. Excpects input as a string
        or something that can be coerced into a string.

        Returns a three amino acid letter code as a string
        on success, or None on failure.

        one_AAcode         Single amino acid code.

    """
    try:
        return (AA1_TO_3[str(i).upper()])
    except KeyError:
        return None

def convert_AA_sequence_1_to_3(sequence):
    returnSequence = []

    for i in sequence:
        converted = convert_AA_1_to_3(i)
        if converted:
            returnSequence.append(converted)

    return returnSequence


from ensemble.biotools.protein.aaGroups import AA1_TO_3

if __name__ == "__main__":
    from Bio import SeqIO
    import string
    import argparse	
    import sys	

    parser = argparse.ArgumentParser(description='Convert a FASTA sequence to a CAMPARI ready .in sequence file. [Contact alex@holehouse.org]')
    parser.add_argument('filename', metavar='filename',  nargs='+',
                        help='.fasta file')

    args = parser.parse_args()
    inputString = args.filename[0]
    sequence = {}
    
    ## READ IN FILE
    handle = open(inputString, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        sequence[record.id] = record.seq
    handle.close()
    
    # CONVERT EACH AA from file and 
    for entry in sequence:
        f = open(str(entry)+".in", 'w')
        f.write("ACE\n")
        for i in sequence[entry]:
            
            converted = convert_AA_1_to_3(i)
            if converted:
                f.write(converted)
                f.write("\n")

            # to avoid silently missing certain residues (e.g. X)
            elif not i == "\n":
                print "WARNING: Tried to convert character '" + str(i) + "' but failed. Skipping..."


        f.write("NME\nEND\n")
        f.close()
