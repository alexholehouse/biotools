#!/usr/bin/python
# klean
#
# Alex Holehouse
# January 2013
# v 0.1


if __name__ == "__main__":
    import string
    import argparse	
    import sys						
    from subprocess import Popen, PIPE

    # test that xsel is installed. 
    try:
        p = Popen(['xsel','-pi'], stdin=PIPE)
    except:
        print " =============================== ERROR! ===============================\n\n          klean requires an external program called xsel.\n\nPlease ensure it is installed and configured to be called from the command\nline correctly. In Ubuntu this is as simple as 'sudo apt-get install xsel'.\nFor more information see http://www.vergenet.net/~conrad/software/xsel/.\n\nTo check it has been installed correctly, type 'xsel --help' in the terminal.\nYou should be greeted with xsel's help information.\n\nFor further help contact alex@holehouse.org"
        exit(1)


    parser = argparse.ArgumentParser(description='Clear an input string of spaces, newlines, tabs etc, and copy to clipboard. [Contact alex@holehouse.org]')
    parser.add_argument('input', metavar='input',  nargs='+',
                        help='input string (in quotation marks)')

    args = parser.parse_args()

    # get input string
    inputString = args.input[0]

    # get characters we care about 
    charSet = 'acgtuACGTU'

    # build and add newstring to the middle mouse click
    p.communicate(input=("".join([x for x in inputString if x in charSet])))
