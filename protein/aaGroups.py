############################################################################
# aaGroups.py
#
# Alex Holehouse (alex@holehouse.org)
# Written January 2013
#
# A note on using these values
#
# Included are the relevant references for values tabulated
# in this file. If you use any of these methods PLEASE cite
# the source data references appropriatly.
#
# CHANGELOG
# ~ Feb 22 2013
# Added AA1_TO_3 group 
#

###############################################


# Amino acids!
AA = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
AA1_TO_3 ={"A":"Ala","C":"Cys","D":"Asp","E":"Glu", "F":"Phe", "G":"Gly", "H":"His","I":"Ile","K":"Lys","L":"Leu","M":"Met","N":"Asn","P":"Pro","Q":"Gln","R":"Arg","S":"Ser","T":"Thr","V":"Val","W":"Trp","Y":"Tyr"}


# ------------------------------------------------------------
# Grouping (based on standard groups)
#





NEG = ("D","E")
POS = ("K", "R", "H")
HYDROPHOBIC = ("A", "V", "L", "I", "M", "F", "W", "P","Y")
POLAR = ("S","T","N","Q","H","C") 
BULKYHYDROPHOBES = ("L","V","I", "F", "Y", "W", "P") # somewhat subjective...
NONPOLAR = ("A", "V", "L", "I", "M", "F", "W", "P","Y", "G")


# side chain chemistry based grouping
AROMATIC = ("F","W", "Y")
AMINES = ("Q", "N")
HYDROXYLS = ("Y","T", "S")


## ADDING ADDITIONAL GROUPS
# You (the user) can define additional sets below, e.g.
# BCAA = ("V","I","L")




# ------------------------------------------------------------
# ENTROPY VALUES FROM 
# " Role of Main-chain Electrostatics, Hydrophobic Effect
# and Side-chain Conformational Entropy in
# Determining the Secondary Structure of Proteins "
#
# F. Avbelj and L. Fele
#
# J. Mol. Biol. (1998) 279, 665-684
#

COIL_ENTROPY = {'A':0,\
                    'C':-0.572,\
                    'D':-1.318,\
                    'E':-1.763,\
                    'F':-0.554,\
                    'G':0,\
                    'H':-0.895,\
                    'I':-0.926,\
                    'K':-1.873,\
                    'L':-0.763,\
                    'M':-1.549,\
                    'N':-1.318,\
                    'P':0,\
                    'Q':-1.763,\
                    'R':-2.120,\
                    'S':-1.695,\
                    'T':-1.618,\
                    'V':-0.541,\
                    'W':-0.909,\
                    'Y':-1.019}

HELIX_ENTROPY = {'A':0,\
                    'C':-0.535,\
                    'D':-0.959,\
                    'E':-1.547,\
                    'F':-0.409,\
                    'G':0,\
                    'H':-0.794,\
                    'I':-0.481,\
                    'K':-1.849,\
                    'L':-0.696,\
                    'M':-1.452,\
                    'N':-1.436,\
                    'P':0,\
                    'Q':-1.547,\
                    'R':-1.991,\
                    'S':-1.686,\
                    'T':-1.363,\
                    'V':-0.172,\
                    'W':-0.633,\
                    'Y':-0.858}

# ------------------------------------------------------------
# HYDROPHOBICITY VALUES FROM 
# Kyte-Doolitle hydrophobicity scale from
#
# "A simple method for displaying the hydrophobic character
#  of a protein"
#
# Jack Kyte and Russell Doolitle
#
# J. Mol. Biol. (1982) 157 105-132
# 
# May add additional scales

KD_HYDROPHOBICITY = {'A':1.8,\
                    'C':2.5,\
                    'D':-3.5,\
                    'E':-3.5,\
                    'F':2.8,\
                    'G':-0.4,\
                    'H':-3.2,\
                    'I':4.5,\
                    'K':-3.9,\
                    'L':3.8,\
                    'M':1.9,\
                    'N':-3.5,\
                    'P':-1.6,\
                    'Q':-3.5,\
                    'R':-4.5,\
                    'S':-0.8,\
                    'T':-1.3,\
                    'V':4.2,\
                    'W':-0.9,\
                    'Y':-1.3}


# ------------------------------------------------------------
# Alpha helix propensity values from 
# "A helix propensity scale based on experimental studies of 
#  peptides and proteins."
#
# C N Pace and J M Scholtz
#
# Biophys J. 1998 July; 75(1): 422-427
#
# Note: Values in kcal mol-1
#
HELIX_PROPENSITY_KCAL = {'A':0,\
                    'C':0.68,\
                    'D':0.69,\
                    'E':0.4,\
                    'F':0.54,\
                    'G':1,\
                    'H':0.66,\
                    'I':0.41,\
                    'K':0.26,\
                    'L':0.21,\
                    'M':0.24,\
                    'N':0.65,\
                    'P':3.16,\
                    'Q':0.39,\
                    'R':0.21,\
                    'S':0.5,\
                    'T':0.49,\
                    'V':0.61,\
                    'W':0.49,\
                    'Y':0.53}

# ------------------------------------------------------------
# Beta sheet propensity values from 
# "Measurement of the beta-sheet-forming pronensities of amino acids"
#
# Minor DL Jr, Kim PS.
#
# Nature. 1994 Feb 17;367(6464):660-3.
#
# Note: Values in kcal mol-1
#       The value for proline here seems somewhat suspect.
#xxxxx
SHEET_PROPENSITY_KCAL = {'A':0,\
                    'C':0.52,\
                    'D':-0.94,\
                    'E':0.01,\
                    'F':0.86,\
                    'G':-1.2,\
                    'H':-0.02,\
                    'I':1.0,\
                    'K':0.27,\
                    'L':0.51,\
                    'M':0.72,\
                    'N':-0.08,\
                    'P':-3,\
                    'Q':0.23,\
                    'R':0.45,\
                    'S':0.7,\
                    'T':1.1,\
                    'V':0.82,\
                    'W':0.54,\
                    'Y':0.96}

# ------------------------------------------------------------
# Alpha normalized helix and beta sheet propensity values from 
# "Dependence of alpha-helical and beta-sheet amino acid propensities \
#  on the overall protein fold type"
#
# Kazuo Fujiwara*, Hiromi Toda and Masamichi Ikeguchi
#
# BMC Structural Biology 2012, 12:18 doi:10.1186/1472-6807-12-18
#
HELIX_PROPENSITY = {'A':1.41,\
                    'C':0.85,\
                    'D':0.82,\
                    'E':1.39,\
                    'F':1.00,\
                    'G':0.44,\
                    'H':0.87,\
                    'I':1.04,\
                    'K':1.17,\
                    'L':1.28,\
                    'M':1.26,\
                    'N':0.73,\
                    'P':0.44,\
                    'Q':1.26,\
                    'R':1.21,\
                    'S':0.76,\
                    'T':0.78,\
                    'V':0.91,\
                    'W':1.07,\
                    'Y':0.98}

SHEET_PROPENSITY = {'A':0.75,\
                    'C':1.36,\
                    'D':0.55,\
                    'E':0.65,\
                    'F':1.4,\
                    'G':0.67,\
                    'H':0.99,\
                    'I':1.79,\
                    'K':0.76,\
                    'L':1.15,\
                    'M':1.01,\
                    'N':0.63,\
                    'P':0.40,\
                    'Q':0.72,\
                    'R':0.85,\
                    'S':0.81,\
                    'T':1.21,\
                    'V':2.00,\
                    'W':1.23,\
                    'Y':1.37}

# ------------------------------------------------------------
# SOLVENT EXPOSED VALUES FROM
# "Suggestions for "safe" residue substitutions in site-directed mutagenesis."
# J Mol Biol. 1991 Feb 20;217(4):721-9.
# Bordo D, Argos P.
#
LOW_EXPOSED = ("I", "C", "L", "V", "F","M", "W","A")
MED_EXPOSED = ("G", "H","P", "S", "T", "Y")
HIGH_EXPOSED = ("R", "D", "N", "E", "Q", "K")
