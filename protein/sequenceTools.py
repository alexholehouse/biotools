# A collection of functions for looking at amino acid sequences
# and their statistical properties.
#
# Note this file requires aaGroups.py for source data.
# PLEASE SEE THE DISCLAIMER IN aaGroups.py before
# using results from these tools.
#

from aaGroups import * # get all our classes 


# ---------------------------------------
#
class sequenceToolsException(Exception):
    """
       Exception class for sequenceTools methods
    """
    pass


# ---------------------------------------
# Nice decorator which santizizes input for all functions
def _convertToUpperCase_sanitize(fn):
    def wrapped(seq, A=None,B=None,C=None,D=None):

        # check it's a string
        try:
            UC_seq = seq.upper()
        except:
            raise sequenceToolsException("Invalid sequence input: Must be a string")
        
        # check we don't have any non-cannonical AAs
        bad = (["".join(x) for x in UC_seq if not x in AA])
        
        if len(bad) > 0:
            raise sequenceToolsException("Invalid sequence input: Contains non-cannonical amino acid: " + str(bad)) 
        
        # this is ugly, but provides decorator polymorphism
        if A==None:
            return fn(UC_seq)
        elif B==None:
            return fn(UC_seq,A)
        elif C==None:
            return fn(UC_seq,A,B)
        elif D==None:
            return fn(UC_seq,A,B,C)
        else:
            return fn(UC_seq,A,B,C,D)
        
    return wrapped



# ---------------------------------------
#
@_convertToUpperCase_sanitize
def calc_TCPR(seq):
    """
       Calculate the total charge per residue (TCPR)
       Return type: <float>

       seq        Protein sequence 
    """
    return(float(calc_TotalCharge(seq))/float(len(seq)))


# ---------------------------------------
#
@_convertToUpperCase_sanitize
def calc_NCPR(seq):
    """
       Calculate the net charge per residue (NCPR)
       Return type: <float>
       
       seq        Protein sequence 
    """
    return(float(calc_NetCharge(seq))/float(len(seq)))


# ---------------------------------------
#
@_convertToUpperCase_sanitize
def calc_ACPR(seq):
    """
       Calculate the aromatic amino acid content per
       residue (ACPR)
       Return type: <float>

       seq        Protein sequence 
    """
    return(float(calc_GroupContent(seq,AROMATIC))/float(len(seq)))


# ---------------------------------------
#
@_convertToUpperCase_sanitize
def calc_HCPR(seq):
    """
       Calculate the hydrophobic amino acid content per
       residue (HCPR)

       seq        Protein sequence 
    """
    return(float(calc_GroupContent(seq,HYDROPHOBIC))/float(len(seq)))


@_convertToUpperCase_sanitize
def calc_PCPR(seq):
    """
       Calculate the polar amino acid content per
       residue (HCPR)

       seq        Protein sequence 
    """
    return(float(calc_GroupContent(seq,POLAR))/float(len(seq)))


# ---------------------------------------
#
@_convertToUpperCase_sanitize
def calc_BHCPR(seq):
    """
       Calculate the bulky hydrophobe content per
       residue (HCPR)

       seq        Protein sequence 
    """
    return(float(calc_GroupContent(seq,BULKYHYDROPHOBES))/float(len(seq)))


# ---------------------------------------
#
@_convertToUpperCase_sanitize
def calc_NetCharge(seq):
    """ Calculate the net charge associated with a sequence at standard pH

        seq        Protein sequence 
    """
  
    charge = 0

    for i in seq:
        if i in NEG:
            charge = charge-1
        elif i in POS:
            charge = charge+1

    return charge

# ---------------------------------------
# ---------------------------------------
# 
@_convertToUpperCase_sanitize
def calc_TotalCharge(seq):
    """ Calculate the total charge associated with a sequence at standard pH

        seq        Protein sequence 
    """

    
    charge = 0

    for i in seq:
        if i in NEG or i in POS:
            charge = charge+1

    return charge

# ---------------------------------------
# ---------------------------------------
# 
@_convertToUpperCase_sanitize
def calc_GroupContent(seq, CLASS):
    """ Calculate the number of residues in a sequence which are in some CLASS, where class is
        a list or set of amino acids (potentially, though not necessarily from aaGroups.py

        seq        Protein sequence 
    """
    content = 0

    for i in seq:
        if i in CLASS:
            content = content+1

    return content

# ---------------------------------------
# ---------------------------------------
# 
@_convertToUpperCase_sanitize
def calc_EntropySum(seq, conf="helix"):
    """ Calculate the entopy associated with sequence based
        on standard entropy values (see aaGroups for more
        information on background).

        seq        Protein sequence 
        conf       Conformation ["helix" or "coil"]
    """

    seqType = -1

    if conf =="helix":
        seqType = HELIX_ENTROPY
    elif conf =="coil":
        seqType = COIL_ENTROPY
    else:
        raise sequenceToolsException("Invalid conformation type for calc_entropySum() funcion.")
    
    ent_sum = 0

    for i in seq:
        try:
            ent_sum = seqType[i]+ent_sum
        except KeyError, e:
            raise sequenceToolsException("Invalid amino acid in sequence")
        
    return ent_sum

# ---------------------------------------
# ---------------------------------------
# 
@_convertToUpperCase_sanitize
def calc_HydrophobicitySum(seq):
    """ Calculate the hydrophobicity associated with sequence based
        on the Kyle-Dootle hydrophobicity scale (see aaGroups for more
        information on background).

        seq        Protein sequence 

    """
    hydr_sum = 0

    for i in seq:
        try:
            hydr_sum = KD_HYDROPHOBICITY[i]+hydr_sum
        except KeyError, e:
            raise sequenceToolsException("Invalid amino acid in sequence")
        
    return hydr_sum

# ---------------------------------------
# ---------------------------------------
# 
@_convertToUpperCase_sanitize
def calc_HydrophobicityScore(seq):
    """ Calculate the hydrophobicity associated with sequence based
        on the Kyle-Dootle hydrophobicity scale (see aaGroups for more
        information on background).

        seq        Protein sequence 

    """
    maxVal = len(seq)*max(KD_HYDROPHOBICITY.values())

    hydr_sum = 0

    for i in seq:
        try:
            hydr_sum = KD_HYDROPHOBICITY[i]+hydr_sum
        except KeyError, e:
            raise sequenceToolsException("Invalid amino acid in sequence")
        
    return hydr_sum/maxVal



@_convertToUpperCase_sanitize
def calc_AntiStructureScore(seq):
    """ The anti-structure score is the summed score
    
        seq        Protein sequence 
    """
    
    score = 0
    
    for i in seq:
        score = SHEET_PROPENSITY[i]+HELIX_PROPENSITY[i]+score
        
    return score

# ---------------------------------------
# ---------------------------------------
# 
@_convertToUpperCase_sanitize
def calc_flankingScore(seq, TARGET, FLANKER, getMax=False, bonusHug=True):
    """ 
        Returns a "flanking score" (between 0 and 1) where a score
        of 1 would represent a pattern of AAs from the TARGET group
        being perfectly flanked by those in the FLANKER group.
        
        e.g. say TARGET = ('A') and FLANKER = ('C') then a perfect
        score would be generated by C[AC]* where AC repeats until 
        the end of the sequence.

        TARGET and FLANKER can be from aaGroups, or can be a list or
        tuple or string of your own choosing.

        Note that this function optimizes for the TARGET being flanked,
        but doesn't care what ELSE the FLANKER residues are doing

        seq        Protein sequence of interest
                   [String]

        TARGET     Group of residues to be treated as
                   the target of flanking [Set, List,
                   Tuple, String]

        FLANKER    Group of residues to be treated as 
                   the flankers [Set,List,Tuple,String]

        getMax     If True just return the max possible
                   score for a sequence of this length 
                   [Boolean]

        bonusHug   If True "TFT" = 6 points while if false
                   "TFT = 2 points. Basically get bonus 
                   points if a target is fully flanked.
    """

    ###################################################
    ## Programatic constants you might want to change 
    ##                                               
    BONUSHUG_SCORE = 5
    
    ##################################################
    
    ## SET the bonus val for a hugging flank
    if bonusHug:
        bval = BONUSHUG_SCORE
        maxMul = 6
    else:
        bval = 1
        maxMul=2

    ## calculate max possible value based on 
    ## sequence length (and bonusHug status)  
    if len(seq) % 2 == 0:
        maxVal = ((len(seq)/2)-1)*maxMul + 1
    else:
        maxVal = (len(seq)/2)*maxMul

    if getMax:
        return maxVal
    
    ## calculate actual flanking score
    flankercount = 0
    for i in xrange(len(seq)):

        # first residue (ugh, ugly code)
        if i == 0:
            if seq[i] in TARGET:
                if seq[i+1] in FLANKER:
                    flankercount = flankercount+1

        else:
            try:
                if seq[i] in TARGET:
                    if seq[i-1] in FLANKER:
                        flankercount = flankercount+1
                        if seq[i+1] in FLANKER:
                            flankercount = flankercount+5
                    elif seq[i+1] in FLANKER:
                        flankercount = flankercount+1

            # comes into play when at final residue, but
            # exception raised AFTER we get flanking for
            # penultimate potential flanker
            except:
                continue

    return flankercount/float(maxVal)

@_convertToUpperCase_sanitize
def calc_maximumResidueSeperation(seq, GROUP, nonDelin=False):
    """
       Returns the largest possible distance between two
       residues which are members of GROUP. 

       There are two ways to think about this, both of which
       can be calculated using this method and toggeling the
       nonDelin variable to True or False (default = False).

       Say we have a sequence AXXXAXAXXAXX and GROUP = ("A").
       The longest distance between two As without another 
       A in between is 3 (A[XXX]A). However the longest distance
       between two As irrespective of other As is 8 
       (A[XXXAXAXX]A).

       To use the former definition (i.e. without an A in
       between them) leave nonDelin=False. To get the biggest
       distance irrespective of anything else set nonDelin
       to True

       seq        Protein sequence of interest
                   [String]
                   
       GROUP      Group from which your determining distance
                  between members [Set,List,Tuple,String]

       nonDelin   See above (A[XX]A or A[XXAXX]A)
   """           
    
    numResiduesFromGroup = 0
    dif = 0
  
    for i in GROUP:
        numResiduesFromGroup = seq.count(i) + numResiduesFromGroup
        
    if numResiduesFromGroup > 1:
        
        counter=0
        pos1 = -1

        for i in seq:
            if i in GROUP:
                if pos1 == -1:
                    pos1 = counter
                else:
                    if counter - pos1 > dif:
                        dif = counter - pos1
                        
                # update the current position
                    if not nonDelin:
                        pos1 = counter
            counter = counter+1

    # dif-1 such that AA=0, AXA=1 etc...
    return dif-1

@_convertToUpperCase_sanitize
def calc_Patterning(seq, G1, G2):

      # determine our starting target
    start = 0
    switches = 0
    for i in seq:
        if i in G1:
            target1 = G2
            target2 = G1
            break
        if i in G2:
            target1 = G1
            target2 = G2
            break

        start = start+1

    seq = seq[start+1:]

    for i in seq:
   
        if i in target1:
            switches = switches+1
            target1,target2 = target2,target1
            
    return switches
