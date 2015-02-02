import time

def get_complement(nucleotide):
    """ Returns the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    complements = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return complements[nucleotide]
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    complement = ""
    for char in dna:
        complement += get_complement(char)
    return (complement)[::-1]
    pass

# print get_reverse_complement('AAAGCGGGCAT')


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop_codons = ["TAG","TAA","TGA"]
    output = ""
    for i in range(0,len(dna),3):
        if(dna[i:i+3] not in stop_codons):
            output += dna[i:i+3]
        else:
            break
    return output
    pass

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    # orfs = [rest_of_ORF(dna)]
    # dna = dna[len(orfs[0]):]
    orfs = []
    while(len(dna) >= 3):
        while(dna[:3] != 'ATG' and len(dna) >= 3):
            dna = dna[3:]
        if(len(dna)<3):
            pass
        to_add = rest_of_ORF(dna)
        if(len(to_add)>=3):
            orfs.append(to_add)
        dna = dna[len(to_add):]

    return orfs
    pass

# print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    return find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])     
    pass


# print find_all_ORFs("ATGCATGAATGTAG")

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    pass

print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")