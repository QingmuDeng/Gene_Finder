# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


completment_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    return completment_dict[nucleotide.upper()]


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
    # Create a empty string to hold the reverse complement
    reverse_complement = ''

    # Add the nucleotide complements to the reverse complement string
    for nucleotide in dna:
        reverse_complement += get_complement(nucleotide)
    return reverse_complement[::-1]


start_codon = codons[3][0]
stop_codon = codons[10]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i = 0
    ORF = ''

    # As long as the current position of the curser is not beyond the dna length,
    # read three as a group until hit a stop codon
    while i < len(dna):
        # Save the current i as the index to start recording 3 necleuctides
        old_i = i

        # Add 3 to i to get the end index
        i += 3

        # Store the 3 necleuctides as one in a temporary variable
        temp = dna[old_i:i:1]

        # check if the codon is a stop codon, return the ORF as is; otherwise
        # keep adding codons to the ORF
        if temp == stop_codon[0] or temp == stop_codon[1] or temp == stop_codon[2]:
            return ORF
        else:
            ORF += temp
    return ORF


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    ORF_list = []

    # base_frame is the base from which to look for a starting codon
    base_frame = dna

    # base_frame_num is the starting index for a start codon
    base_frame_num = base_frame.find(start_codon)

    # Making sure that the frames stay in sync by making sure the ATG location
    # can be divided by three
    while base_frame_num % 3 != 0:
        base_frame = base_frame[3:]
        base_frame_num = base_frame.find(start_codon)

    # Keep looking for an ORF until when there's no more to read or when it
    # would be a frame out of the original sync
    while 1:
        # base_frame_num would be negative when there is no starting codon to
        # be found, in this case, the loop should break and allow the function
        # to return the ORF collected up to this point
        if base_frame_num < 0:
            break

        # Get one full frame from rest_of_ORF() and append it to the ORF_list
        full_frame = rest_of_ORF(base_frame[base_frame_num:])
        ORF_list.append(full_frame)

        # Look for the next ATG index from a base_frame with the previously
        # found full_frame removed from the start
        next_ATG_loc = base_frame[len(full_frame):].find(start_codon)

        # if the index of the next ATG is not divisible by 3 (meaning that
        # the new start codon isn't in sync with the previuos one), then
        # break from the looking for another ORF
        if (next_ATG_loc % 3) != 0:
            break
        else:
            base_frame = base_frame[len(full_frame) + 3:]
            base_frame_num = next_ATG_loc - 3

    return ORF_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORF = []

    # Loop through the possibilities of a starting codon being on the multiples
    # of 0, 1, 2 index
    for i in range(3):
        temp_all_frame = find_all_ORFs_oneframe(dna[i:])

        # In case that there are more than just one frame in one find_all,
        # append individual ones to the list with the for loop
        for y in range(len(temp_all_frame)):
            all_ORF.append(temp_all_frame[y])

    return all_ORF


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
