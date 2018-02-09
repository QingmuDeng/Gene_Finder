# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
import multiprocessing as mp



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
    >>> get_reverse_complement('TGAATGTAG')
    'CTACATTCA'
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
    >>> rest_of_ORF('ATGATAGATTAGGTGCCC')
    'ATGATAGAT'
    """
    i = 0
    ORF = ''

    # As long as the current position of the curser is not beyond the dna length,
    # read three as a group until hit a stop codon
    while i < len(dna):
        #print(dna)
        # Save the current i as the index to start recording 3 necleuctides
        old_i = i
        # Add 3 to i to get the end index
        i += 3
        # store the 3 necleuctides as one in a temporary variable
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
    >>> find_all_ORFs_oneframe('GCATGAATGTAG')
    ['ATG']
    """
    ORF_list = []
    #base_num = dna.find('ATG')
    base_frame = dna
    while 1:
        base_num = base_frame.find('ATG')
        if base_num < 0:
            return ORF_list
        elif base_num % 3 != 0:
            base_frame = base_frame[3:]
        else:
            #base_num
            base_frame = base_frame[base_num:]
            one_frame = rest_of_ORF(base_frame)
            base_frame = base_frame[len(one_frame)+3:]
            ORF_list.append(one_frame)

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
    >>> find_all_ORFs("GCATGAATGTAG")
    ['ATG', 'ATGAATGTAG']
    """
    all_ORF = []
    # Loop through the possibilities of a starting codon being on the multiples
    # of 0, 1, 2 index
    for i in range(3):
#         print(dna[i:])
        temp_all_frame = find_all_ORFs_oneframe(dna[i:])
        # In case that there are more than just one frame in one find_all,
        # append individual ones to the list with the for loop
        if temp_all_frame != None:
            all_ORF += [y for y in temp_all_frame]

    return all_ORF



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG', 'ATGCAT']
    """
    dna_complement = get_reverse_complement(dna)
    original_ORF = find_all_ORFs(dna)
    complement_ORF = find_all_ORFs(dna_complement)
    ORF_both = [i for i in original_ORF]
    ORF_both += [i for i in complement_ORF]

    return ORF_both


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCATGAATGGCATGAATGTAG")
    'ATGCATGAATGGCATGAATGTAG'
    """
    ORF_both = find_all_ORFs_both_strands(dna)
    longest_length = 0
    longest_length_index = None
    for i in range(len(ORF_both)):
        if len(ORF_both[i]) > longest_length:
            longest_length_index = i
            longest_length = len(ORF_both[i])
    if longest_length_index == None:
        return '0'
    else:
        return ORF_both[longest_length_index]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    random.seed()

    # Define an output queue
    output = mp.Queue()

    def rand_dna(dna, output):
        """ Generates a random string of numbers, lower- and uppercase chars. """
        length = len(dna)
        rand_str = ''.join(random.sample(dna, length))
        rand_str = longest_ORF(rand_str)
        output.put(rand_str)

    times_to_run = 1
    if num_trials > 100:
        times_to_run = num_trials / 100
        num_trials = 100

    results = []

    for i in range(int(times_to_run)):
        # Setup a list of processes that we want to run
        shuffle_process = [mp.Process(target=rand_dna, args=(dna, output)) for num in range(num_trials)]

        # Run processes
        for p in shuffle_process:
            p.start()

        # Exit the completed processes
        for p in shuffle_process:
            p.join()

        # Get process results from the output queue
        results.append(max([output.get() for p in shuffle_process], key = len))

    return len(max(results, key = len))


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
        >>> coding_strand_to_AA("ATGCATGAATGGCATGAATGTAG")
        'MHEWHEC'
    """
    AA = ''
    i = 0
    while i < len(dna):
        old_i = i
        i += 3
        if i > len(dna):
            break
        AA += aa_table[dna[old_i:i:1]]
    return AA


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    AA_sequence = []
    threshold = longest_ORF_noncoding(dna, 1500)
    print(threshold)
    both_ORF = find_all_ORFs_both_strands(dna)
    for i in range(len(both_ORF)):
        if threshold < len(both_ORF[i]):
            AA_sequence.append(coding_strand_to_AA(both_ORF[i]))
    print(AA_sequence)
    return AA_sequence

# def main():
    # print("Hello Gene Finder")
from load import load_seq
dna = load_seq("./data/X73525.fa")
gene_finder(dna)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
