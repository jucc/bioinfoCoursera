#!/bin/python
### Functions for manipulating frequency of k-mers (patterns) in genomes (text)
### Following week 1 of Coursera's Bioinformatics Algorithms courser
### https://www.coursera.org/learn/dna-analysis/gradedLti
### Julia

def allKmersFromText(text, k):
    """
    Returns a list with all k-mers (strings with length k) in the text (including repetitions)
    """
    return [ text[i:i+k] for i in range ( len(text) - k + 1 ) ]


def kmersFrequency(text, k):
    """
    Creates a dictionary with the frequency of each k-mer in the text
    """
    freqs = {}
    for kmer in allKmersFromText(text, k):
        if kmer in freqs:
            freqs[kmer] += 1
        else:
            freqs[kmer] = 1
    return freqs


def patternCount(text, pattern):
    """
    Counts how many times a k-mer pattern appears in the text
    """
    return kmersFrequency(text, len(pattern))[pattern]


def mostFrequentKmers(text, k):
    """
    Returns a list with the k-mers with the highest frequency n the text
    """
    freqs = kmersFrequency(text, k)
    maxFreq = max(freqs.values())
    return [k for k,v in freqs.items() if v == maxFreq]


def complement(pattern):
    """
    Finds the complement of a pattern, formed by taking the complement of each nucleotide
    """
    comp = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c' }
    return ''.join(comp[c] for c in pattern)


def reverseComplement(pattern):
    """
    Finds the reverse complement of a pattern, formed by taking the complement of each nucleotide,
    then reversing the resulting string
    """
    return complement(pattern)[::-1]


def patternMatch(pattern, genome):
    """
    Returns the starting indices of all occurrences of pattern in genome
    """
    index = 0
    matches = []
    for word in allKmersFromText(genome, len(pattern)):
        if word == pattern:
            matches.append(index)
        index += 1
    return matches


def slide(window, freqs, k):
    """
    Slides the window by removing the first k-mer from the left and adding a new
    k-mer with the new character added in the right.
    Returns the included k-mer and the freqs dictionary is updated accordingly
    Attention: The input window must start in the previous position (before sliding)
    and finish in the new (after sliding) position, thus, its len is L+1
    """
    # one being excluded from the left...
    kmerExclude = window[0:k]
    freqs[kmerExclude] -= 1

    # ...and one included in the right
    kmerInclude = window[len(window)-k:len(window)]
    if kmerInclude in freqs:
        freqs[kmerInclude] += 1
    else:
        freqs[kmerInclude] = 1
    return kmerInclude


def findClumps(genome, k, L, t):
    """
    Finds which k-mers appear at least t times in any window of lenght L inside
    the genome text
    """
    # figures out the clumps and freqs in the initial window
    freqs = kmersFrequency(genome[0:L], k)
    clumps = set(filter(lambda kmer: freqs[kmer] >= t, freqs))

    # slides the window and checks if the new kmer in the right forms a clump
    for windowStart in range(len(genome) - L + 1):
        kmer = slide(genome[windowStart:windowStart+L], freqs, k)
        if freqs[kmer] >= t:
            clumps.add(kmer)

    return list(clumps)


def patternToNumber(pattern):
    """
    Converts a pattern (base 4 number where 0 = A, 1 = C, 2 = G, 3 = T) into a
    base 10 number.
    """
    base = { 'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3 }
    exp = 0
    number = 0

    for ch in pattern[::-1]:
        number += base[ch] * (4 ** exp)
        exp += 1
    return number


def computeFrequencyArray(text, k):
    """
    This is not a memory efficient structure, as I had already defined a hash for
    keeping the frequencies, but is required in the course
    """
    freqArray = [0 for i in range(4**k)]

    freqs = kmersFrequency(text, k)
    for f in freqs:
        freqArray[patternToNumber(f)] = freqs[f]

    return freqArray
