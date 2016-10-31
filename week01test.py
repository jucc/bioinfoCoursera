from week01functions import *
import unittest

class TestPdfToTable(unittest.TestCase):

    def setUp(self):
        a = 1

    def test_get_all_kmers(self):
        genome = "GATTACA"
        expected = ['GAT', 'ATT', 'TTA', 'TAC', 'ACA']
        result = allKmersFromText(genome, 3)
        self.assertEquals(expected, result)


    def testGetKmersFrequency(self):
        genome = "GCGCG"
        expected = { 'GCG': 2, 'CGC': 1 }
        result = kmersFrequency(genome, 3)
        self.assertEquals(expected, result)


    def testPatternCount(self):
        genome = "GCGCG"
        result = patternCount(genome, 'GCG')
        self.assertEquals(2, result)


    def testFrequentWords(self):
        genome = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        expected = ['CATG', 'GCAT']
        result = mostFrequentKmers(genome, 4)
        self.assertEquals(expected, result)


    def testReverseComplement(self):
        pattern = "AAAACCCGGT"
        expected = "ACCGGGTTTT"
        self.assertEquals(expected, reverseComplement(pattern))


    def testReverseComplementLowercase(self):
        pattern = "aaaacccggt"
        expected = "accgggtttt"
        self.assertEquals(expected, reverseComplement(pattern))


    def testPatternMatch(self):
        genome = "GATATATGCATATACTT"
        pattern = 'ATAT'
        expected = [1, 3, 9]
        result = patternMatch(pattern, genome)
        self.assertEquals(expected, result)


    def testPatternMatchLarge(self):
        try:
            file = open('datasets/Vibrio_cholerae.txt', 'rb')
        except:
            print ("Could not open input file vibrio cholerae")
        genome = file.read()
        file.close()
        expected = [116556, 149355, 151913, 152013, 152394, 186189, 194276,
                    200076, 224527, 307692, 479770, 610980, 653338, 679985,
                    768828, 878903, 985368]
        result = patternMatch('ATGATCAAG', genome)
        self.assertEquals(expected, result)

    #
    # def testSlide(self):
    #     genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
    #     L = 50
    #     freqs =  kmersFrequency(genome, 5)
    #     window = genome[0:L+1]
    #     excluded = genome[0:5]
    #     included = genome[L+1-5:L+1]
    #     fe1 = freqs[excluded]
    #     fi1 = freqs[included]
    #     slide(window, freqs, 4)
    #     fe2 = freqs[excluded]
    #     fi2 = freqs[included]
    #     assertTrue( fe1 - fe2 == 1)
    #     assertTrue( fi2 - fi1 == 1)



if __name__ == '__main__':
    unittest.main()
