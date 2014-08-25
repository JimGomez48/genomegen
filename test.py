import unittest
from unittest import TestCase


class TestGenomeGen(TestCase):
    REF_FILE = 'ref_test.txt'
    PRIV_FILE = 'priv_test.txt'
    ANS_FILE = 'ans_test.txt'
    READS_FILE = 'reads_test.txt'

    def setUp(self):
        self.ref_genome = []
        self.priv_genome = []
        self.answers = []
        with open(self.REF_FILE, mode='r', buffering=1) as ref_file:
            chrom_num = -1
            for line in ref_file:
                if 'test' in line:
                    continue
                if '>' in line:
                    chrom_num += 1
                    continue
                try:
                    self.ref_genome[chrom_num].extend(list(line.rstrip()))
                except StandardError:
                    self.ref_genome.append(list(line.rstrip()))
        with open(self.PRIV_FILE, mode='r', buffering=1) as priv_file:
            chrom_num = -1
            for line in priv_file:
                if 'test' in line:
                    continue
                if '>' in line:
                    chrom_num += 1
                    continue
                try:
                    self.priv_genome[chrom_num].extend(list(line.rstrip()))
                except StandardError:
                    self.priv_genome.append(list(line.rstrip()))

    # @unittest.skip('')
    def test_inserts(self):
        CHROM = 0
        SEQUENCE = 1
        INDEX = 2
        with open(self.ANS_FILE, mode='r', buffering=1) as ans_file:
            # seek up until the INSERT tag
            while 1:
                line = ans_file.readline()
                if not line:
                    self.fail('INSERT flag not found')
                if 'INSERT' in line:
                    line = ans_file.readline()
                    break

            iter = 0
            while line:
                if '>' in line:
                    break
                tokens = line.rstrip().split(',')
                chrom = int(tokens[CHROM])
                ans_insert = str(tokens[SEQUENCE])
                index = int(tokens[INDEX])

                ref_prefix = ''.join(self.ref_genome[chrom][index-5: index])
                ref_suffix = ''.join(self.ref_genome[chrom][index: index+5])
                priv_prefix = ''.join(self.priv_genome[chrom][index - 5: index])
                priv_insert = ''.join(self.priv_genome[chrom][index: index + len(ans_insert)])
                priv_suffix = ''.join(self.priv_genome[chrom][index + len(ans_insert): index + len(ans_insert) + 5])

                # print 'ans: ' + ans_insert + '  priv: ' + priv_insert
                # self.assertEquals(ref_prefix, priv_prefix, 'iter:' + str(iter) + ' ' + ref_prefix + '!=' + priv_prefix)
                self.assertEquals(ans_insert, priv_insert, 'iter:' + str(iter) + ' ' + ans_insert + '!=' + priv_insert)
                # self.assertEquals(ref_suffix, priv_suffix, 'iter:' + str(iter) + ' ' + ref_suffix + '!=' + priv_suffix)

                line = ans_file.readline()
                iter += 1


    @unittest.skip('')
    def test_deletes(self):
        with open(self.REF_FILE, mode='r', buffering=1) as ref_genome, \
            open(self.PRIV_FILE, mode='r', buffering=1) as priv_genome, \
                open(self.ANS_FILE, mode='r', buffering=1) as ans_file:
            # seek into the start of the deletes portion
            while 1:
                ans_line = ans_file.readline().rstrip()
                if '>DELETE' in ans_line:
                    break
                if not ans_line:
                    self.fail("Did not find an DELETE flag in the answer key")

    @unittest.skip('')
    def test_snps(self):
        with open(self.REF_FILE, mode='r', buffering=1) as ref_genome, \
            open(self.PRIV_FILE, mode='r', buffering=1) as priv_genome, \
                open(self.ANS_FILE, mode='r', buffering=1) as ans_file:
            # seek into the start of the snp portion
            while 1:
                ans_line = ans_file.readline().rstrip()
                if '>SNP' in ans_line:
                    break
                if not ans_line:
                    self.fail("Did not find an SNP flag in the answer key")


if __name__ == "__main__":
    unittest.main()