import unittest
from unittest import TestCase


class TestNucleoBaseArray(TestCase):
	REF_FILE = 'ref_test.txt'
	PRIV_FILE = 'priv_test.txt'
	ANS_FILE = 'ans_test.txt'
	READS_FILE = 'reads_test.txt'

	def setUp(self):
		pass

	# @unittest.skip('')
	def test_inserts(self):
		with open(self.REF_FILE) as ref_genome, \
		open(self.PRIV_FILE) as priv_genome, open(self.READS_FILE) as reads:
			pass

	@unittest.skip('')
	def test_deletes(self):
		with open(self.REF_FILE) as ref_genome, \
		open(self.PRIV_FILE) as priv_genome, open(self.READS_FILE) as reads:
			pass

	@unittest.skip('')
	def test_snps(self):
		with open(self.REF_FILE) as ref_genome, \
		open(self.PRIV_FILE) as priv_genome, open(self.READS_FILE) as reads:
			pass


if __name__ == "__main__":
	unittest.main()