import unittest
import os
from Bio import SeqIO


class TestStringMethods(unittest.TestCase):

  def test_existence(self):
      self.assertTrue(os.path.isfile('./test_out/Hyalella_azteca_test/Hyalella_azteca_test_pia3_pep.fasta'))
      self.assertTrue(os.path.isfile('./test_out/Parhyale_hawaiensis_test/Parhyale_hawaiensis_test_pia3_pep.fasta'))

  def test_file(self):
      seq_num = 0
      for seq_record in SeqIO.parse('./test_out/Hyalella_azteca_test/Hyalella_azteca_test_pia3_pep.fasta', "fasta"):
          seq_num += 1
      self.assertTrue(seq_num == 3)
      
      seq_num_2 = 0
      for seq_record in SeqIO.parse('./test_out/Parhyale_hawaiensis_test/Parhyale_hawaiensis_test_pia3_pep.fasta', "fasta"):
          seq_num_2 += 1
      self.assertTrue(seq_num_2 == 2)
      
      
if __name__ == '__main__':
    unittest.main()
