import unittest
import os
from Bio import SeqIO


class TestStringMethods(unittest.TestCase):

  def test_existence(self):
      self.assertTrue(os.path.isfile('./test_out/test_Hyalella_azteca/test_Hyalella_azteca_PIA3_aa.fasta'))

  def test_file(self):
      seq_num = 0
      opsins_type = []
      for seq_record in SeqIO.parse('./test_out/test_Hyalella_azteca/test_Hyalella_azteca_PIA3_aa.fasta', "fasta"):
          opsins_type.append(seq_record.id[0:3])
          seq_num += 1
      print(opsins_type)
      self.assertTrue(seq_num == 3)
      self.assertTrue(sorted(opsins_type) == ['LWS', 'LWS', 'MWS'])
      
      
if __name__ == '__main__':
    unittest.main()
