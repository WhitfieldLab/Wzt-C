import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW

output = open("<save file path>", "w+")

"""
tBLASTn Wzt O9a against translated whole genome sequences,
expect value cutoff is 10e-20 and will take up to 500000 hits
"""

blast_results = NCBIWWW.qblast("tblastn", "nr", "BAA28332.1",  expect=10e-20, hitlist_size=500000) 

output.write(blast_results.read())

output.close()

