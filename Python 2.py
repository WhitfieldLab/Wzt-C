import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

"""
This reads Blast results from xml file and writes the genbank id and HSP
position to an output text file.
"""

blast_results = open("<read file path>", "r")
output = open("<save file path>", "w+")

hsps=[]
blast_records = NCBIXML.parse(blast_results)
for hit in blast_records:
    for alignment in hit.alignments:
        for hsp in alignment.hsps:
            if hsp not in hsps:
                hsps.append(hsp)
                output.write(str(alignment.title.split("|")[3] + ", " +str(hsp.sbjct_start)+"\n"))

    
blast_results.close()
output.close()
