import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

Entrez.email = "<email address>"
Entrez.api_key ="<Entez key>"

"""
This reads in the genbank accession and WztO9a HSP position from text file.
It downloads the genbank file. It then uses the position to identify
the Wzt CDS feature and checks if it is between 325 and 500 amino
acids (975 to 1650 nucleotides). If the length is within this range,
the name of the organism name, the taxid id and the Wzt protein id,
protein description and translation are recorded.
The description of each CDS feature 7500 bp upstream and 7500 bp
downstream of the wzt start position is then read. If any of these
contains a "glyco" key word, then recorded Wzt information is output
to a csv file.
"""

reference_file = open("<read file path>","r") #Contains genome id blast hits and hsp positions
output = open("<save file path>","w+") #Output file for descriptions of extended-Wzt
exceptions = open("<save exceptions path>","w+") #Output file for errors

output.write("count, gb id, description, protein id, length\n") #Output csv headers

genbank_ids=[] #Container for list of genbank ids from input text file "ref_file"
position_wzt_hsp=[] #Container for list of wzt hsp positions from input text file "ref_file"

"""Reads the input file and assigns genbank ids and wzt positions to lists"""
for line in reference_file: 
    genbank_ids.append(str(line.split(",")[0]))
    position_wzt_hsp.append(str(line.split(",")[1])[1:])
   
list_count = -1 #Needed so that info can be obtained from list based on index
neighbourhood = "" #Container for annotation of genes in wzt neighbourhood
organism ="" #Container for organism description
taxid ="" #Container for organism taxid id

"""common words used in annotation for genes in polysaccharide synthesis gene clusters"""
keyword = ['glycosyl', 'glycosyltransferase', 'manno', 'polysaccharide', 'wb', 'lps', 'lipopolysaccharide','wzt', 'tagh', 'glucose', 'rhamnose', 'gt']

for item in genbank_ids: 
    try:
        list_count+=1
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=item) #Fetch genbank file
        seq_record = SeqIO.read(handle, "gb") 
        SeqIO.write(seq_record, "C:/Users/Evan/Desktop/temp_hit.gb", "gb") #Save full genbank file locally (overwritten in next loop)
        seq_record1 = SeqIO.read("C:/Users/Evan/Desktop/temp_hit.gb", "gb") #Read local genbank file
        PSS_cluster = False #Set to false until it is confirmed that Wzt is in a polysaccharide synthesis gene cluster
        wzt_end_position=0 #Container for wzt CDS feature end position
        strand = 0 #Container to record if wzt on forward or reverse strand

        for feature in seq_record1.features: #Scans through each feature in the genbank file

            if feature.type =="source": #Source features contain data related to the organism and its genome
                organism = str(feature.qualifiers["organism"]).split("'")[1] #Records name of the organism
                                
                organism_edit ="" #Container for modified organism name
                for letter in organism: #This replaces ',' with ';' for csv file output
                    if letter != ",":
                        organism_edit+= letter
                    else:
                        organism_edit+= ";"            

            #Find wzt based on the position of HSP start from position_wzt_hsp list lying between the CDS feature start and end position
            elif (int(feature.location.start.position <= int(position_wzt_hsp[list_count])) and int(feature.location.end.position >= int(position_wzt_hsp[list_count]))) and feature.type == "CDS" and abs(feature.location.end.position-feature.location.start.position)<5000:

                """If wzt is between 975 and 1650 nt then the genome id, organism name, taxid id, protein id and translation are recorded in a string"""
                if (975 < abs(feature.location.end.position-feature.location.start.position) < 1650): #Checks if wzt length between 325 and 550 aa
                    wzt_end_position = int(feature.location.end.position) #Records wzt CDS feature end position
                    strand = feature.strand #Records if wzt on forward or reverse strand

                    try:
                        translation = str(feature.qualifiers["translation"]).split("'")[1] #Records wzt amino acid sequence
                    except:
                        translation ="" #If there is a problem obtaining translation then an empty string is recorded
                    genome_id = str(seq_record1.id) #Records the name of the genbank file

                    genome_id_edit = "" #Container for modified genbank file name
                    for letter in genome_id: #This replaces ',' with ';' for csv file output
                        if letter != ",":
                            genome_id_edit+= letter
                        else:
                            genome_id_edit+= ";"
                    try:
                        protien_id = str(feature.qualifiers['protein_id']).split("'")[1] #Records wzt protein id
                    except:
                        protien_id="" #If there is a problem obtaining protein id then an empty string is recorded

                    protein_id_edit = "" #Container for modified protein id
                    for letter in protien_id: #This replaces ',' with ';' for csv file output
                        if letter != ",":
                            protein_id_edit+= letter
                        else:
                            protein_id_edit+= ";"

                   #Wzt_extended_description is a string containing the Wzt information (oranism, taxid, protein id, translation)
                    Wzt_extended_description = (str(list_count)+", "+(genome_id_edit)+", "+(str(organism_edit))+", "+protein_id_edit+ ", "+ str(len(translation)))

        #Finds CDS features on the same strand as wzt and within 7500 nt of wzt and records the annotation of each in the neighbourhood string"""
        for feature in seq_record.features:
            if feature.strand == strand and ((feature.type =="CDS" and feature.location.start.position >= wzt_end_position and feature.location.start.position <=(wzt_end_position+7500)) or (feature.type == "CDS" and feature.location.start.position <= wzt_end_position and feature.location.start.position >=(wzt_end_position-7500))):
                neighbourhood+=(str(feature.qualifiers['product']).split("'")[1]+" | ")
                
        #Scans the neighbourhood string for the presence of a glyco-related keyword"""        
        for word in keyword:
            if word in str(neighbourhood).lower(): #If a keyword is in the neighbourhood...                            
                PSS_cluster = True #then the gene cluster is classified as polysaccharide synthesis cluster.
            neighbourhood=''
        if PSS_cluster: #If wzt is in a polysaccharide synthesis cluster then its description is output in a csv file. 
            output.write(str(Wzt_extended_description) + "\n")
            PSS_cluster = False
              
    except Exception as e:
        print(str(e))
        exceptions.write(str(e)+ " "+str(item)+"\n") #Writes errors to text file

output.close()
reference_file.close()
exceptions.close()

