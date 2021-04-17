"""
Author: Roshani Codipilly
Date last modified: April 16, 2021

Script to run BLAST on Histone 3, a protein important for DNA packaging. 50 sequences of closely related seqs (all
mammals) will be collected as "basic set", 30 seqs from a different clade (birds) will be the "related set".
Then a MSA will be used to align the sequences and a phylogenetic tree will be built to visualize diversity.
Additionally, the related and basic sets will be divided into training and test sets to use used to train and
test an HMM

"""


from Bio.Blast import NCBIWWW                                    #library used to perform BLAST search
from Bio.Blast import NCBIXML                                    #library used to read in BLAST results xml
from Bio import SeqIO                                            #library to read in fasta files
from Bio.Align.Applications import ClustalOmegaCommandline       #library used to perform MSA
from random import sample                                        #library used to randomly split test and training
from Bio import Phylo                                            #library to build the phylogenetic tree
import sys                                                       #library used for exporting the phylo tree as a txt file
from Bio import AlignIO



#Section 1a: running BLAST on histone 3.3 to find sequences for basic set (taxid limited to mammals, taxid:40674)

"""
Running BLAST through python in order to identify closely related mammal species that will make up the 
basic set and the bird species that will make up the related sets
"""

#open and store the seqeunce data for the histone 3.3 protein
hist3_seq_record = next(SeqIO.parse(open('Hist3.3_proteinSeq.txt'), 'fasta'))

#call BLAST on the sequence record, which will be stored in a Seq object
qblast_BasicResults = NCBIWWW.qblast("blastp", "nr", hist3_seq_record.seq, entrez_query='txid40674[ORGN]')

#open Seq object and read into a file
with open('BLAST_BasicResults.xml', 'w') as f:
    blast_BasicResults = qblast_BasicResults.read()
    f.write(blast_BasicResults)


#Section 1b: running BLAST on histone 3.3 to find sequences for related set (taxid limited to birds, taxid:8782 )

#call BLAST on the sequence record, which will be stored in a Seq object
qblast_RelResults = NCBIWWW.qblast("blastp", "nr", hist3_seq_record.seq, entrez_query='txid8782[ORGN]')

#open Seq object and read into a file
with open('BLAST_RelResults.xml', 'w') as f:
    blast_RelResults = qblast_RelResults.read()
    f.write(blast_RelResults)

#Section 2a: Parsing through BLAST results and saving the results to FASTA file - Basic set

#open BLAST results, then parse through them
blast_BasicRecords = NCBIXML.parse(open('BLAST_BasicResults.xml', 'r'))

#set a threshold for the E value
E_value_thres = 0.04
 
#opening a file to write in BLAST results in fasta format to be used for MSA
fasta_BasicSeq = open("BasicSeq.fasta", 'w')

#loops through  matches, displays summary info, then stores seq header and sequence in dictionary
for blast_record in blast_BasicRecords:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_value_thres:                      
                fasta_BasicSeq.write('>%s\n%s\n' % (alignment.title.split()[0], hsp.sbjct))


#important to close file to protect integrity
fasta_BasicSeq.close()

#Section 2b: Parsing through BLAST results and saving results to FASTA file - Related set

#open BLAST results, then parse through them
blast_RelRecords = NCBIXML.parse(open('BLAST_RelResults.xml', 'r'))

#opening a file to write in BLAST results in fasta format to be used for MSA
fasta_relSeq = open("RelSeq.fasta", 'w')

#loops through  matches, displays summary info, then stores seq header and sequence in dictionary
for blast_record in blast_RelRecords:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_value_thres:                     #same threshold as basic set
                fasta_relSeq.write('>%s\n%s\n' % (alignment.title.split()[0], hsp.sbjct))


#important to close to protect integrity
fasta_relSeq.close()

#Section 3: Performing MSA

"""
Multiple sequence alignment performed on the basic and related sets, which will then be written out into an aligned
fasta file. It is important for fasta files to be aligned because the profile HMM will be unable to run if the 
sequences are different lengths
"""

#before doing MSA, basic and related seq sets were first combined into total seq data set on the command line


#defining the input and output files for the total set
in_file = 'totalSeq.fasta'
out_file = 'totalSeq_aligned.fasta'

#the command line command for converting the fasta file to a aligned fasta file
clineBasic = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
clineBasic()

#now that alignment is complete, total data set will now be manually separated back into basic and related aligned datasets


#Section 4: Dividing the aligned files into training and test sets

"""
Splitting aligned fasta files into training and testing sets, which will be used to train and test the
profile HMMs. Training is 80% of set, testing is 20%
"""


#lists to hold align sequences for basic and related sets
BasicAligned_list = list(SeqIO.parse("BasicSeq_aligned.fasta", "fasta"))
RelAligned_list = list(SeqIO.parse("RelSeq_aligned.fasta", "fasta"))

#randomly selecting 80% of sequences to be training set
BasicAligned_train = sample(BasicAligned_list, 40)
RelAligned_train = sample(RelAligned_list, 40)

#create lists of record ids that are already in training set

#for basic set
BasicID_list = []
for i in range(0, len(BasicAligned_train)):
    BasicID_list.append(BasicAligned_train[i].id)

#for related set
RelID_list = []
for i in range(0, len(RelAligned_train)):
    RelID_list.append(RelAligned_train[i].id)

#remaining 20% makes up the test set, which is found by searching for the ids not found in the training set ids

#empty list to store the training set
BasicAligned_test = []

for i in range(0, len(BasicAligned_train)+1):
    if(BasicAligned_list[i].id not in BasicID_list):
        BasicAligned_test.append(BasicAligned_list[i])

#empty list to store the training set
RelAligned_test = []

for i in range(0, len(RelAligned_train)+1):
    if(RelAligned_list[i].id not in RelID_list):
        RelAligned_test.append(RelAligned_list[i])

#write out the training sets into a fasta file
SeqIO.write(BasicAligned_train, "BasicAlignedTrain.fasta", "fasta")
SeqIO.write(RelAligned_train, "RelAlignedTrain.fasta", "fasta")

#write out the test sets into a fasta file
SeqIO.write(BasicAligned_test, "BasicAlignedTest.fasta", "fasta")
SeqIO.write(RelAligned_test, "RelAlignedTest.fasta", "fasta")


#Section 5: Building a phylogenetic tree

"""
Building a phylogenetic tree of both the basic and related data sets in order to gauge the level of diversity 
between all the samples. Since I was interested in the taxonomic diversity, the genus and species associated with 
each seqeuence was first mapped to each accession number and then a phylogenetic tree was generated using Clustalw
"""

#Section 5a: opening BLAST xml results again in order to collect taxomic information for each seq in basic and related

#opening BLAST XML
blast_BasicRecords = NCBIXML.parse(open('BLAST_BasicResults.xml', 'r'))

# set a threshold for the E value
E_value_thres = 0.04

#dictionary to store species taxonomic info to corresponding sequence id
basic_species_dictionary = {}

# loops through  matches, then stores sequence id and taxonomic info in dictionary
for blast_record in blast_BasicRecords:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_value_thres:
                startIndex = alignment.title.find("[")
                endIndex = alignment.title.find("]")
                species = alignment.title[startIndex:endIndex + 1]
                sequence_id = alignment.title.split()[0]
                basic_species_dictionary[sequence_id] = species

#now repeating those same steps for the related data set

# opening BLAST xml
blast_RelRecords = NCBIXML.parse(open('BLAST_RelResults.xml', 'r'))

#dictionary to store sequence id: taxonomic information
rel_species_dictionary = {}

# loops through  matches, then stores sequence id and taxonomic info in dictionary
for blast_record in blast_RelRecords:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_value_thres:
                startIndex = alignment.title.find("[")
                endIndex = alignment.title.find("]")
                species = alignment.title[startIndex:endIndex + 1]
                sequence_id = alignment.title.split()[0]
                rel_species_dictionary[sequence_id] = species

#channging the header of the fasta files to be the taxonomic information instead of the accession number

#combining the two dictionaries together so they have all the information
def merge(dict1, dict2):
    return(dict1.update(dict2))

merge(basic_species_dictionary, rel_species_dictionary)


#now basic species dictionary is updated to contain all names, store in another variable instead
total_specieNames = basic_species_dictionary
total_seq = {}

#opening aligned FASTA file and writing it into total seq dictionary
for record in SeqIO.parse("totalSeq_Aligned.fasta", "fasta"):
    total_seq[record.id] = record.seq


#loop through all the keys in the total seq dictionary and replace it with the taxonomic information instead
#note - you have to loop through a copy of the dictionary in order to iterate through and change the values
for key in total_seq.copy().keys():
    total_seq[total_specieNames[key]] = total_seq.pop(key)


#now write the new dictionary into a FASTA file
with open("totalSeqTax_aligned.fasta", 'w') as f:
    for key, value in total_seq.items():
        f.write(">" + key + '\n' + str(value) + '\n')

#5b - Converting FASTA file to phylip file in order to build tree on PhyML

AlignIO.convert("totalSeqTax_alignedEdited2.fasta", "fasta", "editedTax.phy", "phylip")


