# Histone-3.3-protein-HMM

Project for CS123B at SJSU to build a protein HMM for Histone 3.3

The script is broken into sections and works as follows:
  1. BLASTs Homo Sapiens Histone 3.3 sequence (ref|NP_001005101.1|) and saves results in XML file
  2. Parses through XML file and saves sequences in a multi-sequence FASTA file
  3. Does a MSA using clustalo and saves aligned FASTA file 
  4. Splits aligned file into a training and test set 
  5. Converts fasta files into PHYLIP format so that it could be input into a phylogenetic web tool 
  
  A phylogenetic tree was generated using PhyML and HMM was built using external jar file.
  
  Dependences: 
  Python 3+ 
  ClustalOmegaCommandLine
