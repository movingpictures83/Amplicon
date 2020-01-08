# Extract
# Language: Perl
# Dependency: BioPerl
# Input: TXT
# Output: FASTA
# Tested with: PluMA 1.0,  Perl 5.18

PluMA plugin that looks at sequences in a FASTA file and returns the names of organisms
whose sequences have 16S rRNA genes that generate a product with PCR.  The user must
also specify the forward and reverse primers used for PCR.

The plugin takes as input a TXT file with the following information, each on separate lines:

Input FASTA file name
Output FASTA file name (for PCR fragment)
TXT file (for PCR fragment length)
Forward Primer
Reverse Primer

Names and phylogenetic information of all matching organisms will be output to a TXT file.


