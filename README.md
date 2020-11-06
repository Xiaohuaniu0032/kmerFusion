# kmerFusion
detect DNA fusion from paired-end sequencing using k-mer method

## Method
1. for target gene, get its fasta and make a k-mer list, this k-mer list will store this k-mer's genome location
2. for each read, get its all k-mer sequences (forward and reverse complement) and map it to the `#1` k-mer list
3. 