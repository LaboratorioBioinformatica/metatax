### MetaTax: Meta-classification of taxonomy classification ###

This is the source code of **MetaTax**, a meta-classification algorithm of taxonomy classification. 

The input to MetaTax is a file with reads names (metagenomics) and several taxonomic classification
for those reads that come from different classification tools (such as Kraken, Clark, USEARCH, ...).

The output of MetaTax is a file where each read is given a classification according to MetaTax voting
algorithm. This algorithm will try to classify reads at the lowest level possible, always using the 
majority voting criteria. 

Details about MetaTax algorithm will be released soon. Tests with several datasets has shown that 
at the species level MetaTax achieves a higher sensitivity and F1 score compared to individual 
classification tools and to WEVOTE, an existing voting algorithm based on LCA. 
Precision of MetaTax was higher that Clark-S, though it was lower than WEVOTE.

If you want to test MetaTax, create a file with the following structure:

        reads1_file_path (could be relative or absolute)
        N               (number of classifiers)
        classifier1_name
        classifier1_output_file_path
        ...
        classifierN_name
        classifierN_output_file_path

        ... and so on, for any number of reads files

Currently MetaTax supports output from the following taxonomic classifiers:
* Clark
* Clark-S
* Kraken
* OneCodex
* USEARCH/Blastn\*

\* # Note for USEARCH/Blastn #

It is important to note that for USEARCH/Blastn, the supported format is ``blast6out''. Furthermore, before
calling MetaTax with an output from USEARCH/Blastn, the output must be formatted to be readable. 

To do so, please run the script *process_usearch.py*. It uses the best-match criteria to select the taxonID
for each read from the USEARCH/Blastn output. You can configure the values used to detect the best hit by
editing this script. In the very beginning it defines these values:
* maxE (defaul 1e-4): maximum E-value to consider a hit (i.e. hits with E-value above this limit will not be considered)
* minIdentity (default 80): min percent of identity to consider a hit
* minSizeFactor (default 0.5): min percent of aligment length relative to the mean reads size
* readSize: put here the mean length of your reads

When 2 hits for the same read are considered, we choose the one with the highest identity. If both have the same identity, 
then we choose the one with the lowest E-value. Please feel free to implement your own ``best-hit'' approach by modifying 
the script *process_usearch.py*.

Doubts or comments: (diaztula@ime.usp.br)[mailto:diaztula@ime.usp.br]
