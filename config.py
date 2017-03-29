"""Allowed ranks are those of interest for our study. Feel free to change this list 
according to your needs"""
allowedRank = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

"""This dictionary tell us the name of available tools and their corresponding plugin in the 'plugins' folder. 
If you decide to include more tools, please implement their plugins and put them here"""
plugins = {'clark': 'getTaxonomy_clark', 'clarkS': 'getTaxonomy_clarkS', 'kraken': 'getTaxonomy_kraken', 'usearch': 'getTaxonomy_usearch', 'onecodex': 'getTaxonomy_onecodex', 'centrifuge': 'getTaxonomy_centrifuge'}
