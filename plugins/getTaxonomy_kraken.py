import sys, os, csv, pickle, time
sys.path.append("../")
from getTaxonomyFromEte3 import getTaxonomy as gt

def getTaxonomy(srcFile):
    return gt(srcFile, readNameColumn = 1, taxIdColumn = 2, skipHeader = False)
