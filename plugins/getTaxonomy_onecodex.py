import sys, os, csv, pickle, time
sys.path.append("../")
from getTaxonomyFromEte3 import getTaxonomy as gt

def getTaxonomy(srcFile):
    return gt(srcFile, readNameColumn = 0, taxIdColumn = 1, skipHeader = True)
