import sys, os, csv, pickle, time
from getTaxonomyFromEte3 import getTaxonomy as gt

def getTaxonomy(srcFile):
    return gt(srcFile, readNameColumn = 0, taxIdColumn = 2, skipHeader = False, sep=",")
