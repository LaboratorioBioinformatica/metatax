import sys, os, csv, importlib, argparse, traceback, numpy as np
from os.path import join
from Bio import SeqIO
from ete3 import NCBITaxa
from util import *
from config import *
import multiLevelVoting

myPath = os.path.split( os.path.abspath(__file__) )[0]
sys.path.append(join(myPath, "plugins"))
from getTaxonomyFromEte3 import *


def loadInput(inputFile):
    """Parses input file to load information about reads files, and classifiers output that 
    will be used for metaclassification.

    Parameters:
    - inputFile: text file with the following format:
    	reads1_file_path (could be relative or absolute)
    	N               (number of classifiers)
    	classifier1_name
    	classifier1_output_file_path
    	...
    	classifierN_name
    	classifierN_output_file_path
    	... and so on, for any number of reads files
    Returns:
    - a dictionary with data about the reads. Keys are the reads file names.
    """
    reads        = {}
    reader = open(inputFile, "rt")
    while True:
        try:
            readName = reader.next().strip()
            reads[readName] = []
            n        = int(reader.next().strip())
            for i in range(n):
                classif = {}
                classif['classifName'] = reader.next().strip()
                if classif['classifName'] in plugins.keys():
                    classif['classifData' ] = reader.next().strip()
                    classif['module']       = importlib.import_module(plugins[classif['classifName']])
                else:
                    print "I do not understand output from %s, exiting"%classif['classifName']
                    sys.exit(1)
                reads[readName].append(classif)
        except:
            break
    return reads

# ------------------------------------------------------------------------------- #
# MAIN FUNCTION
# ------------------------------------------------------------------------------- #
if __name__ == "__main__":
    """Main function of the method. It receives an input file with data about reads and their classification 
    by different tools, and produces an output dir with meta-classification.
    """
    # Parsing arguments
    argp = argparse.ArgumentParser()
    argp.add_argument('i', help = 'Input file with information about read files and classification results')
    argp.add_argument('o', help = 'Dir to write the results of meta-classification')
    argp.add_argument('-log', help = 'Write log information (this could create a very big file! (default = 0)', required = False, action="store_true")
    args = argp.parse_args()

    inputFile    = args.i
    outDir       = args.o
    votingMethod = "multilevel" # Future extensions should allow for different voting methods, including the WEVOTE method
    if not os.path.isdir(outDir):
        try:
            os.makedirs(outDir)
        except: 
            traceback.print_exc()
            sys.exit(1)
    LOG = args.log
    writeFullLineage = True
    ncbi = NCBITaxa()
    if LOG:
        logF      = open(join(outDir, "log.txt") , "wt")
    else:
        logF      = open(os.devnull, "w")

    # Get reads information, including reads file name and data from different classifiers
    readsDict = loadInput(inputFile)

    # Analizing each read file and its corresponding classification data
    for readName in readsDict.keys():
        logF.write("<=========================================================>\n")
        logF.write("Analyzing read file: %s\n"%readName)
        logF.write("Classifiers:\n")
        classifiers  = []
        classifNames = []
        for c in readsDict[readName]:
            logF.write("\t%s\t->\t%s\n"%(c['classifName'], c['classifData']))
            classifNames.append(c['classifName'])
            fullLineage = c['module'].getTaxonomy(c['classifData'])
            classifiers.append(fullLineage)

        # Creating output files
        prefix = readName[readName.rfind('/')+1:]
        classifiedF   = open(join(outDir, prefix+"_classified.tsv")  , "wt")
        disagreeF     = open(join(outDir, prefix+"_disagreement.tsv"), "wt")
        naF           = open(join(outDir, prefix+"_NAs.tsv")         , "wt")
        if writeFullLineage:
            classifiedF.write("Read\tToolsN\tTotalClassif\tVotes\tPercentVotes\tFullLineage\n")
        else:
            classifiedF.write("Read\tToolsN\tTotalClassif\tVotes\tPercentVotes\tLineage\n")
        if readName.endswith(".fa") or readName.endswith(".fasta"):
            ext = "fasta"
        else:
            ext = "fastq"
        reads = SeqIO.parse( open(readName, "rU"), ext )
        # For each read
        for r in reads:
            if votingMethod == "multilevel":
                classTree = multiLevelVoting.ClassTree()
            lineageList = []
            usedClassif = []
            rname = r.name.split("/")[0]
            logF.write("<==================================================================>\n")
            logF.write("Classifying read \'%s\'\n"%rname)
            for i, classifier in enumerate(classifiers): # Get classification from each classifier
                try:
                    lineage = classifier[rname]
                    if not isNA(lineage):
                        lineageList.append(lineage)
                        usedClassif.append(classifNames[i])
                except: # Read was not classified
                    pass
            if len(lineageList) > 0: # At least one classification not NA
                # Add each lineage to the classification algorithm
                for i, lineage in enumerate(lineageList):
                    classTree.addClassification(lineage)
                    logF.write("<------------------------------\n")
                    logF.write("Classification according to %s:\n"%(usedClassif[i]))
                    logF.write("%s\n"%(",".join( ["%s:%s"%(k, lineage[k]) for k in lineage.keys()] )))
                    logF.write("------------------------------>\n")
                rank, tid, votes, completeLin = classTree.getClassification()
                logF.write("Final classification by method \'%s\': %s: %s with %i votes\n"%(votingMethod, rank, tid, votes))
                if rank == "NA":
                    naF.write("%s\n"%rname)
                elif rank == "disagreement":
                    disagreeF.write("%s\n"%rname)
                else:
                    if writeFullLineage:
                        ll = list(completeLin.values())
                        names = ncbi.get_taxid_translator(ll)
                        stringList = []
                        for i in range(len(completeLin)):
                            key_i = completeLin.keys()[i]
                            stringList.append("%s|%s|%s"%(key_i, names[ completeLin[key_i] ], completeLin[key_i]))
                        classifiedF.write("%s\t%i\t%i\t%i\t%i\t%s\n"%(rname, len(classifiers), len(lineageList), votes, int(votes*100/len(classifiers)), "\t".join(stringList)))
                    else:
                        name = ncbi.get_taxid_translator([tid])
                        name = name[int(tid)]
                        classifiedF.write("%s\t%i\t%i\t%i\t%i\t%s\n"%(rname, len(lineageList), votes, int(votes*100/len(classifiers)), "|".join([rank, name, tid])))
            else:
                naF.write("%s\n"%rname)
                logF.write("Read without classification: NA\n")
            logF.write("<==================================================================>\n")
        logF.write("<=========================================================>\n")

