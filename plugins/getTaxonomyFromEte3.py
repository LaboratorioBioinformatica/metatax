import sys, os, csv
from ete3 import NCBITaxa
from os.path import join, isfile
from collections import OrderedDict
sys.path.append("../")
import config

def getLineageDict(originalId, ncbi, allowedRank):
    '''
    getLineageDict: returns the full lineage of a taxon ID as an ordered dictionary

    Parameters:
    - originalId: taxon id
    - ncbi: an instance of NCBITaxa() class from ete3 package (to avoid creating this object frequently, it is
    expected that the caller creates one instance of NCBITaxa() and pass it repeatedly to this function).
    - allowedRank: a list of rank names (e.g. ["phylum", ...]) to be considered. Useful to filter out ranks such as 
    "subclass", "strain" etc...

    Returns: an OrderedDict() instance with the full lineage. 

    Remarks: I just found that some taxon ids have been merged and ete3 does not handle them well, e.g 269483 was 
    merged with 482957, but get_lineage(269483) returns [1]. The workaround is to translate taxid to name
    and then back to taxid.
    '''
    nameList = ncbi.get_taxid_translator([originalId])
    if len(nameList) == 0:
        return OrderedDict([("root", "0")])
    name = nameList[originalId]
    idL = ncbi.get_name_translator([name])[name]
    if originalId in idL:
        id = originalId
    else:
        id = idL[-1]
    #if id != originalId:
    #    print originalId, id, name
    lineage = ncbi.get_lineage(id)
    if len(lineage) == 1 and ncbi.get_taxid_translator(lineage)[lineage[0]] == "root": # Only root
        lineageDict = OrderedDict([("root", "0")])
    else:
        lineageDict = OrderedDict()
        for lin in lineage:
            if lin != 0:
                name = ncbi.get_taxid_translator([lin])[lin]
                rank = ncbi.get_rank([lin])[lin]
                if rank in allowedRank:
                    lineageDict[rank] = lin
        if len(lineageDict) == 0:
            lineageDict['root'] = '0'
    return lineageDict


def getTaxonomy(srcFile, readNameColumn, taxIdColumn, skipHeader = True, sep = "\t", secondOption = -1, skipLines = False, log = False):
    '''getTaxonomy: parses a taxonomy file (output from some classifier) and returns a dictionary with the full lineage for each read. 
    
    Paramenters:
    - srcFile: file with the taxonomy classification
    - readNameColumn: column index (zero-based) of the read name
    - taxIdColumn: column index (zer0-based) of the taxon ID
    - skipHeader: whether the file contains a header (true) or not
    - sep: symbol used as separator in the input file (tab, comma, semicolon, etc...)
    - secondOption: index of a second taxon ID to consider. This was included due to Clark-S 
    returning 2 classifications. In the first is NA, then we consider the second (default = -1, no 
    second option available)
    - skipLines: some classifiers output 2 lines for each read (such as MyTaxa). 
    Set this to true if it is the case.
    - log: set to True to see messages in the terminal output

    Returns: a dictionary where read names are the keys, and the values are OrderedDict objects with the full lineage.

    Remarks: This function calls "getLineageDict", which receives a list of allowed ranks.
    '''
    # Validating directory
    if not os.path.isfile(srcFile):
        print "Source file %s not found"%srcFile
        return {}
    ncbi = NCBITaxa()
    taxIds   = []
    if log: print "Reading entries from file : ", f
    fileIn = csv.reader(open(srcFile, "rt"), delimiter=sep)
    if skipHeader:
        fileIn.next()
    skipThisLine = False
    for i, l in enumerate(fileIn):
        if skipLines and skipThisLine:
            skipThisLine = not skipThisLine
            continue
        skipThisLine = not skipThisLine
        taxId = 0
        try:
            taxId = int(l[taxIdColumn])
        except:
            if secondOption != -1:
                try:
                    taxId = int(l[secondOption])
                except:
                    continue
            else:
                continue
        if taxId != 0:
            taxIds.append(taxId)
    if log: print "Total number of entries   : ", len(taxIds)
    singletonList = list(set(taxIds))
    if log: print "Number of unique IDs      : ", len(singletonList)
    lineageDict = dict()
    for i, id in enumerate(singletonList):
        # getting complete lineage
        lineageDict[id] = getLineageDict(id, ncbi, allowedRank = config.allowedRank)
    if log: print "Number of lineages created: ", len(lineageDict.keys())
    # =============================================
    # Populate file entries
    if log: print "Number of lineages loaded: ", len( lineageDict.keys() )
    # Create table with lineage for each read
    if log: print "getting full lineage for %s"%(srcFile)
    fileIn = csv.reader(open(srcFile, "rt"), delimiter=sep)
    if skipHeader:
        fileIn.next()

    outDict = {}
    # fout   = open(os.path.join(outDir, f+"_with_lineage.tsv"), "wt")
    skipThisLine = False
    conut = 0
    for l in fileIn:
        if skipLines and skipThisLine:
            skipThisLine = not skipThisLine
            continue
        skipThisLine = not skipThisLine
        try:
            taxId = int(l[taxIdColumn])
        except:
            if secondOption != -1:
                try:
                    taxId = int(l[secondOption])
                except:
                    continue
            else:
                continue
        if taxId != 0:
            lin = lineageDict[taxId]
        else:
            lin = OrderedDict([("root", "0")])
        cleanReadName = l[readNameColumn].split(" ")[0].replace("@", "").replace(">", "").strip()
        # fout.write("%s\t%s\n"%(l[readNameColumn].split(" ")[0].replace("@", "").replace(">", "").strip(), lin))
        outDict[cleanReadName] = lin
    return outDict

if __name__ == "__main__":
    print getTaxonomy("../../../kraken_db_test/results_classifiers/clark/HiSeq_accuracy_classif_clark.csv", readNameColumn=0, taxIdColumn=2, skipHeader = False, sep = ",")
