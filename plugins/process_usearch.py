import sys, os, subprocess, csv
from bisect import bisect_left
import sqlite3
from PyEntrezId import Conversion

minE          = 1e-4
minIdentity   = 80
minSizeFactor = 0.5
readSize      = 100

class Hit(object):
    def __init__(self, rname, ident, tam, e, taxid):
        self.rname = rname
        self.ident = ident
        self.tam   = tam
        self.e     = e
        self.taxid = taxid
        
def matchCriteria(hit):
    if hit.taxid != -1 and hit.ident >= minIdentity and hit.e <= minE:
        return True
    return False

def getTaxId(an):
    sql = 'SELECT taxid FROM gi_taxid WHERE gi=%s'%an
    c.execute(sql)
    all_rows = c.fetchall()
    if len(all_rows) == 0:
        return -1
    return all_rows[0][0]

def writeFile(fout, hit):
    fout.write("%s\t%s\t%f\t%f\n"%(hit.rname, hit.taxid, hit.ident, hit.e))

def processUsearch(srcFname, outFname):
    ans  = []
    fin  = open(srcFname, "rt")
    fout = open(outFname, "wt")
    bestHit = Hit("", 0, 0, 0.0, "")
    conv = Conversion("diaztula@iq.usp.br")
    if not os.path.isfile("acc_2_taxid.tsv"):
        # First collect all accesion numbers
        accNumbers = []
        for l in fin:
            tupla = l.split("\t")
            if len(tupla) == 12:
                rname = tupla[0].split("/")[0].replace("@", "")
                ident = float(tupla[2])
                e     = float(tupla[10])
                tam   = int(tupla[3])
                an    = tupla[1].split("|")[3].split(".")[0].strip()
                accNumbers.append(an)
        accNumbers = list(set(accNumbers))
        print len(accNumbers)
        d = dict()
        for acc in accNumbers:
            taxid = conv.convert_accession_to_taxid(acc)
            if taxid > 0:
                d[acc] = taxid
        tmpFile = open("acc_2_taxid.tsv", "wt")
        for k in d.keys():
            tmpFile.write("%s\t%s\n"%(k, d[k]))
        return
    else:
        print "carregando dicionario"
        d = dict()
        for l in csv.reader(open("acc_2_taxid.tsv", "rt"), delimiter="\t"):
            d[l[0]] = l[1]
        for l in fin:
            tupla = l.split("\t")                                                                                                                                                                           
            if len(tupla) == 12:
                rname = tupla[0].split("/")[0].replace("@", "")
                ident = float(tupla[2])
                e     = float(tupla[10])
                tam   = int(tupla[3])
                an    = tupla[1].split("|")[3].split(".")[0].strip()
                try:
                    taxid = d[an]
                except:
                    taxid = conv.convert_accession_to_taxid(an)
                    d[an] = taxid

                hit = Hit(rname, ident, tam, e, taxid)
                if taxid == -1:
                    print an
                if taxid == -1:
                    print "No taxid for accession number |%s|"%an, " read = ", hit.rname
                if bestHit.rname == "" and matchCriteria(hit):
                    bestHit = hit
                elif bestHit.rname != hit.rname:
                    if bestHit.rname != "":
                        writeFile(fout, bestHit)
                    if matchCriteria(hit):
                        bestHit = hit
                    else:
                        bestHit = Hit("", 0, 0, 0.0, "")
                else:
                    if hit.ident > bestHit.ident or (hit.ident == bestHit.ident and hit.e < bestHit.e):
                        bestHit = hit
        tmpFile = open("acc_2_taxid.tsv", "wt")
        for k in d.keys():
            tmpFile.write("%s\t%s\n"%(k, d[k]))
    
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: python %s srcFileName outFileName"%sys.argv[0]
        sys.exit(1)
    processUsearch(sys.argv[1], sys.argv[2])
