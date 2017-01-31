import numpy as np, sys
from util import *
import config

EQUAL        = 1
COMPATIBLE   = 2
INCOMPATIBLE = 3

# ----------------------------------------------
class Node(object):
    """Represents a single branch of the taxon ID tree. 
    Args:
	lineage: Ordered dict with the lineage with the format: 'rank: taxonId'.
    	equalC is the number of branches equal to this one.
    	compatC is the number of branches compatible to this one.
    	Note that equalC+compatC is the number of "votes" this branch received, while
	    running the meta-classification.
    """
    def __init__(self, lineage):
        self.lineage = lineage
        self.equalC  = 1
        self.compatC = 0

# ----------------------------------------------
class ClassTree(object):
    """Class for running the Metatax method.
    """
    def __init__(self):
        self.roots    = []
        self.lineages = []


    def compare(self, lineage1, lineage2):
        """Compares two lineages, assuming that the first one has the lowest level 
        (though both lineages can have the lowest level).
        Args:
            lineage1: one lineage having the lowest level.
            lineage2: the other lineage.
        Returns: EQUAL, COMPATIBLE, or INCOMPATIBLE
        """
        if lineage1.keys()[-1] == lineage2.keys()[-1] and lineage1[lineage1.keys()[-1]] == lineage2[lineage2.keys()[-1]]:
            return EQUAL
        big    = lineage1
        bigK   = big.keys()
        small  = lineage2
        smallK = small.keys()
        for k in smallK:
            if not k in bigK: return INCOMPATIBLE
            if big[k] != small[k]: return INCOMPATIBLE
        return COMPATIBLE


    def addClassification(self, lineage):
        """Append a lineage (branch) to the classification algorithm
        Args:
            lineage: branch to be added.
        """
        self.lineages.append(lineage)

    def __addClassification(self, lineage):
        """Adds a new lineage to the classification algorithm. The new branch is compared 
        with all existing branches and, if equal or compatible, the counter of the representative
        branch is increased. It it is not compatible with any branch, it becomes a new representative
        of itself.
        Args:
            lineage: new branch to be added to the algorithm.
        """
        if not isNA(lineage): # Discard NA classification
            if len(self.roots) == 0: # If it is the first lineage, make it a root
                self.roots.append(Node(lineage))
            else:
                # Find if lineage is equal to any axisting root or the largest root compatible with it
                compat = []
                equal  = False
                for i, node in enumerate(self.roots):
                    result = self.compare(node.lineage, lineage)
                    if result == EQUAL:
                        equal = True
                        node.equalC += 1
                        break
                    elif result == COMPATIBLE:
                       compat.append(node)
                if len(compat) == 0:
                    self.roots.append(Node(lineage))
                else:
                    for c in compat:
                        c.compatC += 1

    def getLowestLevel(self, lineages, levels):
        """Lazy method to find the lowest level (e.g. species, genus, ...) present it at least
        one lineage. Assumes that levels are from high to low in the 'levels' list.
        Args: 
            lineages: a list of branches.
            levels: list of leves, sorted by increasing order of levels.
        Returns:
            The lowest level it finds, or None if there is no level.
        """
        for i, rank in enumerate(reversed(levels)):
            for l in lineages:
                if rank in l.keys():
                    return rank, i
        return None, 0

    def pruneLowestLevel(self, lineages):
        """Prune all branches by removing the lowest level.
        Args:
           lineages: list of branches.
        Returns:
            a new list of branches with the lowest level pruned.
        """
        lowestLevel, index = self.getLowestLevel(lineages, config.allowedRank)
        if lowestLevel is None:
            return []
        newLineages = []
        for i in range(len(lineages)):
            l = [(k, v) for k, v in zip(lineages[i].keys(), lineages[i].values()) if k != lowestLevel]
            if len(l) > 0:
                newLineages.append(OrderedDict(l))
        return newLineages


    def getClassification(self):
        """Returns a classification according to the roots. This algorithm is a multilevel algorithm that works as follows: \
        First it tries to classifiy using all levels from all lineages. If no concensous is obtained, it prunes tha lowest level of the \
        taxonomy and tries again. After several iterations, it should either comes up with a classification or prune the highest level and \
        return disagree.
        Return:
            Classification according to meta-tax: lowest rank, corresponding taxon ID, number of votes received
            and full lineage.
        """
        # Sort lineages, so that first ones have the lower levels. Note that lineages (branches) could have missing ranks.
        # Because of this, we first sort by number of ranks and then by levels. Because sort method are stable, we 
        # endup with the larger lineages with the lower levels first.
        # Sort by taxonomy tree length
        l = [len(c) for c in self.lineages]
        idxs = list(np.argsort(l))
        idxs.reverse()
        largerLineages = [self.lineages[i] for i in idxs]
        
        lowestRankList = []
        for l in largerLineages:
            dummy, index = self.getLowestLevel([l], config.allowedRank)
            lowestRankList.append(index)
        idxs = list(np.argsort(lowestRankList))
        sortedLineages = [largerLineages[i] for i in idxs]
        self.lineages = sortedLineages
        # Now lineages are sorted first by number of ranks and second by lowest level
        # This assumption is important for the forthcoming algorithn
        lineages = []
        for i in range(len(self.lineages)):
            l = [(k, v) for k, v in zip(self.lineages[i].keys(), self.lineages[i].values())]
            lineages.append(OrderedDict(l))
        levels = [r for r in reversed(self.lineages[0])]
        levelIdx = 0
        while True:
            self.roots    = []
            for l in lineages:
                self.__addClassification(l)
            rank, tid, votes, completeLin = self.__getClassification()
            if rank == "disagreement":
                # No classification, prune the last level if possible
                newLineages = self.pruneLowestLevel(lineages)
                if len(newLineages) > 0:
                    lineages = newLineages
                else:
                    return "disagreement", "disagreement", 0, None
            else:
                return rank, tid, votes, completeLin
    # --------------------------------------------

    def __getClassification(self):
        """Checks whether there is a winner branch. 

        Returns: If there is a winner branch, returns the lowest rank, 
        lowest taxonID, number of votes to the winner branch, and the full lineage
        as an ordered dict {rank: taxonID}.
        	If there is no winner branch, i.e. there are ties, returns 'disagreement'
        """
        if len(self.roots) == 0:
            return "NA", "0", 0, None
        elif len(self.roots) == 1:
            return self.roots[0].lineage.keys()[-1], self.roots[0].lineage.values()[-1], self.roots[0].equalC + self.roots[0].compatC, self.roots[0].lineage
        # First sort by the sum of equalC + compatC
        l = [r.equalC + r.compatC for r in self.roots]
        idxs = list(np.argsort(l))
        idxs.reverse()
        sortedRoots = [self.roots[i] for i in idxs]
        # Simplest scenario: no ties
        if l.count(l[idxs[0]]) == 1:
            return sortedRoots[0].lineage.keys()[-1], sortedRoots[0].lineage.values()[-1], sortedRoots[0].equalC + sortedRoots[0].compatC, sortedRoots[0].lineage
        else: # We have a tie, then return "disagreement" so the algorithm can prune the tree 
            return "disagreement", "disagreement", 0, None

if __name__ == "__main__":
    from collections import OrderedDict
    
    c1 = OrderedDict( [('superkigdom', 1), ('phylum', 2), ('class', 3), ('order', 4), ('family', 5)] )
    c2 = OrderedDict( [('superkigdom', 1), ('phylum', 2), ('class', 3), ('order', 5)] )
    c3 = OrderedDict( [('superkigdom', 1), ('phylum', 2), ('class', 6)] )
    l = [c1, c2, c3]
    """
    c1 = OrderedDict( [('superkigdom', 1), ('phylum', 2), ('class', 3), ('order', 4)] )
    c2 = OrderedDict( [('superkigdom', 1), ('phylum', 2)] )
    c3 = OrderedDict( [('superkigdom', 1), ('phylum', 5), ('class', 6)] )
    c4 = OrderedDict( [('superkigdom', 1), ('phylum', 5)] )
    l = [c1, c2, c3, c4]
    """
    classTree = ClassTree() # classification data structure
    for c in l:
        classTree.addClassification(c)
    print classTree.getClassification()

