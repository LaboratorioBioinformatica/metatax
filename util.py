from collections import OrderedDict

# ----------------------------------------------
def isNA(lineage):
    """Receives a lineage and returns True if it is NA.
    """
    if "root" in lineage.keys() and lineage["root"] == "0":
        return True
    return False
