#addtodict.py
'''
Scans two dictionaries, A (targetdict) and B (defaultdict).
If an element in B is not in A, the element is added to A.
'''
def addtodict(targetdict,defaultdict):
    tempdict = targetdict
    for key in defaultdict:
        if (key in targetdict) == False:
            tempdict[key] = defaultdict[key]
    return tempdict
