#addtodict.py

def addtodict(targetdict,defaultdict)
    tempdict = targetdict
    for key in defaultdict:
        if (key in targetdict) == False:
            tempdict[key] = defaultdict[key]
    return targetdict
