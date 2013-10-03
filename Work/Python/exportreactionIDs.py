

def exportreactionIDs(model,filename):
    out = []
    for reaction in model.reactions:
        out.append(str(reaction))
    out = "\n".join(out) 
    writeout = open(filename, "a")
    writeout.write(out)
    writeout.close()
    return out
