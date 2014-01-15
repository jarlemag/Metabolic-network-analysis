#loadExpDatafromXML.py

import xml.etree.ElementTree as ET

#from lxml import etree

tree = ET.parse('expdata.xml')
root = tree.getroot()
#simple commands:
#print(root.tag)
#print(root.attrib)
#for child in root:
#print(child.tag,child.attrib)
#root.attrib.keys()
#for element in root:
#    print(element)

def ExpFluxesfromXML(filename,publication_id,reactor_id,experiment_id,vector = False):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/reactor[@id="{react_id}"]/experiment[@id="{exp_id}"]/reactiondata'.format(pub_id = publication_id, react_id = reactor_id, exp_id = experiment_id)
    reactiondata = root.findall(string)
    reactionlist = reactiondata[0].findall('reaction')
    if vector == True:
        expfluxlist = [reaction.get('flux') for reaction in reactionlist]
        return expfluxlist
    else:
        expfluxdict = {reaction.get('id'):float(reaction.get('flux')) for reaction in reactionlist}
    return expfluxdict

expfluxes = ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

def ExpErrorsfromXML(filename,publication_id,reactor_id,experiment_id):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/reactor[@id="{react_id}"]/experiment[@id="{exp_id}"]/reactiondata'.format(pub_id = publication_id, react_id = reactor_id, exp_id = experiment_id)
    reactiondata = root.findall(string)
    reactionlist = reactiondata[0].findall('reaction')
    exp_errordict = {reaction.get('id'):float(reaction.get('error')) for reaction in reactionlist}
    return exp_errordict
    


def ReactionMapfromXML(filename,publication_id,model_id):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/model[@id="{mod_id}"]/link'.format(pub_id=publication_id,mod_id=model_id)
    rawmap = root.findall(string)
    reactionmap = []
    for link in rawmap:
        linkdict = {}
        exprxn_ids = link.findall('./exprxn')
        modrxn_ids = link.findall('./modelrxn')
        exprxnlist = []
        modrxnlist = []
        for element in exprxn_ids:
            exprxn_id = element.get('id')
            if 'coef' in element.attrib:
                exprxn_coef = element.get('coef')
            else:
                exprxn_coef = 1
            exprxndict = {'rxid':exprxn_id,'coef':exprxn_coef}
            exprxnlist.append(exprxndict)
        for element in modrxn_ids:
            modrxn_id = element.get('id')
            if 'coef' in element.attrib:
                modrxn_coef = element.get('coef')
            else:
                modrxn_coef = 1
            modrxndict = {'rxid':modrxn_id,'coef':modrxn_coef}
            modrxnlist.append(modrxndict)
        linkdict['exprxns'] =exprxnlist
        linkdict['modrxns'] =modrxnlist
        reactionmap.append(linkdict)
    return reactionmap


def dictmapToList(dictmap,cobramodel):
    pass

if __name__ == "__main__":
    
    rmap = ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')





