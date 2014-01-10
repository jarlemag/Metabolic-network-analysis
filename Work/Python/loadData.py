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

def ExpFluxesfromXML(filename,publication_id,reactor_id,experiment_id):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/reactor[@id="{react_id}"]/experiment[@id="{exp_id}"]/reactiondata'.format(pub_id = publication_id, react_id = reactor_id, exp_id = experiment_id)
    print(string)
    reactiondata = root.findall(string)
    reactionlist = reactiondata[0].findall('reaction')
    expfluxlist = [reaction.get('flux') for reaction in reactionlist]
    expfluxdict = {reaction.get('id'):reaction.get('flux') for reaction in reactionlist}
    return expfluxdict

expfluxes = ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')


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
            exprxn_coef = element.get('coef')
            exprxntuple = (exprxn_id,exprxn_coef)
            exprxnlist.append(exprxntuple)
        for element in modrxn_ids:
            modrxn_id = element.get('id')
            modrxn_coef = element.get('coef')
            modrxntuple = (modrxn_id,modrxn_coef)
            modrxnlist.append(modrxntuple)
        linkdict['expid'] =exprxnlist
        linkdict['modit'] =modrxnlist
        reactionmap.append(linkdict)
    return reactionmap
    
#rmap = ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')
