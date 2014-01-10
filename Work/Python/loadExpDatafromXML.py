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

def loadExpFluxesfromXML(filename,publication_id,reactor_id,experiment_id):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/reactor[@id="{react_id}"]/experiment[@id="{exp_id}"]/reactiondata'.format(pub_id = publication_id, react_id = reactor_id, exp_id = experiment_id)
    print(string)
    reactiondata = root.findall(string)
    reactionlist = reactiondata[0].findall('reaction')
    expfluxlist = [reaction.get('flux') for reaction in reactionlist]
    expfluxdict = {reaction.get('id'):reaction.get('flux') for reaction in reactionlist}
    return expfluxdict


expfluxes = loadExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')


