#createReactionmap.py


import xml.etree.ElementTree as ET

def createReactionMap(filename,publication_id,model_id,reactionmap = None, debug = False):
    '''
    Writes a reaction map to an .xml file
    '''
    root = ET.Element('xml')
    publication = ET.SubElement(root,'publication',attrib = {'id':publication_id})
    model = ET.SubElement(publication,'model',attrib={'id':model_id})
    if debug == True:
        ET.dump(root)

    if reactionmap is not None:
        for linkdict in reactionmap:
            newlink = ET.SubElement(model,'link')
            for expreactiondict in linkdict['exprxns']:
                new_exprxn = ET.SubElement(newlink,'exprxn',attrib = {'id':expreactiondict['rxid'],'coef':str(expreactiondict['coef'])})
            for modelreactiondict in linkdict['modrxns']:
                if 'coef' in modelreactiondict:
                    new_modrxn = ET.SubElement(newlink,'modelrxn',attrib = {'id':modelreactiondict['rxid'],'coef':str(modelreactiondict['coef'])})
                else:
                    new_modrxn = ET.SubElement(newlink,'modelrxn',attrib = {'id':modelreactiondict['rxid']})
            
    tree = ET.ElementTree(root)
    #tree.write(str(filename))
    tree.write(filename)
    return


def createEmptyMap(publication_id,model_id):
    root = ET.element('xml')
    publication = ET.SubElement(root,'publication',attrib = {'id':publication_id})
    model = ET.SubElement(publication,'model',attrib={'id':model_id})
    tree = ET.ElementTree(root)
    return tree

def createSimpleReactionMap(simple_rmapdict,mapname,filename):
    pass

if __name__ == "__main__":

    publication_id = 'Perrenoud'
    model_id = 'SCHUETZR'
    filename = "testmap.xml"
    createReactionMap(publication_id,model_id,filename)

def addLinksToReactionMap():
    pass



def writeReactionMaptoFile(element_tree,filename):
    import xml.etree.ElementTree as ET
    element_tree.write(str(filename))

def reactionMapWizard():
    pass


def createLinkDict():
    pass



if __name__ == "__main__":
    import loadData as load

    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')

    createReactionMap('rewrittenmap.xml','Perrenoud','SCHUETZR',reactionmap = rmap)
