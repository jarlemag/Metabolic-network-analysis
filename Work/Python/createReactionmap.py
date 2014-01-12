#createReactionmap.py


import xml.etree.ElementTree as ET

def createReactionMap(publication_id,model_id,filename, debug = False):
    '''
    Writes a reaction map to an .xml file
    '''
    root = ET.Element('xml')
    publication = ET.SubElement(root,'publication',attrib = {'id':publication_id})
    model = ET.SubElement(publication,'model',attrib={'id':model_id})
    if debug == True:
        ET.dump(root)
    tree = ET.ElementTree(root)

    tree.write(str(filename))
    return


def createSimpleReactionMap(simple_rmapdict,mapname,filename):
    pass

if __name__ == "__main__":

    publication_id = 'Perrenoud'
    model_id = 'SCHUETZR'
    filename = "testmap.xml"
    createReactionMap(publication_id,model_id,filename)
