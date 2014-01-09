#loadExpDatafromXML.py

import xml.etree.ElementTree as ET

#from lxml import etree

tree = ET.parse('experimentaldata.xml')
root = tree.getroot()

#print(root.tag)
#print(root.attrib)
#for child in root:
#    print(child.tag,child.attrib)
#root.attrib.keys()

#for element in root:
#    print(element)


def getNodesByType(tree,typeID):
    return root.findall(typeID) #Return all child nodes of type "publication"
    


def selectPublications(tree,identity):
    pubs = getNodesByType(tree,'publication')
    hits = []
    for pub in pubs:
        if pub.get('id') == identity: #If the the "id" attribute of the node equals the input
            hits.append(pub) #add the node to the list of hits
    return hits

hits = selectPublications(root,'Perrenoud')

#Simple format (attribute values only):
'''
['Perrenoud, 'Batch','Aerobe'] #
'''


#List with element names ("publication"), attribute names and attribute values:
'''
[Publication  id = "Perrenoud", reactor type = "Batch", experiment id = "aerobe"
'''

#List of lists including element name and attribute:value dictionaries
'''
[['Publication'],{id:"Perrenoud}],["reactor",{type:"Batch"}],["experiment",{id: "aerobe"}]]

'''

def findNodeByAttribute(root,attributes,targetvalues):
    '''
    Finds and returns nodes where the  value of the attribute defined by the argument "field" equals 
    that found in the argument "targets", a list of attribute values.
    '''
    print('findNodeByAttribute:')
    if (type(targetvalues) is str): #If the input target value is not a list (only one value), make it one.
        targetvalues = [targetvalues]
    newnode = root
    print('Starting search for node matching specified attribute values.')
    print('len(targetvaues)',len(targetvalues))
    for i in range(0,len(targetvalues)):
        print('i:',i)
        value = targetvalues[i]
        print('Current targetvalue:',value)
        success = False
        currentnode = newnode
        for j in range(0,len(currentnode)):
            print('j:',j)
            if currentnode[j].get('id') == targetvalues[i]:
                print('Found matching attribute value:',value)
                newnode = currentnode[j]
                print('New current node:',newnode.tag)
                success = True
                
        
        
            
        if success == False:
            print('Failed to find matching node.')
            return
    return newnode

#q = findByAttributeValue(root,'id','Perrenoud')

#w = findByAttributeValue(root,'id',','Batch')
targetvalues = ['Perrenoud','Batch']

z = findNodeByAttribute(root,'id',['Perrenoud','Batch','aerobe'])

print(z)
