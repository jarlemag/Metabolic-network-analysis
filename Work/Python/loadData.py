#loadExpDatafromXML.py

import xml.etree.ElementTree as ET

#simple commands:
#print(root.tag)
#print(root.attrib)
#for child in root:
#print(child.tag,child.attrib)
#root.attrib.keys()
#for element in root:
#    print(element)


class ReactionMap:
    def __init__(self):
        pass


class SplitMap:
    def __init__(self):
        pass


class SplitDefinition:
    def __init__(self):
        pass

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


def ExpSplitsfromXML(filename,publication_id,reactor_id,experiment_id):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/reactor[@id="{react_id}"]/experiment[@id="{exp_id}"]/splitsdata'.format(pub_id = publication_id, react_id = reactor_id, exp_id = experiment_id)
    splitsdata = root.findall(string)
    splitslist = splitsdata[0].findall('split')

    expsplitsdict = {split.get('id'):float(split.get('flux')) for split in splitslist}
    return expsplitsdict
    pass

def ExpSplitsfromTXT(source,destination,publication_id,reactor_id,experiment_id):
    pass


def ExpSplitsToXML(splits,filename):
    pass

def ExpErrorsfromXML(filename,publication_id,reactor_id,experiment_id):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/publication[@id="{pub_id}"]/reactor[@id="{react_id}"]/experiment[@id="{exp_id}"]/reactiondata'.format(pub_id = publication_id, react_id = reactor_id, exp_id = experiment_id)
    reactiondata = root.findall(string)
    reactionlist = reactiondata[0].findall('reaction')
    exp_errordict = {reaction.get('id'):float(reaction.get('error')) for reaction in reactionlist}
    return exp_errordict
    


def SplitMapfromXML(filename,publication_id,model_id):
    tree=ET.parse(filename)
    root=tree.getroot()
    string = './/publication[@id="{pub_id}"]/model[@id="{mod_id}"]/split'.format(pub_id=publication_id,mod_id=model_id)
    rawmap = root.findall(string)
    splitmap = {}
    for element in  rawmap:
        split_id = element.get('id')
        splitdict = {}
        numdict = {}
        denomdict = {}
        numreactions = element.findall('./numerator')[0].findall('./modelrxn')
        denomreactions = element.findall('./denominator')[0].findall('./modelrxn')
        numlist = []
        denomlist = []
        for reaction in numreactions:
            numdict = {}
            for attrib in reaction.attrib:
                numdict[attrib] = reaction.get(attrib)
            numlist.append(numdict)
        for reaction in denomreactions:
            denomdict ={}
            for attrib in reaction.attrib:
                denomdict[attrib] = reaction.get(attrib)
            denomlist.append(denomdict)
        splitdict['numerator'] = numlist
        splitdict['denominator'] = denomlist
        splitmap[split_id] =splitdict
    return splitmap
        
        
    
    return splits

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


class ReactionMap:
    def __init__(self,linkdicts):
        self.links = linkdicts




class ExpData:
    def __init__(self,expfluxdict):
        pass

def dictmapToList(dictmap,cobramodel):
    pass


class SBMLModel:
    def __init__(self,metabolites,reactions):
        self.metabolites = []
        self.reactions = []

    def add_metabolites(self,metabolitelist):
        #Check if metabolites are already in model.
        self.metabolites += metabolitelist #Add metabolites to model
        


class CMetabolite:
    def __init__(self,ident,name,compartment):
        self.id = ident
        self.name = name
        self.compartment = compartment
    def get_reactions(self):
        pass


def ModelfromSBML(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    string = './/listOfSpecies'
    species = root.findall(string)

    metabolitelist = []
    pass

def XMLtoGurobi(filename):
    #Create a Gurobi model directly from an FBA model saved as SBML
    #pseudocode:
    '''
    Load XML file as ElementTree
    Grab node "list of species"
    Make an SBMLmodel class instance.
    (To be continued)
    '''
    pass

if __name__ == "__main__":

    tree = ET.parse('expdata.xml')
    root = tree.getroot()
    
    rmap = ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')

    expfluxes = ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

    expsplits = ExpSplitsfromXML('expdata.xml','Schuetz','Batch','aerobe') #Should create a publication "Schuetz" and move splits there.

    splitmap = SplitMapfromXML('reactionmaps.xml','SCHUETZR','SCHUETZR')





