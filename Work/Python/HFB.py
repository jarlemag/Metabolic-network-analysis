#HFB.py
from collections import defaultdict

def HFBreactions(cobramodel,fluxdict):
    print 'Finding HFB for model ',cobramodel.description
    maxproducerslist = []
    maxconsumerslist = []
    tol = 1e-6
    for metabolite in cobramodel.metabolites:
        reactions = metabolite._reaction
        maxproducer_ids = []
        maxproducer_flux = 0
        maxconsumer_ids = []
        maxconsumer_flux = 0
        for reaction in reactions:
            coef = reaction.get_coefficient(metabolite.id)
            massflux = fluxdict[reaction.id] * coef
            if abs(massflux) > tol: #If flux is non-zero
                if  massflux > maxproducer_flux:   
                    maxproducer_flux = massflux
                    maxproducer_ids = [reaction.id]
                elif massflux == maxproducer_flux:
                    maxproducer_ids.append(reaction.id)
                    
                elif massflux < maxconsumer_flux:
                    maxconsumer_flux = massflux
                    maxconsumer_ids = [reaction.id]
                elif massflux == maxconsumer_flux:
                    maxconsumer_ids.append(reaction.id)
    
            if len(maxproducer_ids) > 0:       
                maxproducerslist += maxproducer_ids
            if len(maxconsumer_ids) > 0:
                maxconsumerslist += maxconsumer_ids
    #print 'maxproducers:',maxproducerslist #DEBUG
    #print 'maxconsumers:',maxconsumerslist #DEBUG
    HFBrxset = set(maxproducerslist).intersection(maxconsumerslist)
   
    return HFBrxset


def listconsumers(metabolite):
    reactions = metabolite.get_reaction()
    

def listreactions(model,metabolite,fluxdict):
    if type(metabolite) == type(''):
        metabolite = model.metabolites.get_by_id(metabolite)
    reactions = metabolite.get_reaction()
    print '\n'
    print 'metabolite:',metabolite.id
    print '{:20s} {:20s} {:20s} {:20s}\n'.format('Reaction ID','Reaction rate','Coefficient','Mass flux')
    for reaction in reactions:
        reactionrate = fluxdict[reaction.id]
        coefficient = reaction.get_coefficient(metabolite.id)
        massflux = reactionrate * coefficient
        print '{:15s} {:15.4f} {:15.1f} {:15.4f}\n'.format(reaction.id,reactionrate,coefficient,massflux)  
    return


def listproducers(model,metabolite_id,fluxdict):
    if type(metabolite_id) == type(''):
        metabolite = model.metabolites.get_by_id(metabolite_id)
    else:
        metabolite = metabolite_id
    metabolite_reactions = metabolite.get_reaction()
    producersdict = {}
    for rx in metabolite_reactions:
        coef = rx.get_coefficient(metabolite_id)
        if coef > 0:
            producersdict[rx.id] = coef * fluxdict[rx.id]
    return producersdict

def reactionInfo(reaction):
    print 'reaction ID:',reaction.id
    print ''
    pass

def maxProducer(model,metabolite_id,fluxdict):
    if type(metabolite_id) == type(''):
        metabolite = model.metabolites.get_by_id(metabolite_id)
    else:
        metabolite = metabolite_id
    metabolite_reactions = metabolite.get_reaction()
    producersdict = {}
    for rx in metabolite_reactions:
        coef = rx.get_coefficient(metabolite_id)
        massflux = coef * fluxdict[rx.id]
        if massflux > 0:
            producersdict[rx.id] = massflux
    maxproducer = max(producersdict, key = producersdict.get)
    return maxproducer


def maxConsumer(model,metabolite_id,fluxdict): #Not giving correct results currently!
    if type(metabolite_id) == type(''):
        metabolite = model.metabolites.get_by_id(metabolite_id)
    else:
        metabolite = metabolite_id
    metabolite_reactions = metabolite.get_reaction()
    consumersdict = {}
    for rx in metabolite_reactions:
        coef = rx.get_coefficient(metabolite_id)
        massflux = coef*fluxdict[rx.id]
        if massflux < 0:
            consumersdict[rx.id] = massflux
    #print 'consumers:',consumersdict
    maxconsumer = min(consumersdict, key = consumersdict.get)
    return maxconsumer


def connectHFBmetabolites(model,fluxdict,reaction_id):
    if type(reaction_id) == type(''):
        reaction = model.reactions.get_by_id(reaction_id)
    else:
        reaction = reaction_id

    connections = defaultdict(list)
    #For every reactant in the reaction
    for reactant in reaction.get_reactants():
        maxconsumed = []
        maxproduced = []
        max_consumer = maxConsumer(model,reactant.id,fluxdict)
        if max_consumer == reaction.id:
            maxconsumed.append(reactant.id)
    #
    for product in reaction.get_products():
        max_producer = maxProducer(model,product.id,fluxdict)
        if max_producer == product.id:
            maxproduced.append(product.id)

    #print '\n'
    #print 'reaction:',reaction.id
    #print 'maxconsumed:',maxconsumed
    #print 'maxproduced:',maxproduced
    for element in maxconsumed:
        for product in maxproduced:
            connections[element].append(product)
    return connections
        
        

def HFBnetwork(cobramodel,fluxdict, grouped = False):
    HFBrxset = HFBreactions(cobramodel,fluxdict)
    #print 'len HFBrxset:',len(HFBrxset)
    #networkdict = {}
    networkdict = defaultdict(list)
    #For every reaction in the HFB
    for reaction_id in HFBrxset:
        #print '-------------'
        #print 'HFB reaction:',reaction_id #DEBUG
        reaction = cobramodel.reactions.get_by_id(reaction_id)

        maxconsumed = []
        maxproduced = []
        #For every reactant in the reaction
        for reactant in reaction.get_reactants():
            reactant_id = reactant.id
            #print 'reactant id:',reactant_id #DEBUG
            reactant = cobramodel.metabolites.get_by_id(reactant_id) 
            metabolite_reactions = reactant.get_reaction() #Get the reactions in which the reactant partakes
            consumersdict = {}
            #For every reaction in which the reactant partakes:
            for rx in metabolite_reactions:
                coef = rx.get_coefficient(reactant_id)
                massflux =  coef * fluxdict[rx.id]
                if massflux < 0: #If the reactant is being consumed in that reaction:
                    consumersdict[rx.id] = massflux #Record the rate of consumption
            maxconsumer =  min(consumersdict, key = consumersdict.get) #Find the highest-consuming reaction for the reactant
         #   print 'max consumer:',maxconsumer #DEBUG
            if maxconsumer == reaction_id: #If the highest-consuming reaction is the current reaction
                maxconsumed.append(reactant_id)
            
        #For every product in the reaction:
        for product in reaction.get_products(): 
            #product = cobramodel.metabolites.get_by_id(product_id)
            product_id = product.id
            metabolite_reactions = product.get_reaction() #Get the reactions in which the product partakes
            producersdict = {}
            #For every reaction in which the product partakes:
            for rx in metabolite_reactions:
                coef = rx.get_coefficient(product_id)
                massflux = coef * fluxdict[rx.id]
                if massflux > 0: #If the product is being produced in that reaction:
                    producersdict[rx.id] = massflux  #Record the rate of production
        maxproducer = max(producersdict, key = producersdict.get) #Find the highest-producing reaction for the product
        #print 'max producer:',maxproducer #DEBUG
        if maxproducer == reaction_id: #If the highst-producing reaction is the current reaction
            maxproduced.append(product_id)

        #print 'max consumed:',maxconsumed
        #print 'maxp produced:',maxproduced
        for metabolite in maxconsumed:
            for product in maxproduced:
                networkdict[metabolite].append(product)
    return networkdict

def HFBtoSIF(networkdict,filename):
    target = open(filename, 'w')
    for reactant in networkdict:
        #print '-----------'
        #print 'reactant:',reactant
        for product in networkdict[reactant]:
            #print 'product:',product
            target.write(reactant)
            target.write(' ')
            target.write('produces')
            target.write(' ')
            target.write(product)
            target.write('\n')
    target.close()
    return


def HFBdetails(model,HFB,fluxdict,hold = False):
    print '\n','HFB details:'
    count = 1
    
    for reaction_id in HFB:
        if hold:
            z = raw_input('Press enter to continue.')
        reaction = model.reactions.get_by_id(reaction_id)
        reactants = reaction.get_reactants()
        products = reaction.get_products()
        for reactant in reactants:
            print 'reactant:',reactant
            print 'maxconsumer:',maxConsumer(model,reactant,fluxdict)
            listreactions(model,reactant,fluxdict)
        for product in products:
            print 'product:',product
            print 'maxproducer:',maxProducer(model,product,fluxdict)
    print 'Number of reactions processed:',len(HFB)
    return

if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file

    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

    cobramodel.optimize(solver='gurobi')

    reactions = cobramodel.reactions
    metabolites = cobramodel.metabolites

    metabolite = metabolites.get_by_id('QH2_c')
    reaction = reactions.get_by_id('cyoABCD')

    HFB = HFBreactions(cobramodel,cobramodel.solution.x_dict)

    #iJO1366b = create_cobra_model_from_sbml_file('../SBML/iJO1366b.xml')

    #iJO1366b.optimize(solver='gurobi')

    #z = HFBreactions(iJO1366b)

    network = HFBnetwork(cobramodel,cobramodel.solution.x_dict)

    HFBtoSIF(network,'test.sif')

    metabolites = cobramodel.metabolites

    fluxdict = cobramodel.solution.x_dict
    m = metabolites.get_by_id("H_e")

    connections = connectHFBmetabolites(cobramodel,fluxdict,reaction.id)

    pgl = reactions.get_by_id('pgl')

    pglconnections = connectHFBmetabolites(cobramodel,fluxdict,pgl.id)

    #listreactions(cobramodel,metabolite,fluxdict)

    #maxConsumer(cobramodel,'QH2_c',fluxdict)

    #maxConsumer(cobramodel,'O2_c',fluxdict)

    #maxProducer(cobramodel,'H_e',fluxdict)

    #listreactions(cobramodel,'H_e',fluxdict)

    #maxProducer(cobramodel,'Q_c',fluxdict)

    #listreactions(cobramodel,'Q_c',fluxdict)

    HFBdetails(cobramodel,HFB,fluxdict, hold = True)

    
