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
                if  massflux > maxproducer_flux:   #If the massflux is greater than the current recorded maximum production flux
                    maxproducer_flux = massflux #Record the current massflux as the new maximum production flux
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
    print 'rpe in maxproducerslist?','rpe' in maxproducerslist
    print 'rpe in maxconsumerslist?','rpe' in maxconsumerslist
    HFBrxset = set(maxproducerslist).intersection(maxconsumerslist)
   
   
    return HFBrxset


def listconsumers(metabolite):
    reactions = metabolite.get_reaction()
    

def listreactions(model,metabolite,fluxdict):
    if type(metabolite) == type(''):
        metabolite = model.metabolites.get_by_id(metabolite)
    reactions = metabolite.get_reaction()
    #print '\n'
    #print 'metabolite:',metabolite.id
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
    networkdict = defaultdict(list)
    #For every reaction in the HFB
    for reaction_id in HFBrxset:
        #print '-------------'
        #print 'HFB reaction:',reaction_id #DEBUG
        reaction = cobramodel.reactions.get_by_id(reaction_id)

        maxconsumed = []
        maxproduced = []
        #For every metabolite in the reaction
        for metabolite in (reaction.get_reactants() + reaction.get_products()):
            metabolite_id = metabolite.id
            #print 'reactant id:',reactant_id #DEBUG
            metabolite = cobramodel.metabolites.get_by_id(metabolite_id) 
            metabolite_reactions = metabolite.get_reaction() #Get the reactions in which the reactant partakes
            consumersdict = {}
            producersdict = {}
            #For every reaction in which the metabolite partakes:
            for rx in metabolite_reactions:
                coef = rx.get_coefficient(metabolite_id)
                massflux =  coef * fluxdict[rx.id]
                if massflux < 0: #If the metabolite is being consumed in that reaction:
                    consumersdict[rx.id] = massflux #Record the rate of consumption
                if massflux > 0: #If the metabolite is being produced in that reaction:
                    producersdict[rx.id] = massflux  #Record the rate of production
         
            maxconsumervalue = min(consumersdict.values())
            maxconsumers = [key for key in consumersdict if consumersdict[key] == maxconsumervalue]  #Find the highest-consuming reaction(s) for the metabolite
            maxproducervalue = max(producersdict.values()) #Find the highest-producing reaction for the metabolite
            maxproducers = [key for key in producersdict if producersdict[key] == maxproducervalue]
                       
            if reaction_id in maxconsumers: #If the the current reaction is a highest-consuming reaction for the metabolite
                maxconsumed.append(metabolite_id) #Add the metabolite to the list of metabolites being maximally consumed by this reaction
            if reaction_id in maxproducers: #If the highst-producing reaction is the current reaction
                maxproduced.append(metabolite_id) #Add the metabolite to the list of metabolites being maximally produced by this reaction
                
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
        print'Reaction:',reaction.id,'\n'
        reactants = reaction.get_reactants()
        products = reaction.get_products()
        for reactant in reactants:
            print 'reactant:',reactant
            print 'maxconsumer:',maxConsumer(model,reactant,fluxdict)
            listreactions(model,reactant,fluxdict)
        for product in products:
            print 'product:',product
            print 'maxproducer:',maxProducer(model,product,fluxdict)
            listreactions(model,product,fluxdict)
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

    print 'Run 1:\n'
    HFB = HFBreactions(cobramodel,cobramodel.solution.x_dict)
    print 'Run 2:\n'
    HFB2 = HFBreactions(cobramodel,cobramodel.solution.x_dict)

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

    #HFBdetails(cobramodel,HFB,fluxdict, hold = True)

    
