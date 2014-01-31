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
        
    HFBrxset = set(maxproducerslist).intersection(maxconsumerslist)
    return HFBrxset

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


def maxminflux(model,metabolite_id,fluxdict, minimize = False ):
    '''
    Find the most producing or most consuming reactions for a given metabolite and flux distribution
    '''
    if type(metabolite_id) == type(''):
        metabolite = model.metabolites.get_by_id(metabolite_id)
    else:
        metabolite = metabolite_id
    metabolite_reactions = metabolite.get_reaction()
    massfluxdict = {}
    for rx in metabolite_reactions:
        coef = rx.get_coefficient(metabolite_id)
        massfluxdict[rx.id] = coef * fluxdict[rx.id]
    if minimize:
        targetvalue = min(massfluxdict.values())
    else:
        targetvalue = max(massfluxdict.values())
    targetkeys = [key for key in massfluxdict if massfluxdict[key] == targetvalue]
    return targetkeys

def connectHFBmetabolites(model,fluxdict,reaction_id):
    if type(reaction_id) == type(''):
        reaction = model.reactions.get_by_id(reaction_id)
    else:
        reaction = reaction_id
    connections = defaultdict(list)
    #For every metabolite partaking in the reaction:
    for metabolite in (reaction.get_reactants()+reaction.get_products()):
        maxconsumed = []
        maxproduced = []
        max_consumers = maxminflux(model,metabolite.id,fluxdict, minimize = True)
        if metabolite.id in max_consumers:
            maxconsumed.append(metabolite.id)
        max_producers = maxminflux(model,metabolite.id,fluxdict)
        if metabolite.id in max_producers:
            maxproduced.append(metabolite.id)
    for reactant in maxconsumed:
        for product in maxproduced:
            connections[reactant].append(product)
    return connections

        
def HFBnetwork(cobramodel,fluxdict, grouped = False):
    HFBrxset = HFBreactions(cobramodel,fluxdict)
    networkdict = defaultdict(list)
    #For every reaction in the HFB
    for reaction_id in HFBrxset:
        reaction = cobramodel.reactions.get_by_id(reaction_id)
        maxconsumed = []
        maxproduced = []
        #For every metabolite in the reaction
        for metabolite in (reaction.get_reactants() + reaction.get_products()):
            metabolite_id = metabolite.id
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
            maxproducervalue = max(producersdict.values()) 
            maxproducers = [key for key in producersdict if producersdict[key] == maxproducervalue] #Find the highest-producing reaction(s) for the metabolite
                       
            if reaction_id in maxconsumers: #If the the current reaction is a highest-consuming reaction for the metabolite:
                maxconsumed.append(metabolite_id) #Add the metabolite to the list of metabolites being maximally consumed by this reaction
            if reaction_id in maxproducers: #If the current reaction is a highest-producing reaction for the metabolite:
                maxproduced.append(metabolite_id) #Add the metabolite to the list of metabolites being maximally produced by this reaction
        for metabolite in maxconsumed:
            for product in maxproduced:
                networkdict[metabolite].append(product)
    return networkdict

def HFBtoSIF(networkdict,filename):
    target = open(filename, 'w')
    for reactant in networkdict:
        for product in networkdict[reactant]:
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

   
    #HFBdetails(cobramodel,HFB,fluxdict, hold = True)

    
