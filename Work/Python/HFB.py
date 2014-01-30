#HFB.py
from collections import defaultdict

def HFBreactions(cobramodel):
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
            massflux = cobramodel.solution.x_dict[reaction.id] * coef
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

    HFBrxset = set(maxproducerslist).intersection(maxconsumerslist)
   
    return HFBrxset

def HFBnetwork(cobramodel):
    HFBrxset = HFBreactions(cobramodel)
    #networkdict = {}
    networkdict = defaultdict(list)
    #For every reaction in the HFB
    for reaction_id in HFBrxset:
        print '-------------'
        print 'HFB reaction:',reaction_id #DEBUG
        reaction = cobramodel.reactions.get_by_id(reaction_id)

        maxconsumed = []
        maxproduced = []
        #For every reactant in the reaction
        for reactant in reaction.get_reactants():
            reactant_id = reactant.id
            print 'reactant id:',reactant_id #DEBUG
            reactant = cobramodel.metabolites.get_by_id(reactant_id) 
            metabolite_reactions = reactant.get_reaction() #Get the reactions in which the reactant partakes
            consumersdict = {}
            #For every reaction in which the reactant partakes:
            for rx in metabolite_reactions:
                coef = rx.get_coefficient(reactant_id)
                if coef < 0: #If the reactant is being consumed in that reaction:
                    consumersdict[rx.id] = coef * cobramodel.solution.x_dict[rx.id] #Record the rate of consumption
            maxconsumer =  max(consumersdict, key = consumersdict.get) #Find the highest-consuming reaction for the reactant
            print 'max consumer:',maxconsumer #DEBUG
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
                if coef > 0: #If the product is being produced in that reaction:
                    producersdict[rx.id] = coef * cobramodel.solution.x_dict[rx.id] #Record the rate of production
        maxproducer = max(producersdict, key = producersdict.get) #Find the highest-producing reaction for the product
        print 'max producer:',maxproducer #DEBUG
        if maxproducer == reaction_id: #If the highst-producing reaction is the current reaction
            maxproduced.append(product_id)
            
        for metabolite in maxconsumed:
            for product in maxproduced:
                networkdict[metabolite].append(product)
    return networkdict

if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file

    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

    cobramodel.optimize(solver='gurobi')

    HFB = HFBreactions(cobramodel)

    #iJO1366b = create_cobra_model_from_sbml_file('../SBML/iJO1366b.xml')

    #iJO1366b.optimize(solver='gurobi')

    #z = HFBreactions(iJO1366b)

    network = HFBnetwork(cobramodel)
