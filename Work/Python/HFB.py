#HFB.py
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import timeit
import time
from time import strftime

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


def countK(metabolite_id,model,fluxdict,sense = 1):
    num = 0
    for reaction in model.metabolites.get_by_id(metabolite_id).get_reaction():
        massflux = reaction.get_coefficient(metabolite_id) * fluxdict[reaction.id]
        if cmp(massflux,0) == cmp(sense,0):
            num +=1
    return num


def getMassFlux(metabolite_id,model,reaction_id,fluxdict):
    return model.reactions.get_by_id(reaction_id).get_coefficient(metabolite_id)*fluxdit[reaction.id]

def fragmentation(metabolite_id,model,fluxdict, sense = 1):
    rx = {}
    for reaction in model.metabolites.get_by_id(metabolite_id).get_reaction():
        massflux = reaction.get_coefficient(metabolite_id) * fluxdict[reaction.id]
        if cmp(massflux,0) == cmp(sense,0):
            rx[reaction.id] = fluxdict[reaction.id]
    if sum(rx.values()) == 0:
        return 0
    else:
        Y = sum([(flux/sum(rx.values()))**2 for flux in rx.values()])
        return Y

def avgfragmentation(model,fluxdict,k,sense = 1):
    values = []
    for metabolite in model.metabolites:
        if countK(metabolite.id,model,fluxdict, sense = sense) == k:
            values.append(fragmentation(metabolite.id,model,fluxdict, sense = sense))
    return np.mean(values)




def plotfragmentation(cobramodel,fluxdict):
    kvals_in = list(set([countK(metabolite.id,cobramodel,fluxdict, sense = 1) for metabolite in cobramodel.metabolites]))
    kvals_out = list(set([countK(metabolite.id,cobramodel,fluxdict, sense = -1) for metabolite in cobramodel.metabolites]))

    yvals_in = [avgfragmentation(cobramodel,fluxdict,kval) for kval in kvals_in]
    yvals_out =[avgfragmentation(cobramodel,fluxdict,kval, sense = -1) for kval in kvals_out]
    fig = plt.figure()
    plt.scatter(kvals_in,yvals_in, c='k',marker='o')
    plt.scatter(kvals_out,yvals_out, c='r',marker='^')
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlim([1,10**3])
    plt.ylim([1,10**2])
    plt.show()



def randomizemedium(cobramodel,filename,percent, limit = -20,immune_metabolites = [],fixuptakes = False):
    #See fig 1, Almaas et al.
    #filename: text-file with list of substrates whose uptake should be allowed (set uptake flux bound to -1000)
    #percent: Percentage of substrates in list for which uptake should be allowed.

    #Generate a list of exchange reactions for which the flux limits should not be changed:
    immune_reactions = []
    for metabolite_id in immune_metabolites:
        metabolite = cobramodel.metabolites.get_by_id(metabolite_id)
        reactions = metabolite.get_reaction()
        immune_reactions +=[x.id for x in reactions if x.boundary == 'system_boundary']
                                

    #Set all exchange reactions to 0 before random selection of allowed fluxes:
    exchange_reactions = [x for x in cobramodel.reactions if len(x.get_reactants()) == 0 or len(x.get_products())==0]
    for rx in exchange_reactions:
        if rx.id not in immune_reactions:
            rx.lower_bound = 0
    
    metabolites = []
    #Read in a list of metabolites for which the exchange flux limits are eligible to be changed:
    with open(filename, 'r') as sourcefile:
        for line in sourcefile:
            if (not line.startswith('#')) and (len(line) > 1):
                metabolite_id = line.split(' ')[0].split('\n')[0].split('M_',1)[1]
                metabolites.append(metabolite_id)

    #Select only those metbolites which are not in the "do not change" list:
    eligible_metabolites = [metabolite_id for metabolite_id in metabolites if metabolite_id not in immune_metabolites]
    #Make a random sample of metabolites for which the exchange reaction limit should be changed:
    chosen_metabolites = random.sample(eligible_metabolites,int(math.floor(percent*0.01*len(eligible_metabolites))))

    #Get the corresponding exchange reactions for the chosen metabolites:
    chosen_reactions = []
    for metabolite_id in chosen_metabolites:
        metabolite = cobramodel.metabolites.get_by_id(metabolite_id)
        reactions = metabolite.get_reaction() 
        chosen_reactions +=[x.id for x in reactions if x.boundary == 'system_boundary']
    #print 'exchange reactions:',exchange_reactions
    #Change the reaction bounds for the chosen exchange reactions:
    for reaction_id in chosen_reactions:
        cobramodel.reactions.get_by_id(reaction_id).lower_bound = limit
        if fixuptakes:
            cobramodel.reactions.get_by_id(reaction_id).upper_bound = limit
        
    return exchange_reactions #Don't need the return value. The main purpose is to shuffle the reaction bounds in the model.
        



def fluxdistributionhistogram(cobramodel,filename,reaction_id,immune_metabolites = [],fixuptakes = False,trials = 20, percentage = 50, frequency = True, normalize = False, binN = 100):
    iterations = 0
    fluxvalues = []
    while iterations < trials:
        iterations +=1
        randomizemedium(cobramodel,filename,percentage,immune_metabolites = immune_metabolites, fixuptakes = fixuptakes)
        cobramodel.optimize()
        if normalize:
            fluxvalues.append(normalizefluxdict(cobramodel.solution.x_dict)[reaction_id])
        else:
            fluxvalues.append(cobramodel.solution.x_dict[reaction_id])
    n, bins = np.histogram(fluxvalues, bins = binN)
    if frequency:
        n = np.true_divide(n,sum(n))
    return n,bins

def plotsinglefluxdistribution(cobramodel,filename,reaction_id, trials = 20, percentage = 50, frequency = True, normalize = False, binN = 100):
    #See fig 4, Almaas et al.

    n, bins = fluxdistributionhistogram(cobramodel,filename,reaction_id,trials = trials, percentage = percentage, frequency = frequency, normalize = normalize, binN = binN)
        
    bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(bins_mean, n)
    print 'len(n):',len(n)
    print 'n:',n
    print 'sum(n):',sum(n)
    if frequency:
        #plt.ylim(0,1)
        plt.ylim(0,0.15)
    else:
        plt.ylim([0,trials])
    #plt.ylim(0,
    plt.show()
    
    pass


def plotfluxdistributions(cobramodel,filename,reaction_ids,positions,x_scales = [],immune_metabolites = [],trials = 20, percentage = 50, frequency = True, normalize = False, binN = 100, fixuptakes = False):
    fig = plt.figure()
    axes = []
    if len(x_scales) == 0:
        x_scales = [None for position in reaction_ids]
    for reaction_id,position,scaling in zip(reaction_ids,positions,x_scales):
        n, bins = fluxdistributionhistogram(cobramodel,filename,reaction_id,fixuptakes = fixuptakes,immune_metabolites = immune_metabolites,trials = trials, percentage = percentage, frequency = frequency, normalize = normalize, binN = binN)
        bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
        axes.append(fig.add_subplot(position))
        axes[-1].scatter(bins_mean, n)
        axes[-1].set_title(reaction_id)
        if frequency:
            #plt.ylim(0,1)
            plt.ylim(0,0.15)
        else:
            plt.ylim([0,trials])
        if scaling is not None:
            plt.xlim(list(scaling))
        #plt.ylim(0,
    plt.show()

def normalizefluxdict(fluxdict):
    '''
    Modify a flux dictionary to contain fluxes corresponding to the flux vector being normalized to unity.
    '''
    ids = []
    values = []
    for entry in fluxdict:
        ids.append(entry)
        values.append(fluxdict[entry])
    normalizedvalues = values / np.linalg.norm(values)
    return dict(zip(ids, normalizedvalues))

if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file
    from cobra import flux_analysis

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
    

    iJR904 = create_cobra_model_from_sbml_file('../SBML/Ec_iJR904.xml')
    iJR904.reactions.get_by_id('BiomassEcoli').objective_coefficient = 1
    bm = iJR904.reactions.get_by_id('BiomassEcoli')

    #Retrieve all exchange reaction
    exchange_reactions = [x for x in iJR904.reactions if len(x.get_reactants()) == 0 or len(x.get_products())==0]
    #alternative: system_boundary_reactions = [x for x in model.reactions if x.boundary == 'system_boundary']
    #See https://groups.google.com/forum/#!topic/cobra-pie/IPIOq30i-js

    iJE660a = create_cobra_model_from_sbml_file('iJE660a_fromMPS.sbml')

    iJE660a.reactions.get_by_id('EX_NH3').lower_bound = -100
    iJE660a.reactions.get_by_id('EX_O2').lower_bound = -20
    iJE660a.reactions.get_by_id('EX_O2').upper_bound = 0
    iJE660a.reactions.get_by_id('ATPM').lower_bound = 7.6
    iJE660a.reactions.get_by_id('ATPM').upper_bound = 7.6
    iJE660a.reactions.get_by_id('EX_CO2').lower_bound = -1000
    iJE660a.reactions.get_by_id('EX_K').lower_bound = -1000 #Potassium
    iJE660a.reactions.get_by_id('EX_PI').lower_bound = -1000
    iJE660a.reactions.get_by_id('EX_SLF').lower_bound = -1000

    #Set maximal carbon source uptake:
    #iJE660a.reactions.get_by_id('EX_GLU').lower_bound = -20
    iJE660a.reactions.get_by_id('EX_SUCC').lower_bound = -20
    

    for rx in exchange_reactions:
        rx.lower_bound = 0
    
    atpm = iJR904.reactions.get_by_id('ATPM')
    atpm.lower_bound = 7.6
    atpm.upper_bound = 7.6
    exo2 = iJR904.reactions.get_by_id('EX_o2_e_')
    exo2.lower_bound = -20
    exo2.upper_bound = 0
    nh4 = iJR904.reactions.get_by_id('EX_nh4_e_')
    nh4.lower_bound = -100

    iJR904.reactions.get_by_id('EX_co2_e_').lower_bound = -1000
    iJR904.reactions.get_by_id('EX_k_e_').lower_bound = -1000 #Potassium
    iJR904.reactions.get_by_id('EX_pi_e_').lower_bound = -1000
    iJR904.reactions.get_by_id('EX_so4_e_').lower_bound = -1000
    


    iJR904.optimize(solver = 'gurobi')

    print 'iJR904 solution:',iJR904.solution

    succCoA = iJR904.solution.x_dict['SUCOAS']

    #normflux = fluxvector / np.linalg.norm(fluxvector)

    #z = HFBreactions(iJO1366b)

    network = HFBnetwork(cobramodel,cobramodel.solution.x_dict)

    HFBtoSIF(network,'test.sif')

    metabolites = cobramodel.metabolites

    fluxdict = cobramodel.solution.x_dict
    #fluxdictGS = iJO1366b.solution.x_dict
    m = metabolites.get_by_id("H_e")

    connections = connectHFBmetabolites(cobramodel,fluxdict,reaction.id)

    pgl = reactions.get_by_id('pgl')

    pglconnections = connectHFBmetabolites(cobramodel,fluxdict,pgl.id)

   
    #HFBdetails(cobramodel,HFB,fluxdict, hold = True)


    '''
    #Plot fragmentation:
    kvals_in = list(set([countK(metabolite.id,cobramodel,fluxdict, sense = 1) for metabolite in cobramodel.metabolites]))
    kvals_out = list(set([countK(metabolite.id,cobramodel,fluxdict, sense = -1) for metabolite in cobramodel.metabolites]))

    yvals_in = [avgfragmentation(cobramodel,fluxdict,kval) for kval in kvals_in]
    yvals_out =[avgfragmentation(cobramodel,fluxdict,kval, sense = -1) for kval in kvals_out]

    kvals_inGS = list(set([countK(metabolite.id,iJO1366b,fluxdictGS, sense = 1) for metabolite in iJO1366b.metabolites]))
    kvals_outGS = list(set([countK(metabolite.id,iJO1366b,fluxdictGS, sense = -1) for metabolite in iJO1366b.metabolites]))

    yvals_inGS = [avgfragmentation(iJO1366b,fluxdictGS,kval) for kval in kvals_inGS]
    yvals_outGS = [avgfragmentation(iJO1366b,fluxdictGS,kval, sense = -1) for kval in kvals_outGS]

    fig = plt.figure()
    #plt.plot(kvals_inGS,yvals_inGS)
    plt.scatter(kvals_inGS,yvals_inGS, c='k',marker='o')
    plt.scatter(kvals_outGS,yvals_outGS, c='r',marker='^')
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlim([1,10**3])
    plt.ylim([1,10**2])
    plt.show()
    '''

    '''
    for metabolite in cobramodel.metabolites:
        if (countK(metabolite.id,cobramodel,fluxdict, sense = 1) +countK(metabolite.id,cobramodel,fluxdict, sense = 0) + countK(metabolite.id,cobramodel,fluxdict, sense = -1)) == len(metabolite.get_reaction()):
            print "it's good."
        else:
            print 'Uh-oh.'
    '''

    time.clock()
    #start = timeit.timeit()
    #plotsinglefluxdistribution(iJR904,'CO2t',trials = 5000, normalize = True)
    #plotsinglefluxdistribution(iJR904,'CO2t',trials = 1000)
    #plotsinglefluxdistribution(iJR904,'CO2t',trials = 1500)
    #plotsinglefluxdistribution(iJR904,'CO2t',trials = 2000)
    #plotsinglefluxdistribution(iJR904,'CO2t',trials = 5000)
    #end = timeit.timeit()
    #print end - start
    print time.clock()
    
    #plotsinglefluxdistribution(iJR904,'TPI',trials = 500, normalize = True)
    #plotsinglefluxdistribution(iJR904,'CO2t',trials = 500, normalize = True)
    #plotsinglefluxdistribution(iJR904,'NADK',trials = 500, normalize = True)
    #plotsinglefluxdistribution(iJR904,'GSNK',trials = 500, normalize = True)

    reaction_ids = ['TPI','NADK','CO2t','GSNK']
    positions = [221,222,223,224]
    #plotfluxdistributions(iJR904,'substrates.txt',reaction_ids,positions,trials = 500, percentage = 50, frequency = True, normalize = True, binN = 100)

    reaction_ids2 = ['TPIA','NADF','R0035','GSKx1']
    x_scales =[(-0.22,-0.08),(3e-6,8e-6),(-0.40,0.1),(0,0.008)]
    #random = randomizemedium(iJE660a,'iJE660a_substrates.txt',50, limit = -1000)

    #plotfluxdistributions(iJE660a,'iJE660a_substrates.txt',reaction_ids2,positions,x_scales = [],trials = 500, percentage = 50, frequency = True, normalize = True, binN = 100)    
    #plotfluxdistributions(iJE660a,'iJE660a_substrates.txt',reaction_ids2,positions,x_scales = x_scales,trials = 500, percentage = 50, frequency = True, normalize = True, binN = 100)

    immune = ['CO2_e','GLU_e','K_e','NH3_e','O2_e','PI_e','SLF_e']
    #plotfluxdistributions(iJE660a,'iJE660a_substrates.txt',reaction_ids2,positions,immune_metabolites = immune,x_scales = x_scales,trials = 500, percentage = 50, frequency = True, normalize = True, binN = 100)

    strftime("%Y-%m-%d %H:%M:%S")
    #plotfluxdistributions(iJE660a,'iJE660a_substrates.txt',reaction_ids2,positions,immune_metabolites = immune,x_scales = x_scales,trials = 2000, percentage = 50, frequency = True, normalize = True, binN = 100)
    plotfluxdistributions(iJE660a,'iJE660a_substrates.txt',reaction_ids2,positions,immune_metabolites = immune,x_scales = x_scales,trials = 500, percentage = 50, frequency = True, normalize = True, binN = 100)
   
