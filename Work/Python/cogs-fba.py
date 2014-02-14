#cogs-fba.py
#Main module file for COGS-FBA
import cobra
import numpy as np
import matplotlib.pyplot as plt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from matplotlib import cm
from matplotlib.patches import Polygon

from mpl_toolkits.mplot3d import Axes3D



#Class definitions:

'''
class Options:
    def __init__(self,model_id = 1,exp_id = 1,verbflag = False,debugmode = False,
                 constraints = 0,objective = 1,optreq = 1,usesplits = 0,FBAonly = 0,
                 iterations = 0, usegurobi = 1,maxwithgurobi = 0, usefmincon = 0,
                 preloadmodels = 0,excludereactions = 0,useweights = 0,makeplots = 0,
                 usescaledfluxes = 0,fluxreport = 1,computecorrcoef =1,normchoice ='euclidean',
                 usealloptions = 0,printfile = 0,usetokens = 1,plotFBA = 0,extractmode ='default',
                 fluxconstrainset = 0,reportsplits = 0,constrainglucose =2,constrainBM = 1,usealtmap = 0
                 )
                self.model_id = model_id
                self.exp_id = exp_id
                self.verbflag = vebflag
                self.debugmode = debugmode
                self.constraints = constraints
                self.objective = objective
                self.optreq = optreq
                self.usesplits = usesplits
                self.FBAonly = FBAonly
                self.iterations = iterations
                self.usegurobi = usegurobi
                self.maxwithgurobi = maxwithgurobi
                self.usefmincon = usefmincon
                self.preloadmodels = preloadmodels
                sef.excludereactions = excludereactions
                self.useweights = useweights
                self.makeplots = makeplots
                self.usescaledfluxes = usescaledfluxes
                self.fluxreport = fluxreport
                self.computecorrcoef = computecorrcoef
                self.normchoice = normchoice
                self.usealloptions = usealloptions
                self.printfile = printfile
                self.usetokens = usetokens
                self.plotFBA =plotFBA
                self.extractmode = extractmode
                self.fluxconstrainset = fluxconstrainset
                self.reportsplits = reportsplits
                self.constrainglucose = constrainglucose
                self.constrainBM = constrainBM
                self.usealtmaps = usealtmaps

'''


class SimpleCobraModel(object):
    def __init__(self,model_id):
        self.id = model_id



def createSimpleCobraModelfromXML(filename):
    pass


def writeModelToSIF(cobramodel,filename):
    target = open(filename,'w')
    for metabolite in cobramodel.metabolites:
        reactions = metabolite._reaction
        for reaction in reactions:
            coef = reaction.get_coefficient(metabolite.id)
            target.write(metabolite.id)
            target.write(' ')
            if coef < 0:
                target.write('reactant')
            else:
                target.write('product')
            target.write(' ')
            target.write(reaction.id)
            target.write('\n')
    target.close()
    return

def findInactiveReactions(cobramodel,optreq = 1, verbose = True):
    '''
    Find reactions in a model which are never active when Z >= optreq*Zmax
    '''
    FVAres = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel,fraction_of_optimum = optreq)
    inactive_set = [key for key in FVAres if (FVAres[key]['minimum'] == 0 and FVAres[key]['maximum'] == 0)]
    if verbose:
        print 'FVA result:'
        print FVAres
    return inactive_set

def plotInactiveReactions(cobramodel,optreqrange = [0,1]):
    optreqs = np.linspace(optreqrange, num = 20)
    inactive_sets = [findInactiveReactions(cobramodel,optreq) for optreq in optreqs ]
    inactive_number = [len(inactive_set) for inactive_set in inactive_sets]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Optimality requirement (Z/Z_max)')
    ax.set_ylabel('# Inactive reactions')
    ax.scatter(optreqs,inactive_number)
    plt.show()

def drawFluxSpace(cobramodel,reactionA,reactionB,):
    FVA_result = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel,fraction_of_optimum = 0,the_reactions=[reactionA])
    Amin = FVA_result[reactionA]['minimum']
    Amax = FVA_result[reactionA]['maximum']
    fluxlinspace = np.linspace(Amin,Amax)
    upper = []
    lower = []
    for flux in fluxlinspace:
        cobramodel.reactions.get_by_id(reactionA).lower_bound = flux
        cobramodel.reactions.get_by_id(reactionA).upper_bound = flux
        FVAres = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel,fraction_of_optimum = 0,the_reactions=[reactionB])
        minimum = FVAres[reactionB]['minimum']
        maximum = FVAres[reactionB]['maximum']

        upper.append([flux,maximum])
        lower.append([flux,minimum])
    vertices = upper + lower[::-1]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.add_patch(Polygon(vertices, closed = False, fill = True))
    ax.set_xlim(0,5)
    ax.set_ylim(0,5)
    plt.show()

    


def z_func(x,y,cobramodel,reactions,boundtolerance = 1e-3):
        #print 'Reactions:',reactions
        cobramodel.reactions.get_by_id(reactions[0]).lower_bound = x - boundtolerance
        cobramodel.reactions.get_by_id(reactions[0]).upper_bound = x + boundtolerance

        cobramodel.reactions.get_by_id(reactions[1]).lower_bound = y - boundtolerance
        cobramodel.reactions.get_by_id(reactions[1]).upper_bound = y + boundtolerance
        cobramodel.optimize(solver='gurobi')
        return cobramodel.solution.f


z_func = np.vectorize(z_func, excluded = ['cobramodel','reactions','boundtolerance'])
z_func.excluded.add(2)
z_func.excluded.add(3)

def PhenotypePhasePlane(cobramodel,reactions,xlimits,ylimits, verbose = False, boundtolerance = 1e-3,noplot = False, wait = False):
    X = np.linspace(xlimits[0],xlimits[1])
    Y = np.linspace(ylimits[0],ylimits[1])
    X, Y = np.meshgrid(X, Y)
    #X = np.ndarray.flatten(X)
    #Y = np.ndarray.flatten(Y)

    Z = z_func(X, Y,cobramodel,reactions)
    print 'len:',len(np.isnan(Z))
    Z[np.isnan(Z)] = 0


##    if verbose:
##        print 'Failed optimizations:',failed

    if verbose:
        print 'Maximum objective value:',Z.max()
    if noplot:
        return [X,Y,Z]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #surf = ax.plot_surface(X, Y, Z)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,
        linewidth=0, antialiased=False)
    ax.set_xlabel(reactions[0])
    ax.set_ylabel(reactions[1])
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlim(xlimits[0],xlimits[1])
    ax.set_ylim(ylimits[0],ylimits[1])
    plt.show()


if __name__ == "__main__":

    SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    xlimits = [-20,0]
    ylimits = [0,20]
    reactions = ['EX_GLC_e','o2']

    #PhenotypePhasePlane(SCHUETZR,reactions,xlimits,ylimits, verbose = True)
    #X,Y,Z = PhenotypePhasePlane(SCHUETZR,reactions,xlimits,ylimits, verbose = True, noplot = True)

    #gx = create_cobra_model_from_sbml_file('../SBML/gxfba_example.xml')
    #A = drawFluxSpace(gx,'R2_lowergly','R5_akgdh')
