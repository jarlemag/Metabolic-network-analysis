#cogs-fba.py
#Main module file for COGS-FBA
import numpy as np
import matplotlib.pyplot as plt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from matplotlib import cm

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

def drawFluxLimits(cobramodel,reactions):
    pass


def z_func(cobramodel,x,y,reactions,boundtolerance = 1e-3):
        cobramodel.reactions.get_by_id(reactions[0]).lower_bound = x - boundtolerance
        cobramodel.reactions.get_by_id(reactions[0]).upper_bound = x + boundtolerance

        cobramodel.reactions.get_by_id(reactions[1]).lower_bound = y - boundtolerance
        cobramodel.reactions.get_by_id(reactions[1]).upper_bound = y + boundtolerance
        cobramodel.optimize(solver='gurobi')
        return cobramodel.solution.f
        
def PhenotypePhasePlane(cobramodel,reactions,xlimits,ylimits, verbose = False, boundtolerance = 1e-3,noplot = False, wait = False):
    X = np.linspace(xlimits[0],xlimits[1])
    Y = np.linspace(ylimits[0],ylimits[1])
    X, Y = np.meshgrid(X, Y)
    #X = np.ndarray.flatten(X)
    #Y = np.ndarray.flatten(Y)
    Z = []

        
    Z = z_func(X, Y)
    points = zip(X,Y)
    if verbose:
        print '# of points:',len(points)
    if wait:
        message = raw_input('Press Enter to continue')
        if message == 'abort':
            print 'Exiting.'
            return
    its = 0
    failed = 0
    for point in points:
        its +=1
        cobramodel.reactions.get_by_id(reactions[0]).lower_bound = point[0] - boundtolerance
        cobramodel.reactions.get_by_id(reactions[0]).upper_bound = point[0] + boundtolerance

        cobramodel.reactions.get_by_id(reactions[1]).lower_bound = point[1] - boundtolerance
        cobramodel.reactions.get_by_id(reactions[1]).upper_bound = point[1] + boundtolerance
        
        cobramodel.optimize(solver='gurobi')
        objval = cobramodel.solution.f
        Z.append(objval)
        if objval is None:
            failed += 1
        if verbose:
            print 'Processing point #',its,'z =',objval

    #if verbose:
    print 'Failed optimizations:',failed

    print 'X:',len(X)
    print 'Y:',len(Y)
    print 'Z:',len(Z)
    if noplot:
        return [X,Y,Z]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #surf = ax.plot_surface(X, Y, Z)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


    
SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
print 'Optiizing:'
SCHUETZR.optimize(solver='gurobi')
print 'solution:',SCHUETZR.solution.f
cobramodel = SCHUETZR
xlimits = [-20,0]
ylimits = [0,20]
reactions = ['EX_GLC_e','o2']

#PhenotypePhasePlane(SCHUETZR,reactions,xlimits,ylimits, verbose = True)
#X,Y,Z = PhenotypePhasePlane(SCHUETZR,reactions,xlimits,ylimits, verbose = True, noplot = True)


X = [-1,-2,-3,-4]
Y = [-1,-2,-3,-4]
X, Y = np.meshgrid(X, Y)
#X = np.ndarray.flatten(X)
#Y = np.ndarray.flatten(Y)
Z = X*Y

fig = plt.figure()
ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
#surf = ax.plot_surface(X, Y, Z, cmap=cm.jet)
surf = ax.plot_surface(X, Y, Z, cmap=cm.RdYlGn)
#fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

