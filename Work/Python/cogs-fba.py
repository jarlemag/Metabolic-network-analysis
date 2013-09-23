#cogs-fba.py
#Main module file for COGS-FBA





#Class definitions:


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
