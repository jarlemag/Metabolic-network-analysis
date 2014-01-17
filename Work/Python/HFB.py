#HFB.py
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
    pass
            

if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file

    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

    cobramodel.optimize(solver='gurobi')

    HFB = HFBreactions(cobramodel)

    iJO1366b = create_cobra_model_from_sbml_file('../SBML/iJO1366b.xml')

    iJO1366b.optimize(solver='gurobi')

    z = HFBreactions(iJO1366b)
