#HFB.py

def HFB(cobramodel):
    maxproducers = {}
    maxconsumers = {}
    tol = 1e-6
    for metabolite in cobramodel.reactions:
        reactions = metabolite._reaction
        maxproducer_id = ''
        maxproducer_flux = 0
        maxconsumer_id = {}
        maxconsumer_flux = ''
        for reaction in reactions:
            coef = reaction.get_coefficient(metabolite.id)
            massflux = cobramodel.solution.x_dict[reaction.id] * coef
            if abs(massflux) > tol: #If flux is non-zero
                if  massflux > maxproducer_flux:   
                    maxproducer_flux = massflux
                    maxproducer_id = reaction.id
                elif massflux < maxconsumer_flux:
                    maxconsumer_flux = massflux
                    maxconsumer_id = reaction.id
            








if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file

    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

    cobramodel.optimize(solver='gurobi')
