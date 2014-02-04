#MPStoCobraModel.py
import cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra import Model, Reaction, Metabolite

def convertMPStoCobraModel(filename):
    keywords = ['NAME','ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA']
    with open(filename, 'r') as sourcefile:


        '''
        for line in sourcefile:
            if line.startswith('NAME'):
                modelname = line.split(' ')[-1]

        '''

        '''
        rows = []
        for line in sourcefile:
            if line.split('\n')[0] == 'ROWS':
                while line.split('\n')[0] not in keywords:
                    rows.append(line)
        '''
        rows = []
        for line in sourcefile:
            linestart = line.split(' ')[0].split('\n')[0]
            if linestart in keywords:
                print 'keyword found:',linestart
                lastkeyword = linestart
            if (lastkeyword == 'NAME') and ("*" not in line):
                modelname = line.split(' ')[-1]
                cobramodel = cobra.Model(modelname)
            if lastkeyword == 'ROWS' and (linestart not in keywords):
                try: 
                    rows.append(line.strip().rsplit(' ',1)[1])
                except MemoryError:
                    print 'Memory error.'
                    return rows
            if lastkeyword == 'COLUMNS' and (linestart not in keywords):
                rdata = line.strip().split()
                reaction_id = rdata[0]
                metabolite_id  = rdata[1]
                coef = float(rdata[2])

                if metabolite_id not in cobramodel.metabolites:
                    new_metabolite = Metabolite(metabolite_id)
                else:
                    new_metabolite = cobramodel.metabolites.get_by_id(metabolite_id)
                    
                if reaction_id not in cobramodel.reactions:
                    new_reaction = Reaction(reaction_id)
                    new_reaction.add_metabolites({new_metabolite:coef})
                    cobramodel.add_reactions(new_reaction)
                else:
                    current_reaction = cobramodel.reactions.get_by_id(reaction_id)
                    current_reaction.add_metabolites({new_metabolite:coef})
                    
    #return rdata  
    return cobramodel


if __name__ == "__main__":
    filename = 'ecoli_iJE660a.mps'

    rows = convertMPStoCobraModel(filename)

    SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
