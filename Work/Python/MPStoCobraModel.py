##MPStoCobraModel.py
##Convert an MPS model file to SBML using CobraPy
##Syntax: MPStoCobraModel.py filename_in filename_out

import cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import write_cobra_model_to_sbml_file
import sys

def convertMPStoCobraModel(filename, verbose = False):
    keywords = ['NAME','ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA']
    with open(filename, 'r') as sourcefile:

        metabolite_list = []
        for line in sourcefile:
            if "*" in line:
                continue #Ignore commmented lines
            linestart = line.split(' ')[0].split('\n')[0]
            if linestart in keywords:
                if verbose:
                    print 'keyword found:',linestart
                lastkeyword = linestart
            if (lastkeyword == 'NAME'):
                modelname = line.split(' ')[-1]
                cobramodel = cobra.Model(modelname)
            #Add metabolites to model:
            if lastkeyword == 'ROWS' and (linestart not in keywords):
                    metabolite_id = line.strip().rsplit(' ',1)[1]
                    row_type = line.strip().rsplit(' ',1)[0].strip()
                    metabolite_list.append(metabolite_id)
                    new_metabolite = Metabolite(metabolite_id)
                    if new_metabolite.id.endswith('xt'):
                        new_metabolite.compartment = 'e'
                    else:
                        new_metabolite.compartment = 'c'
                    if row_type == 'E':
                        cobramodel.add_metabolites([new_metabolite])
                        if verbose:
                            print 'metabolite added:',new_metabolite
            #Construct reactions:
            if lastkeyword == 'COLUMNS' and (linestart not in keywords):
                rdata = line.strip().split()
                reaction_id = rdata[0]
                metabolite_id  = rdata[1]
                coef = float(rdata[2])
                if metabolite_id == 'COST':
                    continue #Ignore the "COST" metabolite
                if metabolite_id not in cobramodel.metabolites:
                    print 'metabolite_id:',metabolite_id
                    raise Exception ('Metabolite not found.')
                else:
                    current_metabolite = cobramodel.metabolites.get_by_id(metabolite_id)
                    
                if reaction_id not in cobramodel.reactions:
                    current_reaction = Reaction(reaction_id)
                    current_reaction.add_metabolites({current_metabolite:coef})
                    cobramodel.add_reactions(current_reaction)
                else:
                    current_reaction = cobramodel.reactions.get_by_id(reaction_id)
                    current_reaction.add_metabolites({current_metabolite:coef})
            if lastkeyword == 'BOUNDS' and (linestart not in keywords):
                data = line.strip().split()
                bound_type = data[0]
                reaction_id = data[2]
                if len(data) > 3:
                    bound_value = float(data[3])
                if bound_type == 'FR':
                    cobramodel.reactions.get_by_id(reaction_id).lower_bound = -1000
                    cobramodel.reactions.get_by_id(reaction_id).upper_bound = 1000
                elif bound_type == 'LO':
                    cobramodel.reactions.get_by_id(reaction_id).lower_bound = bound_value
                elif bound_type == 'UP':
                    cobramodel.reactions.get_by_id(reaction_id).upper_bound = bound_value
                elif bound_type ==  'FX':
                    try:
                        cobramodel.reactions.get_by_id(reaction_id).lower_bound = bound_value
                    except KeyError:
                        print cobramodel.reactions
                    cobramodel.reactions.get_by_id(reaction_id).upper_bound = bound_value
                #print bound_type,reaction_id,bound_values

        #Rename external metabolites
        for metabolite in cobramodel.metabolites:
            if metabolite.id.endswith('xt'):
                if verbose:
                    print 'renaming metabolite ',metabolite.id
                new_metabolite_id = metabolite.id.split('xt')[0] + '_e'
                metabolite.id = new_metabolite_id
                if verbose:
                    print 'new metabolite id:',metabolite.id
        #Combine in/out exchange reactions:
        for reaction in cobramodel.reactions:
            if reaction.id.endswith('xtO'): #If the reaction is an "out" exchange reaction
                metabolite_id = reaction.id.split('xtO')[0]
                in_reaction_id = metabolite_id + 'xtI' #Construct the ID of the corresponding "in" reaction
                if in_reaction_id in cobramodel.reactions: #If a corresponding "in" reaction exists:
                    bound_in = cobramodel.reactions.get_by_id(in_reaction_id).upper_bound #Find the maximal uptake rate defined by the "in" reaction upper bounds
                    reaction.lower_bound = -bound_in #Make the reaction a combined in/out exchange reaction by allowing negative flow (if allowed by the "in" reaction)
                    cobramodel.reactions.get_by_id(in_reaction_id).remove_from_model() #Remove redundant "in" reactions
                reaction.id = 'EX_' + metabolite_id #Rename the now combined exchange reaction
                    
                    
    #return rdata  
    return cobramodel


if __name__ == "__main__":
    
    filename = sys.argv[1]
    print 'filename:',filename
    target = sys.argv[2]
    print 'target',target

    model = convertMPStoCobraModel(filename)
    write_cobra_model_to_sbml_file(model,target)


