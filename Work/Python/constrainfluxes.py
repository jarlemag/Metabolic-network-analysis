#constrainfluxes.py

def constrainfluxes(model,expfluxdict,experrdict,reactionmap,tolerance = 1, fluxconstrainset = None,createTokens = False, verbose = True, debug = False):

    defaulterror = 0
    for linkdict in reactionmap:
        exprxns = linkdict['exprxns']
        modrxns = linkdict['modrxns']
        if len(exprxns) == 1:
            expid = exprxns[0]['rxid']
        else:
            pass #Do something sensible in case the reactionmap entry contains several experimental reactions.
        if fluxconstrainset is not None:
            if expid in fluxconstrainset:
                pass
            else:
                continue

        if len(modrxns) > 1: #If more than one model reaction is listed for a experimental reaction
            token_id = 't{exprxn}'.format(exprxn = expid)
            if token_id in model.reactions:
                print 'Token reaction found in model for experimental reaction ',expid
                modrxnid = token_id
            else:
                print 'Link contains more than model reaction, but no token reaction was found. Skipping.'
        else:
            modrxnid = modrxns[0]['rxid']
            coef = modrxns[0]['coef']

        if expid in expfluxdict:
            expflux = expfluxdict[expid]
            experror = experrdict[expid]
        else:
            print 'No flux value found for experimental reaction ',expid,'. Skipping reaction.'
            continue

        modelreaction = model.reactions.get_by_id(modrxnid)
        if verbose:
            print 'Constraining experimental reaction ',expid

        if debug:
            print 'coef:',coef

        if coef > 0:
            modelreaction.upper_bound = expflux+(experror*tolerance)
            modelreaction.lower_bound = expflux-(experror*tolerance)
        elif coef >0:
            modelreaction.upper_bound = -(expflux+experror*tolerance)
            modelreaction.lower_bound = -(expflux-experror*tolerance)
        else:
            print 'Unexpected reaction coefficient in reaction map. Skipping reaction.'
            continue 

        if verbose:
            print 'New upper bound:',modelreaction.upper_bound
            print 'New lower bound:',modelreaction.lower_bound
    return model

if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file
    import loadData as load

    model = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')
    expfluxdict = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')
    experrdict = load.ExpErrorsfromXML('expdata.xml','Perrenoud','Batch','aerobe')
    
    constrainedmodel = constrainfluxes(model,expfluxdict,experrdict,rmap, debug = True) 
    
