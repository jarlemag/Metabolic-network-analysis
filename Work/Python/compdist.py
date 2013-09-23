#compdist.py
'''
%Computes a distance measure (euclidean or taxicab) (weighed or unweighed) between experimental and computed fluxes/split ratios. 
%Syntax:
%ds = compdist(fluxvector,options,sense)
%Output: Euclidean or taxicab norm.
%Distance can be computed on basis of raw fluxes or split ratios. If using
%split ratios, set options.compsplits = 1;
%options:
%model_id
%exp_id
%normchoice ('euclidean' or 'taxicab')
'''
import scipy.io

def compdist(fluxvector,options,sense)

    #matlab code: load('expdata')
    expdata = scipy.io.loadmat('expdata.mat')
    if ('verbflag' in options) == False:
        options['verbflag'] = 0
    verbflag = options['verbflag']

    if ('model_id' in options) == False:
        options['model_id'] = 1

    if ('exp_id' in options) == False:
        options['exp_id'] = 1

    if ('usescaledfluxes' in options) == False:
        options['usescaledfluxes'] = 0

    if ('useweights' in options) == False:
        options['useweights'] = 0
    if ('normchoice' in options) == False:
        options['normchoice'] = 'euclidean']

    if ('excludereactions' in options) == True:
        excludereactions = options['excludereactions']
    else:
        excludereactions = 0

    if ('debugmode' in options) == False:
        options['debugmode'] = 0

    if ('compsplits' in options) == False:
        options['compsplits'] = 0

    if ('substractBM' in options) == False:
        options['substractBM'] = 0

    exp_id = options['exp_id']
    usescaledfluxes = options['usescaledfluxes']
    useweights = options['useweights']
    debugmode = options['debugmode']
    substractBM = options['substractBM']
    compsplits = options['compsplits']

    if verbflag == 2:
        print 'compdist.m: Excluding following experimental reactions:',excludereactions


        #%Extract/calculate vector with which the experimental will be compared:

            if compsplits == 0:
                Fcomp = extractflux(fluxvector,options) #Use raw fluxes
                if debugmode > 1:
                    print 'compdist.m: Using raw fluxes.'

            elif compsplits == 1:
                Fcomp = computesplits(fluxvector,options['model_id'] #Use split ratios
                if substractBM == 0:
                    if debugmode > 1:
                        print 'compdist.m:Using split ratios fluxes.'
                #matlab code: Fcomp = Fcomp(:,)
                elif substractBM == 1:
                    if debugmode > 1:
                        print 'compdist.m: Using split ratios (BM subtracted).'
                    #Matlab code:Fcomp = Fcomp(:,2)
    if exp_id == 1:
        Fluxvalues = expdata.perrenoud.abs.batch.aerobe.fluxvalues
    else:
        raise ValueError('Data not available.')

#%If split ratios are used, replace experimental raw fluxes with the
#%corresponding split ratios.
    if compsplits == 1:
        Fexp = computesplits(Fexp,0)


#%Update this file to accept both absolute (Perrenoud 2005) and
#%glucose-scaled (Schuetz 2007) data as options.

#Optional: Exclude some reactions:
    if excludereactions != 0:
        for i in range(len(excludereactions)):
            #matlab code:Fexp(excludreactions(i)) = 0

    if debugmode == 1:
        print 'Fexp:',Fexp
        print 'Fcomp:',Fcomp


    if usescaledfluxes == 1:
        scalingfactor = Fexp[0]/Fcomp[0] #Scale to glucose flux
        e = (Fexp/scalingfactor) - Fcomp
    else:
        e = Fexp - Fcomp

    if len(Fexp) != len(Fcomp):
        raise ValueError('compdist.m: Fexpand Fcomp vectors not same length.')
    else:
        n = len(Fexp)

    #matlab code:W = eye(n)
    if useweights == 1:
        #matlab code:sigma = Fluxvalues(:,2)
        #matlab code:inverse_sigma = 1./sigma
        #matlab code:sum_inverse_sigma = sum(inverse_sigma)
        for i in range(n):
            #matlab code:W(i,i) = (1/sigma(i))*(1/sum_inverse_sigma)
            #matlab code:diagonal_sum = sum(diag(W))
            #Computer the sum of the diagonal of W. Should equal 1:
            if diagonal_sum != 1:
                print 'diagonal_sum:',diagonal_sum
                print 'compdist.py: Warning: Sanity check failed. Sum of diagonal elements in W should equal 1.'
    if debugmode > 0:
        print 'W:',W
        print 'compdist.m: normchoice:',options['normchoice']
    if options['normchoice'] == 'euclidean':
        if sense == 1:
            #matlab code: ds = sqrt(e'*w*e)
        elif sense == -1:
            #matlab code: ds = -sqrt(e'*W*e) #Report the negative value of the distance. For use with fmincon.
    if debugmode >1:
        print 'compdist.m: Returned value:',ds
    return ds
                                      
