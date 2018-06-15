import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import os
from copy import copy
import numpy as np
import random as rand

import fswig_hklgen as H
import hkl_model as Mod
import sxtal_model as S

import  bumps.names  as bumps
import bumps.fitters as fitter
from bumps.formatnum import format_uncertainty_pm

#Simple Q learning algorithm to optimize a single parameter
#Will determine the optimal order of measurements to make
#to optimize the given parameter

np.seterr(divide="ignore",invalid="ignore")

#Set data files
DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = None
observedFile = os.path.join(DATAPATH,r"prnio.int")
infoFile = os.path.join(DATAPATH,r"prnio.cfl")

#Read data
spaceGroup, crystalCell, atomList = H.readInfo(infoFile)
# return wavelength, refList, sfs2, error, two-theta, and four-circle parameters
wavelength, refList, sfs2, error = S.readIntFile(observedFile, kind="int", cell=crystalCell)
tt = [H.twoTheta(H.calcS(crystalCell, ref.hkl), wavelength) for ref in refList]
backg = None
exclusions = []

def setInitParams():

    print("Setting parameters...")

    #Make a cell
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem)

    print(error)

    #Define a model
    m = S.Model(tt, sfs2, backg, wavelength, spaceGroup, cell,
                [atomList], exclusions,
                scale=0.06298, error=error, hkls=refList, extinction=[0.0001054])

    #Set a range on the x value of the first atom in the model
    m.atomListModel.atomModels[0].z.range(0,1)

    #Generate list of hkls
    hkls = []
    for r in m.reflections:
        hkls.append(r.hkl)

    return m, hkls


def fit(model):

    print("Fitting problem...")

    #Crate a problem from the model with bumps,
    #then fit and solve it
    problem = bumps.FitProblem(model)
    print(problem.labels())
    fitted = fitter.MPFit(problem)
    x, dx = fitted.solve()

    print(problem.nllf())
    problem.model_update()
    model.update()

    print(problem.nllf())
    return x, dx

#---------------------------------------
#Q learning methods
#---------------------------------------

def learn():

    #Q params
    epsilon = 1
    minEps = 0.01
    epsDecriment = 0.99

    qtable = []

    alpha = .01
    gamma = .9

    maxEpisodes = 1

    model, referenceHkls = setInitParams()
    maxSteps = len(referenceHkls)

    qtable = np.zeros([len(referenceHkls), len(referenceHkls)])    #qtable(state, action)

    for episode in range(maxEpisodes):

        model, remainingHkls = setInitParams()
        refs = model.reflections
        model.reflections = []
        state = 0
        prevDx = None

        for step in range(maxSteps):

	    reward = 0

            guess = rand.random()
            if (guess < epsilon):
                #Explore: choose a random action from the posibilities
                action = rand.choice(remainingHkls)

            else:
                #Exploit: choose best option, based on qtable
                qValue = 0
                for hkl in remainingHkls:
                    if (qtable[referenceHkls.index(state), referenceHkls.index(hkl)] > qValue):
                        qValue = qtable[referenceHkls.index(state), referenceHkls.index(hkl)]
                        action = hkl

            #No repeats
            remainingHkls.remove(action)

            #Find the data for this hkl value and add it to the model
            for reflection in refs:        #TODO, should this be adding hkls not refs?
                if (reflection.hkl == action):
                    model.reflections.append(reflection)
                    model.update()     #may not be necessary
                    break

            print("s, a", state, action)
            print (len(model.reflections))
	    if (step > 1):        #TODO not necessary
                x, dx = fit(model)
                print(model.reflections)
   	        reward -= 1
                if (prevDx != None and dx < prevDx):
                    reward += 1

                print("reward", reward)

                qtable[referenceHkls.index(state), referenceHkls.index(action)] =  qtable[referenceHkls.index(state), referenceHkls.index(action)] + \
                                                                                   alpha*(reward + gamma*(np.max(qtable[referenceHkls.index(state),:])) - \
                                                                                   qtable[referenceHkls.index(state), referenceHkls.index(action)])
                prevDx = dx

            state = action

            if (not remainingHkls):
                break

        #Decriment epsilon to explote more as the model learns
        epsilon = epsilon*epsDecriment
        if (epsilon < minEps):
           epsilon = minEps



learn()
