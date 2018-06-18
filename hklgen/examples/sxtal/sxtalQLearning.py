import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
import  bumps.names  as bumps
import bumps.fitters as fitter
from bumps.formatnum import format_uncertainty_pm
import random as rand

#Simple Q learning algorithm to optimize a single parameter
#Will determine the optimal order of measurements to make
#to optimize the given parameter

np.seterr(divide="ignore",invalid="ignore")

#Set data files
DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,r"hobk_bas.bac")
observedFile = os.path.join(DATAPATH,r"hobk.dat")
infoFile = os.path.join(DATAPATH,r"hobk1.cfl")

#Read collected data
(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]
wavelength = 2.524000
ttMin = 10.010000228881836
ttMax = 89.81000518798828
ttStep = 0.20000000298
exclusions = []
tt, observed, error = H.readIllData(observedFile, "D1B", backgFile)
observedByAlg = [] #currently, the algorithm hasn't measured any datapoints
backg = H.LinSpline(None)
basisSymmetry = copy(symmetry)


def setInitParams():

    print("Setting parameters...")

    #Make a cell
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem)

    #Define a model
    m = Mod.Model(tt, observedByAlg, backg, 1.548048,-0.988016,0.338780, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=6512, scale=59.08, eta=0.0382, zero=0.08416, error=error)


    #Set a range on the x value of the first atom in the model
    m.atomListModel.atomModels[0].x.range(0,1)

    #Generate list of hkls
    hkls = []
    for r in m.reflections:
        hkls.append(r.hkl)


#    m.observed.append(m.reflections[0])
#    m.observed.append(m.reflections[1])
#    m.observed.append(m.reflections[2])
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

#def main():
#    uvw = [1.548048,-0.988016,0.338780]
#    cell = crystalCell
#    H.diffPattern(infoFile=infoFile, uvw=uvw, cell=cell, scale=59.08,
#                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, wavelength = wavelength,
#                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
#                  magnetic=True, info=True, plot=False,
#                  observedData=(tt,observed), base=6512, residuals=True, error=error)
#    print("calling fit")
#    problem = fit()
#    setInitParams()



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
        state = 0

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
            for reflection in model.reflections:
                if (reflection.hkl == action):
                    observedByAlg.append(reflection)
                    #model.observed.append()
                    model.update()     #may not be necessary

            print("s, a", state, action)

	    if (step > 1):
                x, dx = fit(model)
                print("dx", dx)
    	        reward -= 1
                if (dx < prevDx):
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
