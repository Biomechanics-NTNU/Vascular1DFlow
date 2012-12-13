########################################################################################
#                            Vascular1dFlow v0.2
########################################################################################

import time 
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
#sys.path.append(cur+'/Solver')
#from Flow1DSolver import *

sys.path.append(cur+'/NetworkLib')
from classVascularNetwork import VascularNetwork

sys.path.append(cur+'/Solver')
from class1DflowSolver import FlowSolver

sys.path.append(cur+'/UtilityLib')
from modulXML import loadNetworkFromXML

sys.path.append(cur+'/Visualisation')
from class3dVisualisation import Visualisation

import pprint
import matplotlib.pyplot as plt

from optparse import OptionParser
import cPickle

import subprocess

def main():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest='networkName',
                      help="open file with networkName", metavar="FILE")
    parser.add_option("-s", "--save", dest="save", 
                      help="save solution data, 0 = False, 1 = True")
    
    parser.add_option("-n", "--dataNumber", dest='dataNumber',
                      help="number of the solution data (last number in filename), default = 1, max 999", metavar = "FILE")
    
    parser.add_option('-v', '--vizBool', dest='vizBool',
                      help="choose visualisation mode, 0: no visualisation, 1: 2d and 3d, 2: 2d plots, 3: 3d visualisation default = 1")
    
    (options, args) = parser.parse_args()
    
    networkName = 'oldNetwork1'
    if options.networkName != None:
        networkName = options.networkName
    else: 
        print " no networkName passed, visualisation could not be intialized"
        print " use -f networkName to define an file you want to visualize"
        exit()
    filename = str(networkName+'.xml')
    
    dataNumber = '001'
    if options.dataNumber != None:
        dataNumber = options.dataNumber
        if len(dataNumber) == 1:
            dataNumber = '00'+dataNumber
        if len(dataNumber) == 2:
            dataNumber = '0'+dataNumber  
    
    save = False
    if options.save != None:
        if options.save == '0':
            save = False
        elif options.save == '1':
            save = True
        
    vizBool = '1'
    vizBool2d  = True
    vizBool3d = True
    print options.vizBool
    if options.vizBool != None:
        if options.vizBool == '0':
            vizBool2d = False
            vizBool3d = False
        elif options.vizBool == '1':
            vizBool2d = True
            vizBool3d = True
        elif options.vizBool == '2':
            vizBool2d = True
            vizBool3d = False
        elif options.vizBool == '3':
            vizBool2d = False
            vizBool3d = True
      
    # load network from the path!
    Lstart = time.clock()
    vascularNetwork = loadNetworkFromXML(filename=filename) 
    
    if vascularNetwork == None: exit()
    # update the connections of the network, as they prob weren't saved correct
    vascularNetwork.evaluateConnections()
    Lend = time.clock()
    print" needed %1.3f  sec to load" %(Lend-Lstart)
    
    vascularNetwork.initialize()
    
    #create visualisation
    if vizBool3d == True:
        visualisation = Visualisation(vascularNetwork=vascularNetwork)
    
    
    #print '====================================='
    #print '_________ Simulation Context _________'
    #print '%-25s %-20s' % ('Numerical scheme',networkSolver['NumScheme'])
    #print '%-25s %4d' % ('Equation system',int(networkSolver['EqSystem']))
    #print '%-25s %4d' % ('Characteristic system',int(networkSolver['CharSystem']))
    #print '________ Visualisation Context _______'
    
    
    Cstart = time.clock()
    flowSolver = FlowSolver(vascularNetwork)
    P,Q,A = flowSolver.solve()
    print '\n____________ Solver time _____________'
    minutes = int((time.clock()-Cstart)/60.)
    secs = (time.clock()-Cstart)-float(minutes)*60.
    print '%d %s %f %s' % (minutes ,'min',secs ,'sec')   
    
    if save == True:
        #save solution data in c pickle
        #standart solution path
        networkPath = str(cur+"/NetworkFiles/")
        networkDirectory = filename.split('.')[0]
        solutionPath = str(networkPath+networkDirectory+'/')
        # data to save
        solutionDataSave = [{'Pressure': P, 'Flow': Q, 'Area': A, 'WaveSpeed':Q,'Name': str('simulation_'+dataNumber) }]
        #create file with filename in soultionPath
        endName = "_SolutionData_"+dataNumber+".pickle"
        FILE = open(str(solutionPath+vascularNetwork.name+endName),"w")
        # store pickle
        cPickle.dump(solutionDataSave, FILE, protocol=2)
        FILE.close()
        
    if vizBool2d == True:
        string = ' '.join(['python',cur+'/Visualisation/class2dVisualisationSimple.py','-f',vascularNetwork.name, '-n 1','-i 0'])                
        subprocess.call(string, shell=True)
    
    if vizBool3d == True:
        # create solution data
        solutionData = [ {'Pressure': P, 'Flow': Q, 'Area': A, 'Name': str('simulation '+dataNumber) }]
        #set the solution data
        visualisation.setSolutionData(solutionData)
        visualisation.factor = 10
        # set the coefficient for the velocity profile
        visualisation.powerLawCoefficient = 2       
        visualisation.visualize()
        # clear memory
        visualisation.__del__()
        del(visualisation)
        

if __name__ == '__main__':
    main()
