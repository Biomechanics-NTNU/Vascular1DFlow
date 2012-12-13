
try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree

from numpy import sqrt,pi,exp,sin,cos,abs,tan
import numpy as np

import os,sys
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../'+'/NetworkLib')

from classVascularNetwork import VascularNetwork
from classBoundaryConditions import *

from constants import unitsDict
from constants import bcTags
from constants import bcTagsUnits
from constants import bcTagsClassReferences
from constants import bcDuMatrix

def writeNetworkToXML(vascularNetwork,filename = None,networkPath = str(cur+"/../NetworkFiles/")):
    
    #print filename
    if filename == None:
        filename = "oldNetwork1.xml"
    
    networkDirectory = filename.split('.')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))  
    
    try:
        root = etree.Element(filename, id= "1.0", version="2.0")
    except:
        print " Error: path / file does not exist"
        return
    
    
    xmlFile = etree.ElementTree(root)
    
    #SimulationContext
    simContext = etree.SubElement(root, "simulationContext")
    CFL = etree.SubElement(simContext, "CFL")
    CFL.text = str(vascularNetwork.simulationContext['CFL'])
    totalTime = etree.SubElement(simContext, 'totalTime', unit='s')
    totalTime.text = str(vascularNetwork.simulationContext['totalTime'])
    NumScheme = etree.SubElement(simContext, "NumScheme")
    NumScheme.text = str(vascularNetwork.networkSolver['NumScheme'])
    EqSystem = etree.SubElement(simContext, "EqSystem")
    EqSystem.text = str(vascularNetwork.networkSolver['EqSystem'])
    CharSystem = etree.SubElement(simContext, "CharSystem")
    CharSystem.text = str(vascularNetwork.networkSolver['CharSystem'])
 
    
    #BoundaryConditions
    BCs = etree.SubElement(root, "boundaryConditions")
    for vesselId,boundaryInstances in vascularNetwork.boundaryConditions.iteritems():
        
        bc = etree.SubElement(BCs, "boundaryCondition", vessel_id = str(vesselId))
        for boundaryInstance in boundaryInstances:
            
            bcDict = boundaryInstance.__dict__
            bcName = bcDict['name']
            
            bcAttributeParent = bc
            if bcName != 'Rt': 
                bcT= etree.SubElement(bc, str(bcName))
                bcAttributeParent = bcT
            
            for attribute,value in bcDict.iteritems():
                if attribute in bcTags[bcName]:
                    element = etree.SubElement(bcAttributeParent, attribute, unit= bcTagsUnits[attribute])
                    
                    element.text = str(value)

    #global Fluid
    globalFluid = etree.SubElement(root, "globalFluid")
    my = etree.SubElement(globalFluid, "my" , unit='Pa s')
    my.text = str(vascularNetwork.globalFluid['my'])
    rho = etree.SubElement(globalFluid, "rho", unit='kg m-3')
    rho.text = str(vascularNetwork.globalFluid['rho'])
    gamma = etree.SubElement(globalFluid, 'gamma' , unit='')
    gamma.text = str(vascularNetwork.globalFluid['gamma'])
    dlt = etree.SubElement(globalFluid, 'dlt' , unit='')
    dlt.text = str(vascularNetwork.globalFluid['dlt'])
    pref = etree.SubElement(globalFluid, 'pref', unit='Pa')
    pref.text = str(vascularNetwork.globalFluid['pref'])
    
    #vessels
    vessels = etree.SubElement(root, "vessels")
    
    for vessel in vascularNetwork.vessels.itervalues():
        
        vessel_j = etree.SubElement(vessels, "vessel", id = str(vessel.Id),name = str(vessel.name), start_node = str(vessel.start),
                                    end_node = str(vessel.end),leftDaughter =  str(vessel.leftDaughter),
                                    rightDaughter = str(vessel.rightDaughter), angleToMother = str(vessel.angleToMother))
        # creating Grid references
        grid = etree.SubElement(vessel_j,"grid")
        
        geom = etree.SubElement(grid, "geom")
        geom.text = str(vessel.geom) 
        
        N = etree.SubElement(grid, "N")
        scalar = etree.SubElement(N, "scalar" , unit='')
        scalar.text = str(vessel.N)
        interval = etree.SubElement(N, "interval", unit='')
        if vessel.NInterval: interval.text = ' '.join(str(i) for i in vessel.NInterval)
        
        length = etree.SubElement(grid, "length")
        scalar = etree.SubElement(length, "scalar" , unit='m')
        scalar.text = str(vessel.length)
        interval = etree.SubElement(length, "interval", unit='m')
        if vessel.lengthInterval: interval.text = ' '.join(str(i) for i in vessel.lengthInterval)
        
        radiusA = etree.SubElement(grid, "radiusA")
        scalar = etree.SubElement(radiusA, "scalar", unit='m')
        scalar.text = str(vessel.radiusA)
        interval = etree.SubElement(radiusA, "interval", unit='m')
        if vessel.radiusAInterval: interval.text = ' '.join(str(i) for i in vessel.radiusAInterval)
               
        radiusB = etree.SubElement(grid, "radiusB")
        scalar = etree.SubElement(radiusB, "scalar", unit='m')
        scalar.text = str(vessel.radiusB)
        interval = etree.SubElement(radiusB, "interval", unit='m')
        if vessel.radiusBInterval: interval.text = ' '.join(str(i) for i in vessel.radiusBInterval)
        
        # creating solid referneces
        solid = etree.SubElement(vessel_j, "solid")
        
        comp = etree.SubElement(solid, "comp")
        comp.text = str(vessel.comp)
        
        Pfunc = etree.SubElement(solid, "Pfunc")
        Pfunc.text = str(vessel.Pfunc)
        Ps = etree.SubElement(solid, "Ps", unit='Pa')
        Ps.text = str(vessel.Ps)
        As = etree.SubElement(solid, "As", unit='m2')
        As.text = str(vessel.As)
        
        beta = etree.SubElement(solid, "beta")
        scalar = etree.SubElement(beta, "scalar", unit='Pa m-1')
        scalar.text = str(vessel.beta)
        interval = etree.SubElement(beta, "interval", unit='Pa m-1')
        if vessel.betaInterval: interval.text = ' '.join(str(i) for i in vessel.betaInterval)
     
        wallThickness = etree.SubElement(solid, "wallThickness")
        scalar = etree.SubElement(wallThickness, "scalar", unit='m')
        scalar.text = str(vessel.wallThickness)
        interval = etree.SubElement(wallThickness, "interval", unit='m')
        if vessel.wallThicknessInterval: interval.text = ' '.join(str(i) for i in vessel.wallThicknessInterval)
    
        youngModulus = etree.SubElement(solid, "youngModulus")
        scalar = etree.SubElement(youngModulus, "scalar", unit='Pa')
        scalar.text = str(vessel.youngModulus)
        interval = etree.SubElement(youngModulus, "interval", unit='Pa')
        if vessel.youngModulusInterval: interval.text = ' '.join(str(i) for i in vessel.youngModulusInterval) 
        
        # creating fluid refernces
        fluid = etree.SubElement(vessel_j, "fluid")
        if vessel.my is not None:
            my = etree.SubElement(fluid, "my" , unit='Pa s')
            my.text = str(vessel.my)
        if vessel.rho is not None:
            rho = etree.SubElement(fluid, "rho", unit='kg m-3')
            rho.text = str(vessel.rho)
        if vessel.gamma is not None:
            gamma = etree.SubElement(fluid, 'gamma', unit='')
            gamma.text = str(vessel.gamma)
        if vessel.dlt is not None:
            dlt = etree.SubElement(fluid, 'dlt', unit='')
            dlt.text = str(vessel.dlt)
        if vessel.pref is not None:
            pref = etree.SubElement(fluid, 'pref', unit='')
            pref.text = str(vessel.pref)
            

    
    xmlFile.write(networkPath+networkDirectory+'/'+filename,encoding='iso-8859-1',pretty_print = True)
    print " ... file saved"

def unitConversionToSI(unitsDict,value,unit):
        
    if ' ' in unit:
        unit = unit.split(' ')
        for un in range (0,len(unit),1): value = value*unitsDict[unit[un]]
    else: value = float(value)*unitsDict[unit]
    return value
   
def loadNetworkFromXML(filename = None, update = None, networkPath = str(cur+"/../NetworkFiles/")):
    
    # read from file
    if filename == None:
        filename = "oldNetwork1.xml"
        
    networkDirectory = filename.split('.')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        print 'Error: file not in the correct directory'
        return
    
    
    # create vascularNetwork instance
    if update == None: 
        vascularNetwork = VascularNetwork()
        # set name
        vascularNetwork.name = (filename.split('.'))[0]
    else: 
        netBC = {}
        netFluid = {}
        netSimCon = {}
        netSolver = {}
        networkDataDict = {'boundaryConditions': netBC, 'globalFluid':netFluid, 'simulationContext': netSimCon, 'networkSolver': netSolver }
        vesselDataDict = {}
     
    # load the data!
    try:
        parser = etree.XMLParser(encoding='iso-8859-1')
        tree = etree.parse(networkPath+networkDirectory+'/'+filename, parser)
    except:
        print " Error: path / file does not exist or error in file"
        if update == None: 
            return 
        else:
            return {'networkData': networkDataDict, 'vesselData':vesselDataDict}
   
    # create root
    root = tree.getroot()
    
    #SimulationContext
    for simContext in root.findall(".//simulationContext"):
        for data in simContext:
            if data.tag == 'CFL':
                if data.text != 'None' and data.text != None: 
                    if update == None: vascularNetwork.simulationContext['CFL'] = float(eval(data.text))
                    else: netSimCon['CFL'] = float(eval(data.text))
            if data.tag == 'totalTime':
                if data.text != 'None' and data.text != None: 
                    if update == None: vascularNetwork.simulationContext['totalTime'] = unitConversionToSI(unitsDict, float(data.text), data.attrib['unit'])
                    else: netSimCon['totalTime'] = unitConversionToSI(unitsDict, float(data.text), data.attrib['unit'])

            if data.tag == 'EqSystem' or data.tag == 'NumScheme':
                if data.text != 'None' and data.text != None: 
                    if update == None: vascularNetwork.networkSolver[data.tag] = str(data.text)
                    else: netSolver[data.tag] = str(data.text)
            if data.tag == 'CharSystem':
                if data.text != 'None' and data.text != None: 
                    try: charsys = int(data.text)
                    except: raise TypeError('Value of %s in simulationContext not convertable to int' %data.tag )   
                    if update == None: vascularNetwork.networkSolver[data.tag] = charsys
                    else: netSolver[data.tag] = charsys
    
    #load boundary Conditions
    for bcs in root.findall(".//boundaryConditions/boundaryCondition"):
        bcDict = bcs.attrib         
        if 'vessel_id' in bcDict:
            try: vessel_id = int(bcDict['vessel_id'])
            except: raise TypeError('vessel_id %s in boundaryCondition not convertable to int' %bcDict['vessel_id'] )
            typeData2 = {}
            
            boundaryInstances = []
            boundaryIntervals = []
            
            for data in bcs:
                                
                if data.tag in bcTags.keys():
                    boundaryInstance = eval(bcTagsClassReferences[data.tag])()
                    boundaryDataDict = {}
                    boundaryDataDict['name']= data.tag
                    boundaryIntervalDataDict = {}
                    
                    bcItems = bcTags[data.tag]
                                        
                    for bcItem in bcItems:
                        # find scalar Data
                        for value in data.findall(str(".//"+bcItem)):
                            try: floatValue = float(eval(value.text))
                            except: raise TypeError('Value %s in boundaryCondition not convertable to float' %value.text )
                            boundaryDataDict[bcItem] = unitConversionToSI(unitsDict, floatValue, value.attrib['unit'])
                          # find interval Data for polychaos
                        for value in data.findall(".//"+bcItem+"-interval"):
                            if value.text != 'None' and value.text !=  None: 
                                interval = []
                                try:
                                    for number in map(float, value.text.split(' ')):
                                        interval.append(unitConversionToSI(unitsDict, float(number), value.attrib['unit']))
                                except: raise TypeError('Value %s in boundaryCondition Interval not convertable to float' %value.text )  
                                  
                                boundaryIntervalDataDict[bcItem] =  interval
                                boundaryIntervalDataDict['name']= data.tag
                      
                    # exception for reflection coefficient because the value can stand rigth begind the       
                    if data.tag in ['Rt','_Rt']:
                        if type(data.text) != type(None):
                            if '.' in data.text:
                                try: floatValue = float(eval(data.text))
                                except: raise TypeError('Value %s in boundaryCondition not convertable to float' %data.text )
                                
                                # compatibility old xml-version Rt no unit-tag new version <Rt unit=''>..  
                                Rtunit = ''
                                if 'unit' in data.attrib: Rtunit = data.attrib['unit']
                                
                                boundaryDataDict[bcItem] = unitConversionToSI(unitsDict, floatValue, Rtunit)
                                                
                    #check if enough parameters are given:
                    if len(boundaryDataDict) != len(bcItems)+1: 
                        print len(boundaryDataDict) 
                        raise TypeError(str('boundaryCondition not sufficient defined; check if these attributes are given '+str(bcItems)))
                     
                    # set du-matrix correct for type 1
                    if data.tag in bcDuMatrix.keys(): boundaryDataDict['duMatrix'] = bcDuMatrix[data.tag]        
                    # update boundarz condition                           
                    boundaryInstance.update(boundaryDataDict)
                    boundaryInstances.append(boundaryInstance)
                
                    boundaryIntervals.append(boundaryIntervalDataDict)

            if update == None:                 
                vascularNetwork.boundaryConditions[vessel_id] = boundaryInstances
                
                if boundaryIntervalDataDict != {}: 
                    vascularNetwork.boundaryConditionIntervals[vessel_id] = boundaryIntervals
            
            else: netBC[vessel_id] = boundaryInstances
            
    
    #load global Fluid properties
    for my in root.findall(".//globalFluid/my"):
        if str(my.text) != 'None': 
            # find units and convert if needed to SI        
            if update == None: vascularNetwork.globalFluid['my'] = unitConversionToSI(unitsDict, float(eval(my.text)), my.attrib['unit'])
            else: netFluid['my'] = unitConversionToSI(unitsDict,float(eval(my.text)), my.attrib['unit'])
    for rho in root.findall(".//globalFluid/rho"):
        if str(rho.text) != 'None': 
            # find units and convert if needed to SI 
            if update == None: vascularNetwork.globalFluid['rho'] = unitConversionToSI(unitsDict, float(eval(rho.text)), rho.attrib['unit'])
            else: netFluid['rho'] = unitConversionToSI(unitsDict, float(eval(rho.text)), rho.attrib['unit'])
    for gamma in root.findall(".//globalFluid/gamma"):
        if str(gamma.text) != 'None': 
            if update == None: vascularNetwork.globalFluid['gamma'] = float(eval(gamma.text))
            else: netFluid['gamma'] = float(eval(gamma.text))
    for dlt in root.findall(".//globalFluid/dlt"):
        if str(dlt.text) != 'None': 
            if update == None: vascularNetwork.globalFluid['dlt'] = float(eval(dlt.text))
            else: netFluid['dlt'] = float(eval(dlt.text))
        elif str(gamma.text) != 'None': 
            if update == None:  vascularNetwork.globalFluid['dlt'] = (float(eval(gamma.text)) +2.)/(float(eval(gamma.text)) +1.)
            else: netFluid['dlt'] = (float(eval(gamma.text)) +2.0)/(float(eval(gamma.text)) +1.)
                
    for pref in root.findall(".//globalFluid/pref"):
        if str(pref.text) != 'None': 
            if update == None: vascularNetwork.globalFluid['pref'] = unitConversionToSI(unitsDict, float(eval(pref.text)), pref.attrib['unit'])
            else: netFluid['pref'] = unitConversionToSI(unitsDict, float(eval(pref.text)), pref.attrib['unit'])

    # find vessel data
    for vessel_i in root.findall(".//vessel"):
            
            vesselXMLDict = vessel_i.attrib   
            
            vesselID = None
            vesselData = {'name': None, 'start': None, 'end': None, 'leftDaughter': None, 'rightDaughter': None, 'angleToMother': None, 
                          'geom': None, 'length': None, 'lengthInterval': None, 'radiusA': None, 'radiusAInterval': None, 'radiusB': None, 
                          'radiusBInterval': None, 'N': None, 'comp': None, 'Pfunc': None, 'Ps': None, 'As': None, 'wallThickness': None, 
                          'wallThicknessInterval': None, 'youngModulus': None,'youngModulusInterval': None, 'beta': None, 
                          'betaInterval': None,'my': None, 'rho': None, 'gamma': None,'dlt':None , 'pref': None}

            if 'id' in vesselXMLDict and str(vesselXMLDict['id']) != 'None':        
                vesselID = int(vesselXMLDict['id'])
            if 'start_node' in vesselXMLDict and str(vesselXMLDict['start_node']) != 'None':    
                vesselData['start']  = int(vesselXMLDict['start_node'])
            if 'end_node' in vesselXMLDict and str(vesselXMLDict['end_node']) != 'None':      
                vesselData['end']  = int(vesselXMLDict['end_node'])
            if 'name' in vesselXMLDict and str(vesselXMLDict['name']) != 'None':          
                vesselData['name'] = str(vesselXMLDict['name'])
            if 'leftDaughter' in vesselXMLDict and str(vesselXMLDict['leftDaughter']) != 'None':  
                vesselData['leftDaughter']  = int(vesselXMLDict['leftDaughter'])
            if 'rightDaughter' in vesselXMLDict and str(vesselXMLDict['rightDaughter']) != 'None': 
                vesselData['rightDaughter'] = int(vesselXMLDict['rightDaughter'])
            if 'angleToMother' in vesselXMLDict and str(vesselXMLDict['angleToMother']) != 'None': 
                vesselData['angleToMother'] = float(vesselXMLDict['angleToMother'])
            ## old variable names // only version compatibility delete this for next version 0.3
            if 'leftDaugther' in vesselXMLDict and str(vesselXMLDict['leftDaugther']) != 'None':  
                vesselData['leftDaughter']  = int(vesselXMLDict['leftDaugther'])
            if 'rightDaugther' in vesselXMLDict and str(vesselXMLDict['rightDaugther']) != 'None': 
                vesselData['rightDaughter'] = int(vesselXMLDict['rightDaugther'])
            
            # find grid properties
            for grid in vessel_i.findall(".//grid"):
                for data in grid:                     
                    if data.tag == 'geom':
                        if data.text != 'None' and data.text != None: vesselData['geom'] = str(data.text)
                
                    if data.tag == 'length' or data.tag == 'radiusA' or data.tag == 'radiusB' or data.tag == 'N':
                        for value in data.findall(".//scalar"):
                            if value.text != 'None' and value.text !=  None: vesselData[data.tag] = unitConversionToSI(unitsDict, float(eval(value.text)), value.attrib['unit'])
                        for value in data.findall(".//interval"):
                            if value.text != 'None' and value.text !=  None: 
                                interval = []
                                for number in map(float, value.text.split(' ')):
                                    interval.append(unitConversionToSI(unitsDict, float(number), value.attrib['unit']))
                                vesselData[data.tag+'Interval'] =  interval
            
            # find Solid properties  
            for solid in vessel_i.findall(".//solid"):
                for data in solid:   
                                    
                    if data.tag == 'comp':
                        if data.text != 'None' and data.text != None: vesselData['comp']   = str(data.text)
                    if data.tag == 'Pfunc':
                        if data.text != 'None' and data.text != None:
                            if data.text == 'True': data.text = '1'
                            if data.text == 'False': data.text = '0'
                            vesselData['Pfunc']  = bool(int(data.text))
                    if data.tag == 'Ps' or data.tag == 'As':
                        if data.text != 'None' and data.text != None: vesselData[data.tag] = unitConversionToSI(unitsDict, float(eval(data.text)), data.attrib['unit'])
                        
                    if data.tag == 'beta' or data.tag == 'wallThickness' or data.tag == 'youngModulus':
                        for value in data.findall(".//scalar"):
                            if value.text != 'None' and value.text !=  None: vesselData[data.tag] = unitConversionToSI(unitsDict, float(eval(value.text)), value.attrib['unit'])
                        for value in data.findall(".//interval"):
                            if value.text != 'None' and value.text !=  None: 
                                interval = []
                                for number in map(float, value.text.split(' ')):
                                    interval.append(unitConversionToSI(unitsDict, float(number), value.attrib['unit']))
                                vesselData[data.tag+'Interval'] =  interval
                                
            # find fluid properties
            for fluid in vessel_i.findall(".//fluid"):
                for data in fluid:   
                    if data.tag == 'my' or data.tag == 'rho' or data.tag == 'pref' or data.tag == 'gamma' or data.tag == 'dlt':
                        if data.text != 'None' and data.text != None: vesselData[data.tag]   = unitConversionToSI(unitsDict, float(eval(data.text)), data.attrib['unit'])
               
                            
            # pre check if vessels data given is consitent:
            # warning if not
            if vesselData['comp'] == 'Exponential':
                if (vesselData['beta'] or vesselData['Ps'] or vesselData['As']) == None:
                    print "WARNING: not all specifications for %s are met!" %vesselData['comp']
            elif vesselData['comp'] == 'Laplace':
                if (vesselData['wallThickness'] or vesselData['youngModulus'] or vesselData['Ps'] or vesselData['As']) == None:
                    print "WARNING: not all specifications for %s are met!" %vesselData['comp']
            
            
            # calculate area out of radius
            if vesselData['geom'] == 'Uni':
		    if vesselData['radiusA'] != None and vesselData['As'] == None:
		        vesselData['As'] = vesselData['radiusA']**2*pi
		    if vesselData['radiusA'] == None and vesselData['As'] != None:
		        data['radiusA'] = sqrt(data['As']/pi)
	    if vesselData['geom'] == 'Cone':
		    if vesselData['radiusA'] != None and vesselData['radiusB'] != None and vesselData['As'] == None:
		        vesselData['As'] = (vesselData['radiusA']**2+vesselData['radiusB']**2)/2.0*pi
		      
                
            if vesselData['dlt'] == None and vesselData['gamma'] != None:
                vesselData['dlt'] = (vesselData['gamma'] +2.)/(vesselData['gamma'] +1.)
            
            # add vessel to the network
            if update == None: vascularNetwork.addVessel(Id = vesselID, dataDict = vesselData)
            else: vesselDataDict[vesselID] = vesselData
            
    
    if update == None: 
        print " ... loaded successfully"
        return vascularNetwork
    else:
        return {'networkData': networkDataDict, 'vesselData':vesselDataDict}
        
        
def updateNetworkFromXML(filename = None,networkPath = str(cur+"/../NetworkFiles/")):
    
    print ' ... loaded succsessfully updateData'
    updateData = (loadNetworkFromXML(filename=filename, update=True, networkPath = networkPath))
    return updateData
    
