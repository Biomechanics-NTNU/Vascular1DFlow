import csv
from numpy import sqrt,pi

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )

from constants import unitsDict
from constants import bcTags

def writeVesselDataToCSV(vessels,delimiter=';',filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    
    if filename == None:
        filename = "oldNetwork1.csv"
        
    networkDirectory = filename.split('.')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))
    
    tags = ['id','name', 'start', 'end', 'leftDaughter', 'rightDaughter', 'angleToMother', 'geom', 'length','lengthInterval-start',
            'lengthInterval-end', 'radiusA', 'radiusAInterval-start','radiusAInterval-end', 'radiusB', 'radiusBInterval-start',
            'radiusBInterval-end', 'N','NInterval-start','NInterval-end', 'comp', 'Pfunc', 'Ps', 'As', 'wallThickness', 'wallThicknessInterval-start',
            'wallThicknessInterval-end', 'youngModulus','youngModulusInterval-start','youngModulusInterval-end', 'beta', 
            'betaInterval-start','betaInterval-end', 'my', 'rho', 'gamma','dlt', 'pref']
    
    writer = csv.DictWriter(open(networkPath+networkDirectory+'/'+filename,'wb'),tags,delimiter=delimiter)
    
    firstRow = {}
    for item in tags:
        firstRow[item] = item
    
    # write first row
    writer.writerow(firstRow)
    
    #write unit row
    unitRow = {}
    unitRow['id'] = 'unit'
    writer.writerow(unitRow)
    
    data = [] 
    
    for vessel in vessels.itervalues():
        vesselDict = {}
        vesselDict.update( {'id': vessel.Id, 'name': vessel.name, 'start': vessel.start, 'end': vessel.end, 'leftDaughter': vessel.leftDaughter,
                    'rightDaughter': vessel.rightDaughter, 'angleToMother': vessel.angleToMother, 'geom': vessel.geom, 
                    'length': vessel.length, 'radiusA': vessel.radiusA, 'radiusB': vessel.radiusB, 'N': vessel.N, 'comp': vessel.comp, 
                    'Pfunc': vessel.Pfunc, 'Ps': vessel.Ps, 'As': vessel.As,'wallThickness': vessel.wallThickness, 
                    'youngModulus': vessel.youngModulus, 'beta': vessel.beta, 'my': vessel.my, 'rho': vessel.rho, 'gamma': vessel.gamma,
                     'dlt': vessel.dlt, 'pref': vessel.pref})
        if vessel.radiusAInterval != None:
            vesselDict.update({'radiusAInterval-start': vessel.radiusAInterval[0], 'radiusAInterval-end': vessel.radiusAInterval[1]})
        if vessel.radiusBInterval != None:
            vesselDict.update({'radiusBInterval-start': vessel.radiusBInterval[0], 'radiusBInterval-end': vessel.radiusBInterval[1]})
        if vessel.wallThicknessInterval != None:
            vesselDict.update({'wallThicknessInterval-start': vessel.wallThicknessInterval[0],'wallThicknessInterval-end': vessel.wallThicknessInterval[1]})
        if vessel.youngModulusInterval!= None:
            vesselDict.update({'youngModulusInterval-start': vessel.youngModulusInterval[0],'youngModulusInterval-end': vessel.youngModulusInterval[1]})
        if vessel.betaInterval != None:
            vesselDict.update({'betaInterval-start': vessel.betaInterval[0],'betaInterval-end': vessel.betaInterval[1]})
        if vessel.lengthInterval != None:
            vesselDict.update({'lengthInterval-start' : vessel.lengthInterval[0],'lengthInterval-end': vessel.lengthInterval[1]})
        if vessel.NInterval != None:
            vesselDict.update({'NInterval-start' : vessel.NInterval[0],'NInterval-end': vessel.NInterval[1]})
        data.append(vesselDict)
                    
    writer.writerows(data)

def readVesselDataFromCSV(filename= None,delimiter=';',networkPath = str(cur+"/../NetworkFiles/")):
    
    if filename == None:
        filename = "oldNetwork1.csv"
    
    networkDirectory = filename.split('.')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        print ""
        print 'Error: file not in the correct directory / directory does not exists'
        return
    
    if os.path.isfile(networkPath+networkDirectory+'/'+filename) == False:
        print ""
        print "        Error: file does not exists!"
        return {}

    reader = csv.DictReader(open(networkPath+networkDirectory+'/'+filename,'rb'),delimiter=delimiter)

    columUnits = {}
    vesselData = {}
    for row in reader:
        Id = row.pop('id')
        if Id == 'unit':
            columUnits = row
        else:
            Id = int(Id)
            vesselData[Id] = row
    
    strItems   = ['name','geom','comp']
    boolItems  = ['Pfunc']
    floatItems = [ 'angleToMother',  'length', 'radiusA',  'radiusB','Ps', 'As', 'wallThickness', 
                   'youngModulus', 'beta', 'my', 'rho', 'pref','gamma','dlt']
    intItems = ['start', 'end','leftDaughter', 'rightDaughter','N']
    intervallItems = ['NInterval','radiusAInterval','radiusBInterval', 'wallThicknessInterval','youngModulusInterval','betaInterval','lengthInterval']
    
    for data in vesselData.itervalues():
        for item in intervallItems:
            #cast intervall items
            if data[item+'-start'] == '': 
                data[item] = None
            else: 
                # find units and convert if needed to SI
                
                if '#' in data[item+'-start'] or '#' in columUnits[item+'-start']:
                    if '#' not in data[item+'-start']:
                        start = data[item+'-start']
                        nothing,unit = columUnits[item+'-start'].split('#',1)
                    else:
                        start,unit = data[item+'-start'].split('#',1)
                    if ' ' in unit:
                        unit = unit.split(' ')
                        for un in range (0,len(unit),1):
                            start = float(start)*unitsDict[unit[un]]
                    else: start = float(start)*unitsDict[unit]
                else: start = float(data.pop(item+'-start'))
                # find units and convert if needed to SI
                if '#' in data[item+'-end'] or '#' in columUnits[item+'-end']:
                    if '#' not in data[item+'-end']:
                        end = data[item+'-end']
                        nothing,unit = columUnits[item+'-end'].split('#',1)
                    else:
                        end,unit = data[item+'-end'].split('#',1)
                    if ' ' in unit:
                        unit = unit.split(' ')
                        for un in range (0,len(unit),1):
                            end = float(end)*unitsDict[unit[un]]
                    else: end = float(end)*unitsDict[unit]
                else: end = float(data.pop(item+'-end'))
                data[item] = [start,end]
                
                
        for item in strItems:
            if data[item] == '':  
                data[item] = None
            else: data[item] = str(data[item])
        for item in boolItems:
            if data[item] == '':  
                data[item] = None
            else: 
                if data[item] == 'True': data[item] = '1'
                elif data[item] == 'False': data[item] = '0'
                elif data[item] == 'TRUE': data[item] = '1'
                elif data[item] == 'FALSE': data[item] = '0'
                data[item] = bool(int(data[item]))
                
        #cast float items
        for item in floatItems:
            if data[item] == '':  
                data[item] = None
            else:   
                # find units and convert if needed to SI
                if '#' in data[item] or '#' in columUnits[item]:
                    if '#' not in data[item]:
                        value = data[item]
                        nothing,unit = columUnits[item].split('#',1)
                    else:
                        value,unit = data[item].split('#',1)
                    if ' ' in unit:
                        unit = unit.split(' ')
                        for un in range (0,len(unit),1):
                            value = float(value)*unitsDict[unit[un]]
                    else: value = float(value)*unitsDict[unit]
                    data[item] =  value
                    print 'end', item, data[item], unit
                else: data[item] = float(data[item])
                
        # cast int items
        for item in intItems:
            if data[item] == '':  
                data[item] = None
            else:
                # find units and convert if needed to SI
                if '#' in data[item] or '#' in columUnits[item]:
                    if '#' not in data[item]:
                        value = data[item]
                        nothing,unit = columUnits[item].split('#',1)
                    else:
                        value,unit = data[item].split('#',1)
                    if ' ' in unit:
                        unit = unit.split(' ')
                        for un in range (0,len(unit),1):
                            value = float(value)*unitsDict[unit[un]]
                    else: value = float(value)*unitsDict[unit]
                    data[item] = int(value)
                    
                else: data[item] = int(float(data[item]))
        
        # calculate area out of radius
        if data['radiusA'] != None and data['As'] == None:
            data['As'] = data['radiusA']**2*pi
        if data['radiusA'] == None and data['As'] != None:
            data['radiusA'] = sqrt(data['As']/pi)
            
    return {'vesselData':vesselData}    



def writeBCToCSV(boundaryConditions,delimiter=';',filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    
    if filename == None:
        filename = "oldNetwork1BC.csv"
    
    print filename
    networkDirectory = filename.split('BC')[0]
    print networkDirectory
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))
    
    superTag = ['id','type', 'arg1', 'arg2', 'arg3', 'arg4', 'arg5','arg6']
  
    writer = csv.DictWriter(open(networkPath+networkDirectory+'/'+filename,'wb'),superTag,delimiter=delimiter)
    
    firstRow = {}
    emptyRow = {}
    for item in superTag:
        firstRow[item] = item
        emptyRow[item] = ''
    # write first row
    writer.writerow(firstRow)
    writer.writerow(emptyRow)
   
    ## sort BC conditions after types
    typeSort = {}
    for id,BCdicts in boundaryConditions.iteritems():
        print BCdicts
        if BCdicts == {}: BCdicts = {'None':''}
        for type in BCdicts.keys():
            if type not in typeSort.keys() : typeSort[type] = [id]
            else: typeSort[type].append(id)
    
    
    for bcType,ids in typeSort.iteritems():
        args = bcTags[bcType]
        bcTagsRow = {'type':bcType}
        count = 2
        for arg in args:
            bcTagsRow[superTag[count]] = arg
            count = count+1
        writer.writerow(bcTagsRow)          
        #write unit row
        unitRow = {}
        unitRow['id'] = 'unit'
        unitRow['type']= bcType
        writer.writerow(unitRow)
        for id in ids:
            bcValueRow = {}
            if bcType == 'None':
                args = []
                bcValueRow['type']='None'
            else:
                args = boundaryConditions[id][bcType]
                bcValueRow['type']= bcType
            bcValueRow['id'] = id
            count = 2
            for arg in args:
                bcValueRow[superTag[count]] = arg
                count = count+1
            writer.writerow(bcValueRow)     
        writer.writerow(emptyRow)
    
def readBCFromCSV(delimiter=';',filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    
    if filename == None:
        filename = "oldNetwork1BC.csv"
    
    print filename
    networkDirectory = filename.split('BC')[0]
    print networkDirectory
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))
       
    reader = csv.DictReader(open(networkPath+networkDirectory+'/'+filename,'rb'),delimiter=delimiter)
    import pprint as pprint
    
    # create boundary Dictiotionary withs strings
    columUnits = {}
    TypeDescriptions = {}
    BCconditionData = {}
    for row in reader:
        Id = row.pop('id')
        bcType = row.pop('type')
        
        if Id == 'unit':
            columUnits[bcType] = row.values()
        elif Id == '':
            TypeDescriptions[bcType] = row.values()
        else:
            Id = int(Id)
            if bcType != 'None':
                if Id not in BCconditionData.keys(): BCconditionData[Id] = { bcType: row.values()}
                else: BCconditionData[Id].update({bcType: row.values()})
  
    for bcCondition in BCconditionData.itervalues():
        for bcType,args in bcCondition.iteritems():
            
            newArgs = []  
            typeUnits = columUnits[bcType]
            typeDesc = TypeDescriptions[bcType]
            typeDescS = bcTags[bcType]
            print typeDesc
            print typeDescS
            for arg,unitStr,desc,bcTag in zip(args,typeUnits,typeDesc,typeDescS):
                if bcTag == '':break
                print bcTag
                if arg != '':  
                    print arg,unitStr,desc,bcTag
                    # find units and convert if needed to SI
                    if '#' in arg or '#' in unitStr:
                        if '#' not in arg:
                            value = arg
                            nothing,unit = unitStr.split('#',1)
                        else:
                            value,unit = arg.split('#',1)
                        if ' ' in unit:
                            unit = unit.split(' ')
                            for un in range (0,len(unit),1):
                                value = float(value)*unitsDict[unit[un]]
                        else: value = float(value)*unitsDict[unit]
                        newArgs.append(float(value))
                        
                    else: newArgs.append(float(arg))
                else: newArgs.append(None)
                
            if newArgs == []: newArgs.append(None)
            bcCondition[bcType] = newArgs
            
    #pprint.pprint(BCconditionData) 
    
    return BCconditionData  
  