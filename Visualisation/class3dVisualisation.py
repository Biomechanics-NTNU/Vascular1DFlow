import sys,os,gc
from numpy import pi, cos, sin, empty, linspace, ravel, array, hstack, sqrt, delete, append, column_stack, ones
import time

try: # new mayavi version

    from tvtk.api import tvtk
    from mayavi.modules.glyph import Glyph 
    from mayavi.modules.surface import Surface
    from mayavi.sources.vtk_data_source import VTKDataSource
    
    from mayavi.api import Engine
    from pyface.api import GUI
    from pyface.timer.api import Timer

except: #old mayavi version
    from enthought.tvtk.api import tvtk
    from enthought.mayavi.modules.glyph import Glyph 
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    
    from enthought.mayavi.api import Engine
    from enthought.pyface.api import GUI
    from enthought.pyface.timer.api import Timer

# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../NetworkLib')
from classVascularNetwork import VascularNetwork

import pprint
import thread

import subprocess

## Visualisation class for vascular networks, defined with classVascularNetwork

class Visualisation(object):
    ''' Class to visualize vascularNetworks in 3D
        using vtk and mayavi to visualize the vascular network and
        the area over time. The simulation results data as pressure and Flow
        are mapped on the surface of the vessels'''
    
    def __init__(self, vascularNetwork = None):
        '''
        Constructor of the Visualisation
        '''
        #dictionarys to save the src and actors
        self.vizActors = {}
        self.vizSources = {}
        self.vizVTKGrids = {}
        
        self.velArrowGrid = {}
        self.velArrowSrc = {}
        self.glyph = {}
        
        self.vesselLengths = {}
        self.vesselN = {}
        
        # Solution data
        # is numpyarray within numpyarrays 
        # control as 
        # AreaSolution[ vesselID ] [ point in time ]
        # scalarSolution[ vesselID ] [ point in time ]
        self.areaSolution = None
        self.scalarSolution = None
        self.scalarNames = None
        self.scalarPosition = 0
        self.solutionName = None
        self.solutionPosition = 0
        
        # set min max for look uptable
        self.SourceMin = 0
        self.SourceMax = 0
        
        #dimensions | 2 radius (inner outer) | 24 points on circle 
        # | third dimension  N = gridpoints in z axis: not fixed, vessel property
        self.dims = [1, 24]
    
        # vascularNetwork
        self.vascularNetwork = vascularNetwork 
        
        # engine and scene for visualisation
        self.engine = Engine()
        self.scene = None
        
        # to increase the radius change
        # radius(t) = radius0 + deltaRadius * factor
        self.factor = 50
        
        # Movie 
        self.movieBool = False 
        self.movNumber = 0
        self.movPath = None
        
        # update control
        self.moveWall = True
        self.moveWallOff = False
        self.pauseBool = True
        self.calls = 0
        self.currentTime = 0
        self.endTime = None 
        self.updateTimeStep = 1
        
        # bool for velocity profile
        self.vecBool = False
        self.powerLawCoefficient = 2.0
        
        self.static = False
        
        # 2 d plot figures
        self.figspace = None
        self.figtime = None
        self.plotSpace = False
        self.picked = None
        
        #### init this stuff prob move to another function
        self.initialize()

    def __del__(self):  
        '''
        Class Destructor
        '''
        del(self.vascularNetwork)
        print " visualisation removed"
        
#-------------------## public methods##-------------------#

    
    def initialize(self):
        '''
        Initialize and create a 3D representation of the network
        '''
        if self.vascularNetwork != None:
            
            gc.enable()
            
            for key,vessel in self.vascularNetwork.vessels.iteritems():
                self.vesselLengths[key] = vessel.length
                self.vesselN[key] = vessel.N
            
            self.scene = self.engine.new_scene(name = self.vascularNetwork.name)
            
            self.engine.scenes[0].name = str(self.vascularNetwork.name+' - '+'paused')
            
            self.scene.scene.background = (0.0,0.0,0.0)
            
            # start the engine
            self.engine.start()
            # set cameraposition
            self.scene.scene.camera.azimuth(90)
            self.scene.scene.camera.roll(90)

            # create the initial network
            self.createVizData()
            
            # add the initial network
            self.scene.scene.disable_render = True
            for Id in self.vizSources.iterkeys():
                self.engine.add_source(self.vizSources[Id],self.scene)
                self.engine.add_module(self.vizActors[Id], obj=self.vizSources[Id])
                self.engine.add_source(self.velArrowSrc[Id], self.scene)
                self.glyph[Id].actor.actor.visibility = self.vecBool # set bool for visualize the arrows
                self.engine.add_module(self.glyph[Id], obj=self.velArrowSrc[Id])
            # refit the camera position
            self.scene.scene.reset_zoom()
            self.scene.scene.disable_render = False
            
            if self.static == False:
                # add interactor
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_N)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_M)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_B)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_V)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_C)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_X)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_Y)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_J)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_H)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_U)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_G)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_komma)
                self.scene.scene.interactor.add_observer("KeyPressEvent",self.key_press_punkt)
                #self.scene.mayavi_scene.on_mouse_pick(self.picker_callback)
                self.engine.scenes[0].on_mouse_pick(self.picker_callback)

            print " ..3d visualisation initialized!"
           
    
    def setVascularNetwork(self,vascularNetwork):
        '''
        Set the vascular Network if not set in constructor
        and initialize the network
        '''
        self.vascularNetwork = vascularNetwork
        self.initialize() 
        
    def setSolutionData(self,solutionData):
        '''
        Process the solutionData, calculates the velocity and set all needed varibles.
        solutionData is a list with dictionarys containing the solution one or several simulation solutions.
           solutionData = [ {'Pressure': pressureData, 'Flow': flowData, 'Area': areaData, 'Name': 'solution1' }, ...]
        '''
        #initialize
        self.areaSolution = []
        self.scalarSolution = []
        self.scalarNames = []
        self.solutionName = []
        
        for solution in solutionData:
            
            # save the solution name
            currentSolutionName = solution.pop('Name')
            
            #save the scalars of the solution
            currentScalars = []
            currentScalarNames = []
            
            # transfer pressure to mmHg
            Psol = {}
            for vesselId,Ptemp in solution['Pressure'].iteritems():
            	Psol[vesselId] = Ptemp/133.32
            
            currentScalars.append(Psol)
            currentScalarNames.append('Pressure')
            currentScalars.append(solution.pop('Flow'))
            currentScalarNames.append('Flow')
            currentScalars.append(solution.pop('Area'))
            currentScalarNames.append('Area')
            
            # calculate velocities 
            velocity = []
            for i in range (0,len(currentScalars[1]),1):
                vel_t = array([currentScalars[1][i][0]/currentScalars[2][i][0]])
                                
                for t in range(1,len(currentScalars[1][i]),1):
                    v = currentScalars[1][i][t]/currentScalars[2][i][t]
                    vel_t = append(vel_t,[v],axis = 0)
                velocity.append(vel_t)
            
            # set the velocity    
            currentScalars.append(velocity)
            currentScalarNames.append('Velocity')
            
            # check if additional solution data was passed:
            for key,value in solution.iteritems():
                currentScalars.append(value)
                currentScalarNames.append(key)
            
            # save the data in class variables
            self.areaSolution.append(currentScalars[2])
            self.scalarSolution.append(currentScalars)
            self.scalarNames.append(currentScalarNames)
            self.solutionName.append(currentSolutionName)
        
        del(currentScalars)
        del(currentScalarNames)
        del(currentSolutionName)

        self.endTime = len(self.areaSolution[self.solutionPosition][0])
        self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'paused'+' - '+self.solutionName[self.solutionPosition])
            
        print ' .. solution data initialized!'
        
        
    def visualize(self):
        '''
        Start the visualisation of the simulation result.
        Set up Upgrade method and starts the main loop.
        Use only if the simulation result was set!
        '''
        # update the lookuptables and set them
        minL = []
        maxL = []
       
        for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
            minL.append(array.min())
            maxL.append(array.max())
        self.SourceMin = min(minL)
        self.SourceMax = max(maxL)

        for actor in self.vizActors.itervalues():
            actor.module_manager.scalar_lut_manager.data_range = [self.SourceMin, self.SourceMax]
        #self.vizActors[0].module_manager.scalar_lut_manager.data_name = self.scalarNames[self.solutionPosition][self.scalarPosition]
        #self.vizActors[0].module_manager.scalar_lut_manager.show_scalar_bar = True

        timer = Timer(0.03, self.update)
        gui = GUI()
        gui.busy = True
        gui.start_event_loop()

    def visualizeStatic(self):
        gui = GUI()
        self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'static')
        gui.start_event_loop()        
    
            
    def picker_callback(self, picker):
        '''
        Mouse pick callback, opens a 2D time plot if a vessel is picked; If the plotSpace is activated key "g"
        the animated space plot of the vessel is loaded aswell
        '''
        printOut = ' .. world - picked'
        if len(picker.actors) > 0:
            picked = picker.actors[0]
            for key,value in self.vizActors.iteritems():
                if value.actor.actor._vtk_obj == picked._vtk_obj:
                    
                    self.picked = self.vascularNetwork.vessels[key].name
                    printOut = str(' .. '+self.vascularNetwork.name +' - '+ self.picked + ' - picked')
                    
                    string = ' '.join(['python',cur+'/class2dVisualisationSimple.py','-f',self.vascularNetwork.name, '-n 1','-i',str(key)])                
                    subprocess.call(string, shell=True)
#                    if self.static == False:
#                        
#                        if self.figspace is not None and type(self.figspace) is not int: 
#                            self.figspace.clf()
#                        if self.figtime is not None:
#                            self.figtime.clf()
#        
#                        networkPlot = NetworkPlot([self.scalarSolution[self.solutionPosition][0][key]],
#                                                  [self.scalarSolution[self.solutionPosition][1][key]],
#                                                  [self.scalarSolution[self.solutionPosition][2][key]])
#                        
#                        #self.figspace,self.figtime =  thread.start_new_thread(networkPlot, ({name: self.vascularNetwork.getSimulationDictionary()[name]},
#                        #                               0,1,False,True,True,self.vascularNetwork.simulationContext['totalTime'],True))
#                        self.figspace,self.figtime = networkPlot(VesselNetwork = {self.picked: self.vascularNetwork.getSimulationDictionary()[key]},
#                                                        cutTime=0,CUT=1,plotVelocity=False,plotSpace=self.plotSpace,plotTime=True,
#                                                       totaltime = self.vascularNetwork.simulationContext['totalTime'], scaleplot=True)                 
        print printOut
            
    def key_press_G(self,obj,event):
        '''
        Key event G: enable / disable 2D space plot if vessel is picked
        '''
        key = obj.GetKeyCode()
        if key=='G':
            if self.plotSpace  == False: 
                self.plotSpace  = True
                print " .. 2D space plot enabled"
            else: 
                self.plotSpace = False
                print " .. 2D space plot disabled"
        elif key=='g':
            if self.figtime is not None:
                print " .. time plot for ",self.picked," saved"
                path = str(cur+'/../network_files/'+self.vascularNetwork.name+'/2dTimePlots/')
                if not os.path.exists(path):
                    os.makedirs(path) 
                self.figtime.savefig(str(path+self.vascularNetwork.name+'-'+self.solutionName[self.solutionPosition]+'-'+self.picked+'.png'))

    def key_press_N(self,obj,event):
        '''
        Key event n: next lookup table - change the lookup and scalar mapping
        '''
        key = obj.GetKeyCode()
        if key=='n' or key=='N':
            self.scalarPosition = self.scalarPosition + 1
            if self.scalarPosition == len(self.scalarNames[self.solutionPosition]):
                self.scalarPosition = 0
            minL = []
            maxL = []
            for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
                minL.append(array.min())
                maxL.append(array.max())
            self.SourceMin = min(minL)
            self.SourceMax = max(maxL)
            
            
    def key_press_J(self,obj,event):
        '''
        Key event j: last lookup table - change the lookup and scalar mapping
        '''
        key = obj.GetKeyCode()
        if key=='j' or key=='J':
            self.scalarPosition = self.scalarPosition - 1
            if self.scalarPosition < 0:
                self.scalarPosition = len(self.scalarNames[self.solutionPosition]) - 1
            # recalculate the range of the lookuptable ## Note: change this in further version to save memory -> (list)
            minL = []
            maxL = []
            for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
                minL.append(array.min())
                maxL.append(array.max())
            self.SourceMin = min(minL)
            self.SourceMax = max(maxL)
            
            
    def key_press_H(self,obj,event):
        '''
        Key event h: next solution data - show next soultion data set
        '''
        key = obj.GetKeyCode()
        if key=='h' or key=='H':
            self.solutionPosition = self.solutionPosition + 1
            if self.solutionPosition == len(self.solutionName):
                self.solutionPosition = 0
            # recalculate the range of the lookuptable ## 
            self.scalarPosition = 0
            minL = []
            maxL = []
            for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
                minL.append(array.min())
                maxL.append(array.max())
            self.SourceMin = min(minL)
            self.SourceMax = max(maxL)
    
            self.endTime = len(self.areaSolution[self.solutionPosition][0])
            
            #set name of the solution
            if self.pauseBool == False: self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'paused'+' - '+self.solutionName[self.solutionPosition])
            else: self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'running'+' - '+'t = 0.000' +' - '+self.solutionName[self.solutionPosition])
    
            #restart simulation
            self.currentTime = 0
            
    def key_press_U(self,obj,event):
        '''
        Key event u: last solution data - show last soultion data set
        '''
        key = obj.GetKeyCode()
        if key=='u' or key=='U':
            self.solutionPosition = self.solutionPosition - 1
            if self.solutionPosition < 0:
                self.solutionPosition = len(self.solutionName) - 1
            # recalculate the range of the lookuptable ## 
            self.scalarPosition = 0
            minL = []
            maxL = []
            for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
                minL.append(array.min())
                maxL.append(array.max())
            self.SourceMin = min(minL)
            self.SourceMax = max(maxL)
    
            #set time
            self.endTime = len(self.areaSolution[self.solutionPosition][0])
            
            #set name of the solution
            if self.pauseBool == False: self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'paused'+' - '+self.solutionName[self.solutionPosition])
            else: self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'running'+' - '+'t = 0.000' +' - '+self.solutionName[self.solutionPosition])
    
            #restart simulation
            self.currentTime = 0
    
            
    def key_press_M(self,obj,event):
        '''
        Key event m: start saving png for movie
        '''
        key = obj.GetKeyCode()
        if key=='m' or key=='M':
            # create movie directory
            if self.movPath == None:
                self.movPath = str(cur+'/../NetworkFiles/'+self.vascularNetwork.name+'/movieTemplateData/')
                if not os.path.exists(str(self.movPath)):
                    os.makedirs(str(self.movPath)) 
            
            if self.movieBool == False: 
                self.movieBool = True
                print " .. start creating movie files"
            else: 
                self.movieBool = False
                self.movNumber = 0
                print " .. stop creating movie files"
        
    def key_press_B(self,obj,event):
        '''
        Key event B: pause/break the simulation
        '''
        key = obj.GetKeyCode()
        if key=='b' or key=='B':
            if self.pauseBool == False: 
                self.pauseBool = True
                self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'paused'+' - '+self.solutionName[self.solutionPosition])
            else: 
                self.pauseBool = False
                self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'running'+' - '+'t = 0.000' +' - '+self.solutionName[self.solutionPosition])
    
    
    def key_press_V(self,obj,event):
        '''
        Key event V: open vessels and show vectorfields
        '''
        key = obj.GetKeyCode()
        if key=='v' or key=='V':
            if self.vecBool == False: 
                self.vecBool = True
                for glyph in self.glyph.itervalues(): glyph.actor.actor.visibility = self.vecBool
                self.scalarPosition = 3
                
                minL = []
                maxL = []
                for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
                    minL.append(array.min())
                    maxL.append(array.max())
                    
                self.SourceMin = min(minL)
                self.SourceMax = max(maxL)
                
                self.vizActors[0].module_manager.scalar_lut_manager.show_scalar_bar = False
                self.glyph[0].module_manager.scalar_lut_manager.data_range = [0, self.SourceMax]
                self.glyph[0].module_manager.scalar_lut_manager.show_scalar_bar = True
                self.glyph[0].module_manager.scalar_lut_manager.data_name = self.scalarNames[self.solutionPosition][self.scalarPosition]
                
                print " .. velocity profile enabled"
            else: 
                self.vecBool = False
                for glyph in self.glyph.itervalues(): glyph.actor.actor.visibility = self.vecBool
                print " .. velocity profile disabled"
                self.scalarPosition = 3
                minL = []
                maxL = []
                for array in self.scalarSolution[self.solutionPosition][self.scalarPosition].itervalues(): 
                    minL.append(array.min())
                    maxL.append(array.max())
                self.SourceMin = min(minL)
                self.SourceMax = max(maxL)
                self.vizActors[0].module_manager.scalar_lut_manager.data_name = self.scalarNames[self.solutionPosition][3]
                self.vizActors[0].module_manager.scalar_lut_manager.show_scalar_bar = True
                
                
    def key_press_C(self,obj,event):
        '''
        Key event C: stop/start simulation of the wall movement
        '''
        key = obj.GetKeyCode()
        print key
        if key=='c' or key=='C':
            if self.moveWall  == False: 
                self.moveWall  = True
                print " .. wall movement enabled"
            else: 
                self.moveWallOff = True
                print " .. wall movement disabled"
                
    def key_press_X(self,obj,event):
        '''
        Key event X: increase the radius factor
        '''
        key = obj.GetKeyCode()
        if key=='x' or key=='X':
            self.factor = self.factor + 1
            print " .. radius factor increased to: ",self.factor
    
    def key_press_Y(self,obj,event):
        '''
        Key event Y: decrease the radius factor
        '''
        key = obj.GetKeyCode()
        if key=='y' or key=='Y' or key=='z' or key=='Z':
            if self.factor > 0:
                self.factor = self.factor - 1
            print " .. radius factor decreased to: ",self.factor
    
    def key_press_komma(self,obj,event):
        '''
        Key event ,: increase the timestep
        '''
        key = obj.GetKeyCode()
        if key==',' or key==',':
            self.updateTimeStep = self.updateTimeStep + 1
            print " .. update timestep increased to: ",self.updateTimeStep
    
    def key_press_punkt(self,obj,event):
        '''
        Key event ,: decrease the timestep
        '''
        key = obj.GetKeyCode()
        if key=='.' or key=='.':
            if self.updateTimeStep > 1:
                self.updateTimeStep = self.updateTimeStep - 1
            print " .. update timestep decreased to: ",self.updateTimeStep
            
    

#-------------------## private methods##-------------------#

    def update(self):
        '''
        update function: updates the scene,actors every time it is called
        '''
        if self.endTime == None:
            self.endTime = len(self.areaSolution[self.solutionPosition][0])
            
        if self.pauseBool == False:
            
            
            Lstart = time.clock()
            
            ###start the update  
            # stop the rendering while updating grids and the piplines
            # without mayavi is rendering during the update 
            # --> updateTime is 10x higher
            # --> all vessels are updated at the same time
            self.scene.scene.disable_render = True
            
            # to accelerate fist calculat data then apply
            TempColor = {}
            TempPoints = {}
            TempVelPoints = {}
            TempVelColor = {}
            #create new data
            LstartD = time.clock()
            for Id in self.vesselLengths.iterkeys():
                
                if self.moveWall == True:
                    
                    radiusAr = sqrt(self.areaSolution[self.solutionPosition][Id][0]/pi)+(sqrt(self.areaSolution[self.solutionPosition][Id][self.currentTime]/pi)-sqrt(self.areaSolution[self.solutionPosition][Id][0]/pi))*self.factor # calculate the radius from the area    
                    if self.moveWallOff == True:
                        radiusAr = sqrt(self.areaSolution[self.solutionPosition][Id][0]/pi)
                    # create points for the vessels
                    TempPoints[Id] = self.generate_vertices(r_array = radiusAr, length = self.vesselLengths[Id], N = self.vesselN[Id])        
                    # create points for the velocity arrows
                    TempVelPoints[Id] = self.generate_velocity_points(length = self.vesselLengths[Id], N = self.vesselN[Id], r_array= radiusAr)
                    
                # set size and color of the velocity arrows 
                if self.vecBool == True:
                    ## here the 3 must refer to the velocity!! in scalarSolution 
                    self.scalarPosition = 3
                    TempVelColor[Id] = self.generate_velocity_colorScalars(velArray = self.scalarSolution[self.solutionPosition][self.scalarPosition][Id][self.currentTime], N = self.vesselN[Id])
                    
                    if Id == 1:
                        print self.scalarSolution[self.solutionPosition][self.scalarPosition][Id][self.currentTime][1]
                    #set the scalar color array to velocity
                    colorAr = self.scalarSolution[self.solutionPosition][self.scalarPosition][Id][0]  
                else: 
                    # set the scalar color array to the current choosen
                    colorAr = self.scalarSolution[self.solutionPosition][self.scalarPosition][Id][self.currentTime] 
                # create scalar color map
                TempColor[Id]= ravel(self.generate_colorScalars(colorAr, N = self.vesselN[Id]))  
            
            LendD = time.clock()
            #print ' time:', LendD - LstartD
            
            for Id,grid in self.vizVTKGrids.iteritems():
                # modify grid of the vizVTKGrid
                # create scalar color map
                grid.point_data.scalars = TempColor[Id]
                grid.point_data.scalars.name = 'scalars' 
                # create update data if wall movement is enabled
                if self.moveWall == True:
                    grid.points = TempPoints[Id]
                    
                # set lut range (necessary)
                self.vizActors[Id].module_manager.scalar_lut_manager.data_range = [self.SourceMin, self.SourceMax]
            # set lut (for 1 actor is enough
            
            
            
            # update velocity vectors
            if self.vecBool == True:
                for Id,grid in self.velArrowGrid.iteritems():
                    grid.point_data.scalars = TempVelColor[Id][1]
                    grid.point_data.scalars.name = 'scalars' 
                    grid.point_data.vectors = TempVelColor[Id][0]
                    grid.point_data.vectors.name = 'vectors'
                    
                    self.glyph[Id].glyph.scale_mode = 'scale_by_scalar'
                    self.glyph[Id].glyph.glyph.range = [0, self.SourceMax]     
                    self.glyph[Id].glyph.color_mode = 'color_by_scalar'       
                    self.glyph[Id].module_manager.scalar_lut_manager.data_range = [0, self.SourceMax]
                    self.glyph[Id].module_manager.scalar_lut_manager.use_default_range = False
                    self.glyph[0].module_manager.scalar_lut_manager.show_scalar_bar = True
                    self.glyph[0].module_manager.scalar_lut_manager.data_name = self.scalarNames[self.solutionPosition][self.scalarPosition]
                    
                    # update position if they change (wallmovement)
                    if self.moveWall == True:
                        grid.points = TempVelPoints[Id]
            #else:
                #self.vizActors[0].module_manager.scalar_lut_manager.data_name = self.scalarNames[self.solutionPosition][self.scalarPosition]
                #self.vizActors[0].module_manager.scalar_lut_manager.show_scalar_bar = True
            
            if self.moveWallOff == True:
                self.moveWall = False
                self.moveWallOff = False
            
            dt = (self.vascularNetwork.simulationContext['totalTime']/self.endTime)*self.currentTime
            self.engine.scenes[0].name = str(self.vascularNetwork.name +' - '+'running'+' - '+'t = %0.4f' %dt+' - '+self.solutionName[self.solutionPosition])
    
            # control the current timeStep
            self.currentTime = self.currentTime+self.updateTimeStep
            # restart simulation if end reached
            if self.currentTime >= self.endTime:
                self.currentTime = 0
                
            # if no movie:
            if self.movieBool == False:
                #start the rendering again
                self.scene.scene.disable_render = False
            #save pictures for the movies
            else: 
                #thread.start_new_thread(self.scene.scene.save_png,(('../data/temp/'+'movieTemplateData'+str(self.movNumber).zfill(4)+'.png'),  ))
                Lstart = time.clock()
                self.scene.scene.save_png(self.movPath+'movieTemplateData'+str(self.movNumber).zfill(4)+'.png')
                Lend = time.clock()
                print" needed %1.6f  sec to save the picture" %((Lend-Lstart)) 
                self.movNumber = self.movNumber + 1
            
            Lend = time.clock()
            #print" update time %1.6f -- framerate: %1.2f -- increase in framerate: %1.2f " %((Lend-Lstart), 1.0/(Lend-Lstart), (1.0/(Lend-Lstart-LendD+LstartD))-(1.0/(Lend-Lstart))) 
            
 
    def generate_vertices(self,length,N,r_array = None,r_init = None):
        '''
        Calculates vertices of an vessel
        either based on radius array with lenght N
        or radius init = [radius start, radius end]
        '''
        N = int(N)
        if r_array == None:
            r_initA = r_init[0]
            r_initB = r_init[1]
            r_array = linspace(r_initA, r_initB, N)
        
        # inner wall (not necesarry!)   
        #ri_array = r_array*array([0.9])
        #r = vstack((ri_array,r_array)).transpose()
        r = r_array.transpose()
        
        if self.vecBool == False : 
            thetaFactor = 2
        else: 
            thetaFactor = 1
        
        theta = linspace(pi/2.0, pi/2.0+thetaFactor*pi, self.dims[1])
        z = linspace(0, length, N)
        
        ## create Points
        aLength = self.dims[0]*self.dims[1]
        points = empty([aLength*N,3])
        
        start = 0
        for z_plane in range(0,len(z),1):
            end = start+aLength
            plane_points = points[start:end]
            plane_points[:,0] = (cos(theta)*r[z_plane]).ravel() # for inner wall add r[z_plane,:][:,None]
            plane_points[:,1] = (sin(theta)*r[z_plane]).ravel()
            plane_points[:,2] = z[z_plane]
            start = end
        return points

    def generate_colorScalars(self,colorArray, N):
        '''
        map the values of the colorArray to the points of the grid
        '''
        # for loop
        #colors = array([])
        #for z_plane in range(0,int(N),1):
        #    s_z = linspace(colorArray[z_plane] ,colorArray[z_plane], (self.dims[0]*self.dims[1]))
        #    colors = hstack((colors,s_z))
        
        # matrix calculation
        colors = column_stack(ones(((self.dims[0]*self.dims[1]),int(N)))*colorArray).ravel()
       
        return colors
    
    def generate_velocity_points(self, length, N , r_array = None,r_init = None, numberOfVec = 11):
        '''
        calculates the points for the velocity arrows of the velocity profil
        '''
        #numberOfVec should be odd
        
        if r_array == None:
            r_initA = r_init[0]
            r_initB = r_init[1]
            r_array = linspace(r_initA, r_initB, N)
        
        dz = length/N
        z = linspace(dz, length-dz, (N-2))
        
        ## create Points
        aLength = numberOfVec
        points = empty([aLength*(N-2),3])
        
        start = 0
        for z_plane in range(0,len(z),1):
            end = start+aLength
            y_plane = linspace(- r_array[z_plane]*0.9, r_array[z_plane]*0.9, numberOfVec)
            plane_points = points[start:end]
            plane_points[:,0] = 0
            plane_points[:,1] = y_plane
            plane_points[:,2] = z[z_plane]
            start = end
        return points
    
    def generate_velocity_colorScalars(self, N, velArray = None, numberOfVec=11):
        '''
        Map the values of the velocity calculated of the mean velocity to the points of the grid
        '''
        
        if velArray == None:
            velC = linspace(0.001,0.005,int((N-2)))
        else:
            velC = delete(velArray,[0,int(N-1)])
                  
        profileA = linspace(0.9,0.0,6)
        profileB = linspace(0,0.9,6)
        prof = append(profileA,profileB)
        profile = delete(prof, 6)
        
        #cast N
        colorVectors = empty([numberOfVec*(N-2),3])
        colorScalars = array([])
        count = 0
        for z in range(0,int(N-2),1):
            for i in range (0,numberOfVec,1):
                #velocity profile using powerlaw with n = self.powerLawCoefficient
                vel = velC[z]*((self.powerLawCoefficient+2.0)/self.powerLawCoefficient)*(1.0-profile[i]**self.powerLawCoefficient)
                
                # vectors are for direction of the arrows
                t = 1.0
                if vel != 0:
                    t = vel/abs(vel)
                # scalares determine size and color
                vec_i = array([0,0,t*abs(vel)])
                colorVectors[count] = vec_i
                colorScalars = hstack((colorScalars,[abs(vel)]))
                count = count + 1
            
                    
        
        return [colorVectors,colorScalars]
            
    
    def createVizData(self):
        '''
        Creates the actors and sources of the vascular network
        parsing through the network as binary tree
        using connection of mother and daughter vessels
        '''
        ### INIT
        viz = []
        vessels = self.vascularNetwork.vessels
        root = self.vascularNetwork.root[0]
        
        ### root vessel
        
        ## find data
        #find initial rotation
        rootRot = 0.0
        if vessels[root].angleToMother != None:
            rootRot = vessels[root].angleToMother 
        
        #find inital RadiusA
        dRadiusA = 0.05
        if vessels[root].radiusA != None:
            dRadiusA = vessels[root].radiusA
            
        #find inital RadiusB
        dRadiusB = dRadiusA
        if vessels[root].radiusB != None:
            dRadiusB = vessels[root].radiusB
        
        #find length of root
        rLength = 1.0
        if vessels[root].length != None:
            rLength = vessels[root].length
        
        ## create visualisation Data
        # create the data points 
        points =  self.generate_vertices(r_init=[dRadiusA,dRadiusB],length = rLength, N = self.vesselN[root])
        # create scalar color map
        col = linspace(0.0, 0.0, self.vesselN[root])
        color = self.generate_colorScalars(col, N = self.vesselN[root])
        # generate structured grid, data source
        self.vizVTKGrids[root] = tvtk.StructuredGrid(dimensions=(self.dims[1], self.dims[0], self.vesselN[root]))
        #set data points
        self.vizVTKGrids[root].points = points
        #set scalar color map
        self.vizVTKGrids[root].point_data.scalars = ravel(color)
        self.vizVTKGrids[root].point_data.scalars.name = 'scalars' #self.scalarNames[self.scalarPosition]
        
        # create datasource for mayavi
        self.vizSources[root] = VTKDataSource(data = self.vizVTKGrids[root])
  
        # create surface actor for the datas ource
        self.vizActors[root] = Surface()
        self.vizActors[root].actor.actor.position=0,0,0
        self.vizActors[root].actor.actor.rotate_x(rootRot)
    
        ### create the velocity vector field 
        # creat gird for the velocity arrows
        self.velArrowGrid[root] = tvtk.StructuredGrid(dimensions=(11, 1, self.vesselN[root]))
        #set data points
        self.velArrowGrid[root].points = self.generate_velocity_points(self.vesselLengths[root], self.vesselN[root], r_init=[dRadiusA,dRadiusB])
        #set scalar color map
        color = self.generate_velocity_colorScalars(self.vesselN[root])
        self.velArrowGrid[root].point_data.vectors = color[0].copy()
        self.velArrowGrid[root].point_data.vectors.name = 'vectors'
        self.velArrowGrid[root].point_data.scalars = color[1].copy()
        self.velArrowGrid[root].point_data.scalars.name = 'scalars'
        # create source
        self.velArrowSrc[root] = VTKDataSource(data = self.velArrowGrid[root])
        # create vectors data
        self.glyph[root] = Glyph()
        self.glyph[root].glyph.scale_mode = 'scale_by_scalar'
        self.glyph[root].glyph.color_mode = 'color_by_scalar'
        
        self.glyph[root].glyph.glyph_source.glyph_source = self.glyph[root].glyph.glyph_source.glyph_dict['arrow_source']
        self.glyph[root].glyph.glyph_source.glyph_position = 'tail'
        self.glyph[root].glyph.glyph.scale_factor = self.vesselLengths[root] / self.vesselN[root] * 4 / 7
        
        
        # set position
        self.glyph[root].actor.actor.position = 0,0,0
        self.glyph[root].actor.actor.rotate_x(rootRot)
    
        viz.append(root)
    
        ## set rest of the network
        while len(viz) != 0:
        # Current left Daughter
            currentVessel = viz.pop(0)
    
            #find values of the current mother
            moPos = self.vizActors[currentVessel].actor.actor.position
            moRot = self.vizActors[currentVessel].actor.actor.orientation[0]
            moRotY = self.vizActors[currentVessel].actor.actor.orientation[1]
            
            moLength = 1.0
            if vessels[currentVessel].length != None:
                moLength = vessels[currentVessel].length
            if vessels[currentVessel].radiusB != None:
                moRadiusB = vessels[currentVessel].radiusB
            elif vessels[currentVessel].radiusA != None:
                moRadiusB = vessels[currentVessel].radiusA
            else:
                moRadiusB = 0.05
            
            #find Daughters
            rightDaughter = None
            rightDaughter = vessels[currentVessel].rightDaughter 
            leftDaughter = None
            leftDaughter = vessels[currentVessel].leftDaughter 
            
            # create left Daughter vessel visualisation
            if leftDaughter is not None:
                        
                #find length of Daughter
                dLength = 1.0
                if vessels[leftDaughter].length != None:
                    dLength = vessels[leftDaughter].length
                
                #find inital RadiusA
                dRadiusA = moRadiusB
                if vessels[leftDaughter].radiusA != None:
                    dRadiusA = vessels[leftDaughter].radiusA
                
                #find inital RadiusB
                dRadiusB = dRadiusA
                if vessels[leftDaughter].radiusB != None:
                    dRadiusB = vessels[leftDaughter].radiusB
                
                #set rotation
                
                if rightDaughter is not None:
                    if vessels[leftDaughter].angleToMother is not None:
                        drot =  -moRotY+(1+2/180.0*moRotY)*((1+2/180.0*moRotY)*-vessels[leftDaughter].angleToMother + moRot) 
                    else:
                        drot = -moRotY+(1+2/180.0*moRotY)*((1+2/180.0*moRotY)*-30.0 + moRot)
                else:
                    drot = moRot
                
                #set z_position
                pos_z = moLength*cos(moRot*pi/180)*(1+2/180.0*moRotY)
                #set y_position
                if moRot != 0.0:
                    pos_y = -moLength*sin(moRot*pi/180)
                else:
                    pos_y = 0.0
                # apply position changes
                dpos = moPos[0],moPos[1]+pos_y,moPos[2]+pos_z
                
                ## create visualisation Data
                # create the data points 
                points =  self.generate_vertices(r_init=[dRadiusA,dRadiusB],length = dLength,N = self.vesselN[leftDaughter])
                # create scalar color map
                col = linspace(0.0, 0.0, self.vesselN[leftDaughter])
                color = self.generate_colorScalars(col, N = self.vesselN[leftDaughter])
                # generate structured grid, data source
                self.vizVTKGrids[leftDaughter] = tvtk.StructuredGrid(dimensions=(self.dims[1], self.dims[0], self.vesselN[leftDaughter]))
                #set data points
                self.vizVTKGrids[leftDaughter].points = points
                #set scalar color map
                self.vizVTKGrids[leftDaughter].point_data.scalars = ravel(color.copy())
                self.vizVTKGrids[leftDaughter].point_data.scalars.name = 'scalars' #self.scalarNames[self.scalarPosition]
                
                # create datasource for mayavi
                self.vizSources[leftDaughter] = VTKDataSource(data = self.vizVTKGrids[leftDaughter])
                
                # create surface actor for the datas ource
                self.vizActors[leftDaughter] = Surface()
                self.vizActors[leftDaughter].actor.actor.position = dpos
                self.vizActors[leftDaughter].actor.actor.rotate_x(drot)
                
                ## set arrows
                self.velArrowGrid[leftDaughter] = tvtk.StructuredGrid(dimensions=(11, 1, self.vesselN[leftDaughter]))
                #set data points
                self.velArrowGrid[leftDaughter].points = self.generate_velocity_points(self.vesselLengths[leftDaughter], self.vesselN[leftDaughter], r_init=[dRadiusA,dRadiusB])
                #set scalar color map
                color = self.generate_velocity_colorScalars(self.vesselN[leftDaughter])
                self.velArrowGrid[leftDaughter].point_data.vectors = color[0].copy()
                self.velArrowGrid[leftDaughter].point_data.vectors.name = 'vectors'
                self.velArrowGrid[leftDaughter].point_data.scalars = color[1].copy()
                self.velArrowGrid[leftDaughter].point_data.scalars.name = 'scalars'
                # create source
                self.velArrowSrc[leftDaughter] = VTKDataSource(data = self.velArrowGrid[leftDaughter])
                # create vectors data
                self.glyph[leftDaughter] = Glyph()
                self.glyph[leftDaughter].glyph.scale_mode = 'scale_by_scalar'
                self.glyph[leftDaughter].glyph.color_mode = 'color_by_scalar'
                self.glyph[leftDaughter].glyph.glyph_source.glyph_source = self.glyph[leftDaughter].glyph.glyph_source.glyph_dict['arrow_source']
                self.glyph[leftDaughter].glyph.glyph_source.glyph_position = 'tail'
                self.glyph[leftDaughter].glyph.glyph.scale_factor = self.vesselLengths[leftDaughter] / self.vesselN[leftDaughter] * 4 / 7
                
                # set position
                self.glyph[leftDaughter].actor.actor.position = dpos
                self.glyph[leftDaughter].actor.actor.rotate_x(drot)
                
                
                #check for children
                if vessels[leftDaughter].leftDaughter is not None:
                    viz.append(leftDaughter)
                
                # right Daughter (only if left Daughter)
                if rightDaughter is not None:
                    
                    #find length of Daughter
                    dLength = 1.0
                    if vessels[rightDaughter].length != None:
                        dLength = vessels[rightDaughter].length
                    
                    #find inital RadiusA
                    dRadiusA = moRadiusB
                    if vessels[rightDaughter].radiusA != None:
                        dRadiusA = vessels[rightDaughter].radiusA
 
                    #find inital RadiusB
                    dRadiusB = dRadiusA
                    if vessels[rightDaughter].radiusB != None:
                        dRadiusB = vessels[rightDaughter].radiusB
                    
                    # set values          
                    
                    if vessels[rightDaughter].angleToMother is not None:
                        drot = -moRotY+(1+2/180.0*moRotY)*( (1+2/180.0*moRotY)*vessels[rightDaughter].angleToMother + moRot)
                    else:
                        drot =  -moRotY+(1+2/180.0*moRotY)*((1+2/180.0*moRotY)*30.0 + moRot) 
                                   
                    # set z position
                    pos_z = moLength*cos(moRot*pi/180)*(1+2/180.0*moRotY) # truns with 180 degree
                    #set rotation
                    if moRot != 0:
                        pos_y = -moLength*sin(moRot*pi/180)   
                    else:
                        pos_y = 0.0
                    dpos = moPos[0],moPos[1]+pos_y,moPos[2]+pos_z
                    
                    ## create visualisation Data
                    # create the data points 
                    points =  self.generate_vertices(r_init=[dRadiusA,dRadiusB],length = dLength,N = self.vesselN[rightDaughter])
                    # create scalar color map
                    col = linspace(0.0, 0.0, self.vesselN[rightDaughter])
                    color  = self.generate_colorScalars(col, N = self.vesselN[rightDaughter])
                    # generate structured grid, data source
                    self.vizVTKGrids[rightDaughter] = tvtk.StructuredGrid(dimensions=(self.dims[1], self.dims[0], self.vesselN[rightDaughter]))
                    #set data points
                    self.vizVTKGrids[rightDaughter].points = points
                    #set scalar color map
                    self.vizVTKGrids[rightDaughter].point_data.scalars = ravel(color.copy())
                    self.vizVTKGrids[rightDaughter].point_data.scalars.name =  'scalars' #self.scalarNames[self.scalarPosition]
                    
                    # create datasource for mayavi
                    self.vizSources[rightDaughter] = VTKDataSource(data = self.vizVTKGrids[rightDaughter])
                    
                    # create surface actor for the datas ource
                    self.vizActors[rightDaughter] = Surface()
                    self.vizActors[rightDaughter].actor.actor.position = dpos
                    self.vizActors[rightDaughter].actor.actor.rotate_x(drot)
                    
                    
                    self.velArrowGrid[rightDaughter] = tvtk.StructuredGrid(dimensions=(11, 1, self.vesselN[rightDaughter]))
                    #set data points
                    self.velArrowGrid[rightDaughter].points = self.generate_velocity_points(self.vesselLengths[rightDaughter], self.vesselN[rightDaughter], r_init=[dRadiusA,dRadiusB])
                    #set scalar color map
                    color = self.generate_velocity_colorScalars(self.vesselN[rightDaughter])
                    self.velArrowGrid[rightDaughter].point_data.vectors = color[0].copy()
                    self.velArrowGrid[rightDaughter].point_data.vectors.name = 'vectors'
                    self.velArrowGrid[rightDaughter].point_data.scalars = color[1].copy()
                    self.velArrowGrid[rightDaughter].point_data.scalars.name = 'scalars'
                    # create source
                    self.velArrowSrc[rightDaughter] = VTKDataSource(data = self.velArrowGrid[rightDaughter])
                    # create vectors data
                    self.glyph[rightDaughter] = Glyph()
                    self.glyph[rightDaughter].glyph.scale_mode = 'scale_by_scalar'
                    self.glyph[rightDaughter].glyph.color_mode = 'color_by_scalar'
                    self.glyph[rightDaughter].glyph.glyph_source.glyph_source = self.glyph[rightDaughter].glyph.glyph_source.glyph_dict['arrow_source']
                    self.glyph[rightDaughter].glyph.glyph_source.glyph_position = 'tail'
                    self.glyph[rightDaughter].glyph.glyph.scale_factor = self.vesselLengths[rightDaughter] / self.vesselN[rightDaughter] * 4 / 7
                    # set position
                    self.glyph[rightDaughter].actor.actor.position = dpos
                    self.glyph[rightDaughter].actor.actor.rotate_x(drot)
                
                    
                    # check for children
                    if vessels[rightDaughter].leftDaughter is not None:
                        viz.append(rightDaughter)     
            
    
