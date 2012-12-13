import matplotlib.pyplot as plt   
import numpy as np
from optparse import OptionParser
import cPickle

import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(cur+'/../UtilityLib')
from modulXML import loadNetworkFromXML 

sys.path.append(cur+'/../NetworkLib')
from classVascularNetwork import VascularNetwork 

import gtk,sys

class PyApp(gtk.Window):
    def __init__(self,nNodes,filenames,vesselName,updateVesselId):
        super(PyApp,self).__init__()
        self.connect('destroy',gtk.main_quit)
        self.set_size_request(680,670)
        self.set_position(gtk.WIN_POS_CENTER)
        self.set_title('2d Visualisation - '+vesselName)
        
        #gridpoints
        self.N = nNodes
        self.filenames = filenames
        self.pixBuf = []
        
        scale = gtk.HScale()
        scale.set_range(0,self.N-1)
        scale.set_increments(1,1)
        scale.set_digits(0)
        scale.set_size_request(400,30)
        scale.connect('value-changed',self.on_changed)
        
        self.load_pixbufs()
        self.val = 0
        self.image = gtk.Image()
        self.image.set_from_pixbuf(self.pixBuf[self.val])

        vbox = gtk.VBox(False,1)
        
        valignIm = gtk.Alignment(0,1 , 1, 0)
        valignIm.add(self.image)
        
        valign = gtk.Alignment(0, 1, 1, 0)
        valign.add(scale)
        
        vbox.pack_start(valignIm)
        vbox.pack_start(valign)
        
        self.add(vbox)         
        
        if updateVesselId == True:
            gtk.timeout_add(3,update_vesselId)
        
        self.show_all()
        
    def load_pixbufs(self):
        try: 
            for filename in self.filenames:
                a = gtk.gdk
                pixbuf = a.pixbuf_new_from_file(filename)
                pixbuf = pixbuf.scale_simple(680, 640, a.INTERP_BILINEAR)

                self.pixBuf.append(pixbuf)
            
        except Exception, e:
            print "reading Error"
            print e.massage
            sys.exit(1)
    
    def on_changed(self,widget):
        self.val = int(widget.get_value())
        self.image.set_from_pixbuf(self.pixBuf[self.val])
        
    def on_image_resize(self, widget, event, window):
        allocation = widget.get_allocation()
        if self.temp_height != allocation.height or self.temp_width != allocation.width:
            self.temp_height = allocation.height
            self.temp_width = allocation.width
            a = gtk.gdk
            self.pixBuf[self.val] = self.pixBuf[self.val].scale_simple(allocation.width, allocation.height, a.INTERP_BILINEAR)
            widget.set_from_pixbuf(self.pixBuf[self.val])

    def update_vesselId(self):
        f = open('tempVesselId.txt',r)
        self.vesselId = int(f.read())

def main2Dplot():
        
    parser = OptionParser()
    parser.add_option("-f", "--file", dest='networkName',
                      help="open file with networkName", metavar="FILE")
    parser.add_option("-i", "--vesselId", dest="vesselId", 
                      help="id of the vessel to visualize")
    parser.add_option("-n", "--dataNumber", dest='dataNumber',
                      help="number of the solution data (last number in filename), default = 1, max 999", metavar = "FILE")
    (options, args) = parser.parse_args()
    
    networkName = 'oldNetwork2'
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
    
    vesselId = 0
    updateVesselId = False
    if options.vesselId != None:
        vesselId = int(options.vesselId)
    else:
        updateVesselId = True
    
    vascularNetwork = None
    vascularNetwork = loadNetworkFromXML(filename=filename)
    if vascularNetwork == None: exit()
    vascularNetwork.evaluateConnections()
    vascularNetwork.initialize()
    
    if vesselId not in vascularNetwork.vessels.keys():
        "No vessel with given id in current vascular network"
        sys.exit(1)
    vesselName = vascularNetwork.vessels[vesselId].name

    solutionData = []
    endName1 = "_SolutionData_"
    endName2 = ".pickle"
    networkPath = str(cur+'/../'+"/NetworkFiles/")
    filePath= ''.join([networkPath,networkName,'/',networkName,endName1,dataNumber,endName2])
    try:
        FILE = open(filePath,"rb")
        solutionData = cPickle.load(FILE)
    except:
        print "no solution file with given data number"
        exit()
        
    pRef = vascularNetwork.globalFluid['pref']
    totalTime = vascularNetwork.simulationContext['totalTime']
    
    save = True
    
    saveDirectory = ''.join([networkPath,networkName,'/2dTimePlots/'])
    
    for dir in [saveDirectory]: 
        if not os.path.exists(dir): os.makedirs(dir) 

    nodeID = 0

    nNodes = int(vascularNetwork.vessels[vesselId].N)
    
    filenames = []
    
    Pmax = np.max(solutionData[0]['Pressure'][vesselId])/133.32
    Pmin = np.min(solutionData[0]['Pressure'][vesselId])/133.32
    Qmax = np.max(solutionData[0]['Flow'][vesselId])*1e6*60
    Qmin = np.min(solutionData[0]['Flow'][vesselId])*1e6*60
    
    for n in range(0,nNodes):
    
        Psol = solutionData[0]['Pressure'][vesselId][:,[n]]
        Qsol = solutionData[0]['Flow'][vesselId][:,[n]]
    
        timeNormal = np.linspace(0,totalTime,len(Psol))
        
        fig = plt.figure(1, edgecolor='k',dpi=80)
        fig.subplots_adjust(hspace=0.5)   
        fig.set_figwidth(8.27)
        fig.set_figheight((11.69/3)*2)
        ax = plt.subplot(2,1,1)    
        ax.plot(timeNormal,Psol/133.32,color='b' ,linestyle = '-',label='total', linewidth = 1.5) 
        ax.set_ylabel('Pressure [mmHg]')
        ax.set_xlabel('Time [t]')
        ax.set_xlim(0,totalTime)
        ax.set_ylim(Pmin,Pmax)
        
        ax2 = plt.subplot(2,1,2)    
        ax2.plot(timeNormal,Qsol*1e6*60,color='r' ,linestyle = '-',label='total', linewidth = 1.5) 
        ax2.set_ylabel('Flow [ml/min]')
        ax2.set_xlabel('Time [t]')
        ax2.set_xlim(0,totalTime)
        ax2.set_ylim(Qmin,Qmax)
        
        savePath = ''.join([saveDirectory,str(vesselId).zfill(2),str(n).zfill(2),'.png'])
        filenames.append(savePath)
        plt.savefig(savePath)
        plt.clf()

    PyApp(nNodes,filenames,vesselName,updateVesselId)
    
    gtk.main()

if __name__ == '__main__':
    main2Dplot()