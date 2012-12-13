print ""
print "Check for necessary moduls and 3thrd party libarys needed for vascular1DFlow and vnc"
print ""
print "check for matplotlib ..."
installed = "                         .. installed."
try:
    import matplotlib.pyplot 
    print installed
except:
    print "IMPORT ERROR; no version of matplotlib.pyplot found"

print "check for  numpy ..."
try:
    import numpy
    print installed
except:
    print "IMPORT ERROR; no version of numpy found"

print "check for  scipy ..."
try:
    import scipy
    print installed
except:
    print "IMPORT ERROR; no version of scipy found"

print "check for  mayavi2 ..."
try:
    import mayavi
    print installed
except:
    print "IMPORT ERROR; no version of mayavi2 found"

print "check for  tvtk ..."
try:
    import tvtk
    print installed
except:
    print "IMPORT ERROR; no version of tvtk found"

print "check for  lxml2 ..."
try:
    import lxml
    print installed
except:
    print "IMPORT ERROR; no version of lxml2 found"

print "check for  pydot ..."
try:
    import pydot
    print installed
except:
    print "IMPORT ERROR; no version of pydot found"

print "check for  gtk ..."
try:
    import gtk
    print installed
except:
    print "IMPORT ERROR; no version of gtk found" 
