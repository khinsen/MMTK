# This example demonstrates using MMTK and the Scientific
# Python VRML2 (VRML97) visualisation module for drawing a
# protein from several different angles using the camera object.
#
# The camera object is not available in the VRML, VPython
# and VMD visualisation modules.
#
# The new space filling van der Waals representation (model='vdw')
# is also demonstrated.  [Introduced in June 2004, should be
# available in the next release after Scientific Python 2.4.6]
#
# Written by Peter Cock
# p.j.a.cock [at] warwick.ac.uk
# Molecular Organisation and Assembly in Cells (MOAC) Doctoral
# Training Centre, University of Warwick, CV4 7AL, UK
#
# last revision: 2004-06-09
#
from MMTK import *
from MMTK.Proteins import Protein

# Import the graphics module. You can substitute other graphics
# module names to make the example use that module, but not
# everything supports the camera object, used for view points
from Scientific.Visualization import VRML2; visualization_module = VRML2
#In fact VRML, VPython and VMD do not support cameras:-
#from Scientific.Visualization import VPython; visualization_module = VPython
#from Scientific.Visualization import VRML; visualization_module = VRML
#from Scientific.Visualization import VMD; visualization_module = VMD

print "Loading protein file"
#This can be any name from the MMTK database, or a PDB file.
protein = Protein('insulin')
#Note that by default this looks for all the atoms.  Some PDB files
#don't include the hydrogens, in which case MMTK will guess them.
#Some PDB files don't even have the positions of the side chains,
#and MMTK will not be able to guess their locations.

print "Finding centre of mass and moment of inertia"
center, inertia = protein.centerAndMomentOfInertia()

print "Creating view points"
#Use the centre of mass as the centre of the picture, with the camera
#pulled back from it.
#
#X and Y are the axis of the screen
#Y is in (-ve) and out (+ve) of the screen
#
#I am rotating the camera about the y-axis (which is vertical in screen)
#in units of pi/2 in order to give four views (Front, Right, Back, Left)
distance_away = 8.0
front_cam = visualization_module.Camera(position=[center[0],center[1],center[2]+distance_away],
                    description="Front")
right_cam = visualization_module.Camera(position=[center[0]+distance_away,center[1],center[2]],
                    orientation=(Vector(0, 1, 0),3.14159*0.5),
                    description="Right")
back_cam = visualization_module.Camera(position=[center[0],center[1],center[2]-distance_away],
                    orientation=(Vector(0, 1, 0),3.14159),
                    description="Back")
left_cam = visualization_module.Camera(position=[center[0]-distance_away,center[1],center[2]],
                    orientation=(Vector(0, 1, 0),3.14159*1.5),
                    description="Left")

#Other visualisation methods ("models") you might like to try:-
model_name = 'vdw'
#model_name = 'tube'
#model_name = 'ball_and_stick'

print "Creating " + model_name + " graphic of Protein"
#You can also specify colour for the whole molecule if you want,
#otherwise the individual atom colours are used (Oxygen = Red,
#Hydrogen = White, Sulphur = Yellow etc)
graphics = protein.graphicsObjects(graphics_module = visualization_module,
                                   model = model_name)

print "Creating arrows below protein"
#In order to help visualise what the camera objects are doing, I
#am adding some arrows in the x-z planes above and below the molecule.
#
#The red arrows point towards the "front" camera, green towards the
#"right" camera, blue towards the "back" camera and yellow towards the
#"left" camera.
#
d = 2.0 #distance from centre of mass (above or below)
l = 2.0 #length of each arrow
graphics.append(visualization_module.Arrow(center + Vector(0,-d,0),
                center + Vector(0,-d,l), 0.1,
                material=visualization_module.DiffuseMaterial('red')))
graphics.append(visualization_module.Arrow(center + Vector(0,-d,0),
                center + Vector(l,-d,0), 0.1,
                material=visualization_module.DiffuseMaterial('green')))
graphics.append(visualization_module.Arrow(center + Vector(0,-d,0),
                center + Vector(0,-d,-l), 0.1,
                material=visualization_module.DiffuseMaterial('blue')))
graphics.append(visualization_module.Arrow(center + Vector(0,-d,0),
                center + Vector(-l,-d,0), 0.1, 
                material=visualization_module.DiffuseMaterial('yellow')))

print "Creating arrows above protein"
graphics.append(visualization_module.Arrow(center + Vector(0,+d,0),
                center + Vector(0,+d,l), 0.1,
                material=visualization_module.DiffuseMaterial('red')))
graphics.append(visualization_module.Arrow(center + Vector(0,+d,0),
                center + Vector(l,+d,0), 0.1,
                material=visualization_module.DiffuseMaterial('green')))
graphics.append(visualization_module.Arrow(center + Vector(0,+d,0),
                center + Vector(0,+d,-l), 0.1,
                material=visualization_module.DiffuseMaterial('blue')))
graphics.append(visualization_module.Arrow(center + Vector(0,+d,0),
                center + Vector(-l,+d,0), 0.1, 
                material=visualization_module.DiffuseMaterial('yellow')))

#print "Viewing front of " + model_name + " model..."
#visualization_module.Scene(graphics, cameras=[front_cam]).view()

#print "Viewing right of " + model_name + " model..."
#visualization_module.Scene(graphics, cameras=[right_cam]).view()

#print "Viewing back of " + model_name + " model..."
#visualization_module.Scene(graphics, cameras=[back_cam]).view()

#print "Viewing left of " + model_name + " model..."
#visualization_module.Scene(graphics, cameras=[left_cam]).view()

print "Viewing " + model_name + " model with all four cameras..."
visualization_module.Scene(graphics, cameras=[front_cam,right_cam,back_cam,left_cam]).view()
#Your VRML2 viewing program should let you switch between the
#four cameras (and show you their names "Front", "Right", etc).
#
#For example, on Windows, using GLView 4.4 this is done using
#the menu "Camera","Viewpoints" (or shortcut key F6)
#Available here: http://home.snafu.de/hg/

print "Done"
