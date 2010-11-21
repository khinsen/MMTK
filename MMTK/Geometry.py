# This module defines some geometrical objects in 3D-space.
#
# Written by Konrad Hinsen
#

"""
Elementary geometrical objects

There are essentially two kinds of geometrical objects: shape objects
(spheres, planes, etc.), from which intersections can be calculated,
and lattice objects, which define a regular arrangements of points.
"""

__docformat__ = 'restructuredtext'

from Scientific.Geometry import Vector
from Scientific import N


# Error type
class GeomError(Exception):
    pass

# Small number
eps = 1.e-16

#
# The base class
#
class GeometricalObject3D(object):

    """
    3D shape object

    This is an abstract base class. To create 3D objects,
    use one of its subclasses.
    """
    
    def intersectWith(self, other):
        """
        :param other: another 3D object
        :returns: a 3D object that represents the intersection with other
        """
	if self.__class__ > other.__class__:
	    self, other = other, self
	try:
	    f, switch = _intersectTable[(self.__class__, other.__class__)]
	    if switch:
		return f(other, self)
	    else:
		return f(self, other)
	except KeyError:
	    raise GeomError("Can't calculate intersection of " +
			     self.__class__.__name__ + " with " +
			     other.__class__.__name__)

    def volume(self):
        """
        :returns: the volume of the object
        :rtype: float
        """
        raise NotImplementedError

    def hasPoint(self, point):
        """
        :param point: a point in 3D space
        :type point: Scientific.Geometry.Vector
        :returns: True of the point lies on the surface of the object
        :rtype: bool
        """
	return self.distanceFrom(point) < eps

    # subclasses that enclose a volume should override this method
    # a return value of None indicates "don't know", "can't compute",
    # or "not implemented (yet)".
    def enclosesPoint(self, point):
        """
        :param point: a point in 3D space
        :type point: Scientific.Geometry.Vector
        :returns: True of the point is inside the volume of the object
        :rtype: bool
        """
        return None

_intersectTable = {}

#
# Boxes
#
class Box(GeometricalObject3D):

    """
    Rectangular box aligned with the coordinate axes
    """

    def __init__(self, corner1, corner2):
        """
        :param corner1: one corner of the box
        :type corner1: Scientific.Geometry.Vector
        :param corner2: the diagonally opposite corner
        :type corner2: Scientific.Geometry.Vector
        """
        c1 = N.minimum(corner1.array, corner2.array)
        c2 = N.maximum(corner1.array, corner2.array)
        self.corners = c1, c2

    def __repr__(self):
        return 'Box(' + `Vector(self.corners[0])` + ', ' \
               + `Vector(self.corners[1])` + ')'

    __str__ = __repr__

    def volume(self):
        c1, c2 = self.corners
        return N.multiply.reduce(c2-c1)

    def hasPoint(self, point):
        c1, c2 = self.corners
        min1 = N.minimum.reduce(N.fabs(point.array-c1))
        min2 = N.minimum.reduce(N.fabs(point.array-c2))
        return min1 < eps or min2 < eps

    def enclosesPoint(self, point):
        c1, c2 = self.corners
        out1 = N.logical_or.reduce(N.less(point.array-c1, 0))
        out2 = N.logical_or.reduce(N.less_equal(c2-point.array, 0))
        return not (out1 or out2)

    def cornerPoints(self):
        (c1x, c1y, c1z), (c2x, c2y, c2z) = self.corners
        return [Vector(c1x, c1y, c1z),
                Vector(c1x, c1y, c2z),
                Vector(c1x, c2y, c1z),
                Vector(c2x, c1y, c1z),
                Vector(c2x, c2y, c1z),
                Vector(c2x, c1y, c2z),
                Vector(c1x, c2y, c2z),
                Vector(c2x, c2y, c2z)]

#
# Spheres
#
class Sphere(GeometricalObject3D):

    """
    Sphere
    """

    def __init__(self, center, radius):
        """
        :param center: the center of the sphere
        :type center: Scientific.Geometry.Vector
        :param radius: the radius of the sphere
        :type radius: float
        """
	self.center = center
	self.radius = radius

    def __repr__(self):
        return 'Sphere(' + `self.center` + ', ' + `self.radius` + ')'
    __str__ = __repr__

    def volume(self):
	return (4.*N.pi/3.) * self.radius**3

    def hasPoint(self, point):
        return N.fabs((point-self.center).length()-self.radius) < eps

    def enclosesPoint(self, point):
        return (point - self.center).length() < self.radius

#
# Cylinders
#
class Cylinder(GeometricalObject3D):

    """
    Cylinder
    """

    def __init__(self, center1, center2, radius):
        """
        :param center1: the center of the bottom circle
        :type center1: Scientific.Geometry.Vector
        :param center2: the center of the top circle
        :type center2: Scientific.Geometry.Vector
        :param radius: the radius of the cylinder
        :type radius: float
        """
        self.center1 = center1            # center of base
        self.center2 = center2            # center of top
        self.radius = radius
        self.height = (center2-center1).length()

    def volume(self):
        return N.pi*self.radius*self.radius*self.height

    def __repr__(self):
        return 'Cylinder(' + `self.center1` + ', ' + `self.center2` + \
               ', ' + `self.radius` + ')'
    __str__ = __repr__

    def hasPoint(self, point):
        center_line = LineSegment(self.center1, self.center2)
        pt = center_line.projectionOf(point)
        if pt is None:
            return 0
        return N.fabs((point - pt).length() - self.radius) < eps

    def enclosesPoint(self, point):
        center_line = LineSegment(self.center1, self.center2)
        pt = center_line.projectionOf(point)
        if pt is None:
            return 0
        return (point - pt).length() < self.radius

#
# Planes
#
class Plane(GeometricalObject3D):

    """
    2D plane in 3D space
    """

    def __init__(self, *args):
        """
        :param args: three points (of type Scientific.Geometry.Vector)
                     that are not collinear, or a point in the plane and
                     the normal vector of the plane
        """
	if len(args) == 2:   # point, normal
	    self.normal = args[1].normal()
	    self.distance_from_zero = self.normal*args[0]
	else:                # three points
	    v1 = args[1]-args[0]
	    v2 = args[2]-args[1]
	    self.normal = (v1.cross(v2)).normal()
	    self.distance_from_zero = self.normal*args[1]

    def __repr__(self):
        return 'Plane(' + str(self.normal*self.distance_from_zero) + \
               ', ' + str(self.normal) + ')'
    __str__ = __repr__

    def distanceFrom(self, point):
	return abs(self.normal*point-self.distance_from_zero)

    def projectionOf(self, point):
        return point - (self.normal*point-self.distance_from_zero)*self.normal

    def rotate(self, axis, angle):
	point = rotatePoint(self.distance_from_zero*self.normal, axis, angle)
	normal = rotateDirection(self.normal, axis, angle)
	return Plane(point, normal)

    def volume(self):
	return 0.

#
# Infinite cones
#
class Cone(GeometricalObject3D):

    """
    Cone
    """

    def __init__(self, center, axis, angle):
        """
        :param center: the center (tip) of the cone
        :type center: Scientific.Geometry.Vector
        :param axis: the direction of the axis of rotational symmetry
        :type axis: Scientific.Geometry.Vector
        :param angle: the angle between any straight line on the cone
                      surface and the axis of symmetry
        :type angle: float
        """
	self.center = center
	self.axis = axis.normal()
	self.angle = angle

    def __repr__(self):
        return 'Cone(' + `self.center` + ', ' + `self.axis` + ',' + \
               `self.angle` + ')'
    __str__ = __repr__

    def volume(self):
	return None

#
# Circles
#
class Circle(GeometricalObject3D):

    """
    2D circle in 3D space
    """

    def __init__(self, center, normal, radius):
        """
        :param center: the center of the circle
        :type center: Scientific.Geometry.Vector
        :param normal: the normal vector of the circle's plane
        :type normal: Scientific.Geometry.Vector
        :param radius: the radius of the circle
        :type radius: float
        """
	self.center = center
	self.normal = normal
	self.radius = radius

    def planeOf(self):
        return Plane(self.center, self.normal)

    def __repr__(self):
        return 'Circle(' + `self.center` + ', ' + `self.normal` + \
               ', ' + `self.radius` + ')'
    __str__ = __repr__

    def volume(self):
        return 0.

    def distanceFrom(self, point):
        plane = self.planeOf()
        project_on_plane = plane.projectionOf(point)
        center_to_projection = project_on_plane - self.center
        if center_to_projection.length() < eps:
            return 0
        closest_point = self.center + self.radius*center_to_projection.normal()
        return (point - closest_point).length()

#
# Lines
#
class Line(GeometricalObject3D):

    """
    Line
    """

    def __init__(self, point, direction):
        """
        :param point: any point on the line
        :type point: Scientific.Geometry.Vector
        :param direction: the direction of the line
        :type direction: Scientific.Geometry.Vector
        """
	self.point = point
	self.direction = direction.normal()

    def distanceFrom(self, point):
        """
        :param point: a point in space
        :type point: Scientific.Geometry.Vector
        :returns: the smallest distance of the point from the line
        :rtype: float
        """
	d = self.point-point
	d = d - (d*self.direction)*self.direction
	return d.length()

    def projectionOf(self, point):
        """
        :param point: a point in space
        :type point: Scientific.Geometry.Vector
        :returns: the orthogonal projection of the point onto the line
        :rtype: Scientific.Geometry.Vector
        """
	d = self.point-point
	d = d - (d*self.direction)*self.direction
	return point+d

    def perpendicularVector(self, plane):
        """
        :param plane: a plane
        :type plane: Plane
        :returns: a vector in the plane perpendicular to the line
        :rtype: Scientific.Geometry.Vector
        """
        return self.direction.cross(plane.normal)

    def __repr__(self):
        return 'Line(' + `self.point` + ', ' + `self.direction` + ')'
    __str__ = __repr__

    def volume(self):
	return 0.


class LineSegment(Line):

    def __init__(self, point1, point2):
        Line.__init__(self, point1, point2 - point1)
        self.point2 = point2

    def __repr__(self):
        return 'LineSegment(' + `self.point` + ', ' + `self.point2` + ')'
    __str__ = __repr__

    def distanceFrom(self, point):
        pt = self.projectionOf(point)
        if pt is not None:
            return (pt - point).length()
        d1 = (self.point - point).length()
        d2 = (self.point2 - point).length()
        return min(d1, d2)

    def projectionOf(self, point):
        d = self.point-point
        d = d - (d*self.direction)*self.direction
        pt = point+d
        if self.isWithin(pt):
            return pt
        return None

    def isWithin(point):
        v1 = point - self.point
        v2 = point - self.point2
        if abs(v1 * v2) < eps:
            return 0
        return not Same_Dir(v1, v2)

#
# Intersection calculations
#
def _addIntersectFunction(f, class1, class2):
    switch = class1 > class2
    if switch:
	class1, class2 = class2, class1
    _intersectTable[(class1, class2)] = (f, switch)


# Box with box

def _intersectBoxBox(box1, box2):
    c1 = N.maximum(box1.corners[0], box2.corners[0])
    c2 = N.minimum(box1.corners[1], box2.corners[1])
    if N.logical_or.reduce(N.greater_equal(c1, c2)):
        return None
    return Box(Vector(c1), Vector(c2))

_addIntersectFunction(_intersectBoxBox, Box, Box)

# Sphere with sphere

def _intersectSphereSphere(sphere1, sphere2):
    r1r2 = sphere2.center-sphere1.center
    d = r1r2.length()
    if d > sphere1.radius+sphere2.radius:
	return None
    if d+min(sphere1.radius, sphere2.radius) < \
                             max(sphere1.radius, sphere2.radius):
	return None
    x = 0.5*(d**2 + sphere1.radius**2 - sphere2.radius**2)/d
    h = N.sqrt(sphere1.radius**2-x**2)
    normal = r1r2.normal()
    return Circle(sphere1.center + x*normal, normal, h)

_addIntersectFunction(_intersectSphereSphere, Sphere, Sphere)

# Sphere with cone

def _intersectSphereCone(sphere, cone):
    if sphere.center != cone.center:
	raise GeomError("Not yet implemented")
    from_center = sphere.radius*N.cos(cone.angle)
    radius = sphere.radius*N.sin(cone.angle)
    return Circle(cone.center+from_center*cone.axis, cone.axis, radius)

_addIntersectFunction(_intersectSphereCone, Sphere, Cone)

# Plane with plane

def _intersectPlanePlane(plane1, plane2):
    if abs(abs(plane1.normal*plane2.normal)-1.) < eps:
	if abs(plane1.distance_from_zero-plane2.distance_from_zero) < eps:
	    return plane1
	else:
            return None
    else:
	direction = plane1.normal.cross(plane2.normal)
	point_in_1 = plane1.distance_from_zero*plane1.normal
	point_in_both = point_in_1 - (point_in_1*plane2.normal -
				      plane2.distance_from_zero)*plane2.normal
	return Line(point_in_both, direction)

_addIntersectFunction(_intersectPlanePlane, Plane, Plane)

# Circle with plane

def _intersectCirclePlane(circle, plane):
    if abs(abs(circle.normal*plane.normal)-1.) < eps:
	if plane.hasPoint(circle.center):
	    return circle
	else:
	    return None
    else:
	line = plane.intersectWith(Plane(circle.center, circle.normal))
	x = line.distanceFrom(circle.center)
	if x > circle.radius:
	    return None
	else:
	    angle = N.arccos(x/circle.radius)
	    along_line = N.sin(angle)*circle.radius
	    normal = circle.normal.cross(line.direction)
	    if line.distanceFrom(circle.center+normal) > x:
		normal = -normal
	    return (circle.center+x*normal-along_line*line.direction,
		    circle.center+x*normal+along_line*line.direction)
	    
_addIntersectFunction(_intersectCirclePlane, Circle, Plane)

#
# Rotation
#
def rotateDirection(vector, axis, angle):
    s = N.sin(angle)
    c = N.cos(angle)
    c1 = 1-c
    try:
        axis = axis.direction
    except AttributeError:
        pass
    return s*axis.cross(vector) + c1*(axis*vector)*axis + c*vector

def rotatePoint(point, axis, angle):
    return axis.point + rotateDirection(point-axis.point, axis, angle)

#
# Lattices
#

#
# Lattice base class
#
class Lattice(object):

    """
    General lattice

    Lattices are special sequence objects that contain vectors
    (points on the lattice) or objects that are constructed as
    functions of these vectors. Lattice objects behave like
    lists, i.e. they permit indexing, length inquiry, and iteration
    by 'for'-loops. Note that the lattices represented by these
    objects are finite, they have a finite (and fixed) number
    of repetitions along each lattice vector.

    This is an abstract base class. To create lattice objects,
    use one of its subclasses.
    """

    def __init__(self, function):
	if function is not None:
	    self.elements = map(function, self.elements)

    def __getitem__(self, item):
	return self.elements[item]

    def __setitem__(self, item, value):
	self.elements[item] = value

    def __len__(self):
	return len(self.elements)

#
# General rhombic lattice
#
class RhombicLattice(Lattice):

    """
    Rhombic lattice
    """

    def __init__(self, elementary_cell, lattice_vectors, cells,
                 function=None, base=None):
        """
        :param elementary_cell: a list of points in the elementary cell
        :param lattice_vectors: a list of lattice vectors. Each lattice
                                vector defines a lattice dimension (only
                                values from one to three make sense) and
                                indicates the displacement along this
                                dimension from one cell to the next.
        :param cells: a list of integers, whose length must equal the number
                      of dimensions. Each entry specifies how often a cell is
                      repeated along this dimension.
        :param function: a function that is called for every lattice point with
                         the vector describing the point as argument. The return
                         value of this function is stored in the lattice object.
                         If the function is 'None', the vector is directly
                         stored in the lattice object.
        """
	if len(lattice_vectors) != len(cells):
	    raise TypeError('Inconsistent dimension specification')
        if base is None:
            base = Vector(0, 0, 0)
	self.dimension = len(lattice_vectors)
	self.elements = []
	self.makeLattice(elementary_cell, lattice_vectors, cells, base)
	Lattice.__init__(self, function)

    def makeLattice(self, elementary_cell, lattice_vectors, cells, base):
	if len(cells) == 0:
	    for p in elementary_cell:
		self.elements.append(p+base)
	else:
	    for i in range(cells[0]):
		self.makeLattice(elementary_cell, lattice_vectors[1:],
				 cells[1:], base+i*lattice_vectors[0])
	    
#
# Bravais lattice
#
class BravaisLattice(RhombicLattice):

    """
    Bravais lattice

    A Bravais lattice is a special case of a general rhombic lattice
    in which the elementary cell contains only one point.
    """

    def __init__(self, lattice_vectors, cells, function=None, base=None):
        """
        :param lattice_vectors: a list of lattice vectors. Each lattice
                                vector defines a lattice dimension (only
                                values from one to three make sense) and
                                indicates the displacement along this
                                dimension from one cell to the next.
        :param cells: a list of integers, whose length must equal the number
                      of dimensions. Each entry specifies how often a cell is
                      repeated along this dimension.
        :param function: a function that is called for every lattice point with
                         the vector describing the point as argument. The return
                         value of this function is stored in the lattice object.
                         If the function is 'None', the vector is directly
                         stored in the lattice object.
        """
	cell = [Vector(0,0,0)]
	RhombicLattice.__init__(self, cell, lattice_vectors, cells,
                                function, base)

#
# Simple cubic lattice
#
class SCLattice(BravaisLattice):

    """
    Simple Cubic lattice

    A Simple Cubic lattice is a special case of a Bravais lattice
    in which the elementary cell is a cube.
    """

    def __init__(self, cellsize, cells, function=None, base=None):
        """
        :param cellsize: the edge length of the elementary cell
        :type cellsize: float
        :param cells: a list of integers, whose length must equal the number
                      of dimensions. Each entry specifies how often a cell is
                      repeated along this dimension.
        :param function: a function that is called for every lattice point with
                         the vector describing the point as argument. The return
                         value of this function is stored in the lattice object.
                         If the function is 'None', the vector is directly
                         stored in the lattice object.
        """
	lattice_vectors = (cellsize*Vector(1., 0., 0.),
                           cellsize*Vector(0., 1., 0.),
                           cellsize*Vector(0., 0., 1.))
	if type(cells) != type(()):
	    cells = 3*(cells,)
	BravaisLattice.__init__(self, lattice_vectors, cells, function, base)

#
# Body-centered cubic lattice
#
class BCCLattice(RhombicLattice):

    """
    Body-Centered Cubic lattice

    A Body-Centered Cubic lattice has two points per elementary cell.
    """

    def __init__(self, cellsize, cells, function=None, base=None):
        """
        :param cellsize: the edge length of the elementary cell
        :type cellsize: float
        :param cells: a list of integers, whose length must equal the number
                      of dimensions. Each entry specifies how often a cell is
                      repeated along this dimension.
        :param function: a function that is called for every lattice point with
                         the vector describing the point as argument. The return
                         value of this function is stored in the lattice object.
                         If the function is 'None', the vector is directly
                         stored in the lattice object.
        """
        cell = [Vector(0,0,0), cellsize*Vector(0.5,0.5,0.5)]
	lattice_vectors = (cellsize*Vector(1., 0., 0.),
                           cellsize*Vector(0., 1., 0.),
                           cellsize*Vector(0., 0., 1.))
	if type(cells) != type(()):
	    cells = 3*(cells,)
	RhombicLattice.__init__(self, cell, lattice_vectors, cells,
                                function, base)

#
# Face-centered cubic lattice
#
class FCCLattice(RhombicLattice):

    """Face-Centered Cubic lattice

    A Face-Centered Cubic lattice has four points per elementary cell.
    """

    def __init__(self, cellsize, cells, function=None, base=None):
        """
        :param cellsize: the edge length of the elementary cell
        :type cellsize: float
        :param cells: a list of integers, whose length must equal the number
                      of dimensions. Each entry specifies how often a cell is
                      repeated along this dimension.
        :param function: a function that is called for every lattice point with
                         the vector describing the point as argument. The return
                         value of this function is stored in the lattice object.
                         If the function is 'None', the vector is directly
                         stored in the lattice object.
        """
        cell = [Vector(0,0,0),
                cellsize*Vector(  0,0.5,0.5),
                cellsize*Vector(0.5,  0,0.5),
                cellsize*Vector(0.5,0.5,  0)]
	lattice_vectors = (cellsize*Vector(1., 0., 0.),
                           cellsize*Vector(0., 1., 0.),
                           cellsize*Vector(0., 0., 1.))
	if type(cells) != type(()):
	    cells = 3*(cells,)
	RhombicLattice.__init__(self, cell, lattice_vectors, cells,
                                function, base)
