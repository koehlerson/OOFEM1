from pview3d import *
from math import sqrt, pi
from numpy.linalg import norm
from numpy import cross, array, count_nonzero, log10, floor, abs, sum


class MyViewer:

    def __init__(self, structure):
        self.Structure = structure
        self.__symbolScale = 1
        self.__displacementScale = 1
        self.__Viewer = Viewer()
        self.__forceScale = 1

    def draw_elements(self):
        """draw elements as a 3d cylinder"""
        cs = CylinderSet()
        cs.color = 'black'
        for i in self.Structure.Elements:
            a = i.get_area()
            r = sqrt(a/pi)
            cs.add_cylinder(i.get_node1().get_position(), i.get_node2().get_position(), r)
        self.__Viewer.add_object(cs)

    def draw_constraints(self, symbolscale = 1):
        """draw the constraints in an engineering manner"""
        self.__symbolScale = symbolscale
        for i in self.Structure.Nodes:
            constraint = i.get_constraint()
            for j in range(3):
                if constraint.is_free(j) != True:
                    cone = Cone()
                    if j == 0:
                        cone.direction = self.__symbolScale * array([1, 0, 0])
                        cone.center = i.get_position() - array([cone.height, 0, 0])
                    elif j == 1:
                        cone.direction = self.__symbolScale * array([0, 1, 0])
                        cone.center = i.get_position() - array([0, cone.height, 0])
                    elif j == 2:
                        cone.direction = self.__symbolScale * array([0, 0, 1])
                        cone.center = i.get_position() - array([0, 0, cone.height])
                    cone.height = self.__symbolScale*1
                    cone.radius = self.__symbolScale*0.5
                    cone.resolution = self.__symbolScale*40
                    cone.color = 'black'
                    self.__Viewer.add_object(cone)

    def draw_element_forces(self):
        """draw the element forces as rectangulars attached to the elements. compression red, tension blue"""
        ps = PolygonSet()
        for i in self.Structure.Elements:
            force = i.compute_force()
            if force < 0:
                color = -1
            elif force == 0:
                color = 1
            else:
                color = 2
            x1 = i.get_node1().get_position()
            x2 = i.get_node2().get_position()
            d = (x1 - x2) / norm(x1 - x2, 2)
            n1 = cross(d, array([0, 0, 1]))
            p = cross(n1, d)
            s1 = x1 + force * p * self.__forceScale
            s2 = x2 + force * p * self.__forceScale
            ps.insert_vertex(x1[0], x1[1], x1[2], color)
            ps.insert_vertex(s1[0], s1[1], s1[2], color)
            ps.insert_vertex(s2[0], s2[1], s2[2], color)
            ps.insert_vertex(x2[0], x2[1], x2[2], color)
            ps.polygon_complete()
            self.__Viewer.add_object(ps)
        ps.smooth = True
        ps.banded = True
        ps.color_by_data = True
        ps.outlines_visible = False
        ps.contour_lines_visible = False
        ps.create_colors()

    def draw_nodal_forces(self):
        """draws all nodal forces as arrows. compression/tension forces red/blue"""
        cs = CylinderSet()
        for i in self.Structure.Nodes:
            force = i.get_force()
            force = array([force.get_component(0), force.get_component(1), force.get_component(2)])
            if count_nonzero(force) != 0:
                force_norm = force / norm(force, 2)
                cone = Cone()
                cone.direction = force/norm(force, 2)
                cone.color = 'blue'
                cone.height = self.__symbolScale * 1
                cone.center = i.get_position() - cone.height*cone.direction
                cs.add_cylinder(i.get_position() - force_norm * cone.height, i.get_position() - (force_norm * 5)
                                * self.__symbolScale, self.__symbolScale * 0.1)
                self.__Viewer.add_object(cone)
                if sum(force) > 0:
                    cs.color = 'blue'
                    cone.color = 'blue'
                else:
                    cs.color='red'
                    cone.color='red'
            self.__Viewer.add_object(cone)
        self.__Viewer.add_object(cs)

    def draw_displacements(self):
        """draws the current configuration in 'antiquewhite' """
        cs = CylinderSet()
        cs.color = 'antiquewhite'
        for i in self.Structure.Elements:
            a = i.get_area()
            r = sqrt(a / pi)
            u1 = i.get_node1().get_displacement()
            u2 = i.get_node2().get_displacement()
            pos1 = i.get_node1().get_position() + (u1 * self.__displacementScale)
            pos2 = i.get_node2().get_position() + (u2 * self.__displacementScale)
            cs.add_cylinder(pos1, pos2, r)
        self.__Viewer.add_object(cs)

    def set_visible(self):
        """let the viewer run"""
        self.__Viewer.run()

    def set_displacement_scale(self, scale = 0):
        """change the displacement scale"""
        if scale == 0:
            scale = floor(log10(abs(self.Structure.u[-1])))
            scale = abs(scale)
            scale = 10 ** (scale.max() - 1.6)
            self.__displacementScale = scale
        else:
            self.__displacementScale = scale

    def set_force_scale(self, scale = 0):
        """change the force scale"""
        forces = []
        for i in self.Structure.Elements:
            if i.compute_force() == 0:
                forces.append(1)
            else:
                forces.append(i.compute_force())
        if scale == 0:
            scale = floor(log10(abs(forces)))
            scale = abs(scale)
            scale = 10 ** (scale.max()*-1)
            self.__forceScale = scale
        else:
            self.__forceScale = scale
