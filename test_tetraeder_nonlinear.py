from Structure import *
from MyViewer import *
from math import pi, sqrt, pow
from numpy import linspace

t = linspace(0, 1, 20)
structure = Structure(t)

lb = 15.0
r = 457.2/2000.0
t = 10.0/1000.0
a = pi * (pow(r,2) - pow(r-t, 2))
e = 2.1e11
c1 = Constraint(False, False, False)
c2 = Constraint(True, True, False)
f = Force(0.0, -20.0e3, -100.0e3)
n1 = structure.add_node(0.0, 0.0, lb * sqrt(2.0 / 3.0))
n2 = structure.add_node(0.0, lb / sqrt(3.0), 0)
n3 = structure.add_node(-lb / 2.0, -lb / sqrt(12.0), 0.0)
n4 = structure.add_node(lb / 2.0, -lb / sqrt(12.0), 0.0)
n1.set_force(f)
n2.set_constraint(c1)
n3.set_constraint(c1)
n4.set_constraint(c2)
structure.add_element(CrisfieldTruss(e, a, n1, n2))
structure.add_element(CrisfieldTruss(e, a, n1, n3))
structure.add_element(CrisfieldTruss(e, a, n1, n4))
structure.add_element(CrisfieldTruss(e, a, n2, n3))
structure.add_element(CrisfieldTruss(e, a, n3, n4))
structure.add_element(CrisfieldTruss(e, a, n4, n2, -1e12))
u = structure.nonlinear_solve()
structure.print_structure()
structure.print_results()
v = MyViewer(structure)
v.draw_elements()
v.set_force_scale()
v.draw_element_forces()
v.draw_constraints()
v.draw_nodal_forces()
v.set_displacement_scale()
v.draw_displacements()
v.set_visible()
print(structure.u)