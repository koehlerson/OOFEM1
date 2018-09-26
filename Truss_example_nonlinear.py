from Structure import *
from MyViewer import *
from numpy import linspace

t = linspace(0, 1, 20)
structure = Structure(t)
r = 0.2e-1
a = 0.01
e = 2.1e11

c1 = Constraint(False, False, False)
c2 = Constraint(False, False, True)
c3 = Constraint(True, False, True)
c4 = Constraint(True, False, False)
c5 = Constraint(False, True, False)
f1 = Force(1e4, -20.0e3, -100.0e2)
f2 = Force(0, 0, -100e1)
length = 10.
height = 10
width = 10.

n0 = structure.add_node(0., 0., height)
n1a = structure.add_node(length*(3/5.), -width, height*(3/5.))
n2a = structure.add_node(length, 0., height*(3/5.))
n3a = structure.add_node(length*(3/5.), width, height*(3/5.))
n1b = structure.add_node(-length*(3/5.), -width, height*(3/5.))
n2b = structure.add_node(-length, 0., height*(3/5.))
n3b = structure.add_node(-length*(3/5.), width, height*(3/5.))
n1 = structure.add_node(0, width*1.5, 0.)
n2 = structure.add_node(0, -width*1.5, 0.)
n3 = structure.add_node(length*1.2, width*0.7, 0.)
n4 = structure.add_node(-length*1.2, -width*0.7, 0.)
n5 = structure.add_node(-length*1.2, width*0.7, 0.)
n6 = structure.add_node(length*1.2, -width*0.7, 0.)

n0.set_force(f2)
n1a.set_force(f2)
n2a.set_force(f2)
n3a.set_force(f1)
n1b.set_force(f1)
n2b.set_force(f2)
#n3b.set_force(f1)
n1.set_constraint(c1)
n2.set_constraint(c1)
n3.set_constraint(c1)
n4.set_constraint(c1)
n5.set_constraint(c1)
n6.set_constraint(c1)

structure.add_element(CrisfieldTruss(e,a,n1a,n0))
structure.add_element(CrisfieldTruss(e,a,n2a,n0))
structure.add_element(CrisfieldTruss(e,a,n3a,n0))
structure.add_element(CrisfieldTruss(e,a,n1b,n0))
structure.add_element(CrisfieldTruss(e,a,n2b,n0))
structure.add_element(CrisfieldTruss(e,a,n3b,n0))
structure.add_element(CrisfieldTruss(e,a,n1a,n2a))
structure.add_element(CrisfieldTruss(e,a,n2a,n3a))
structure.add_element(CrisfieldTruss(e,a,n3a,n3b))
structure.add_element(CrisfieldTruss(e,a,n1b,n2b))
structure.add_element(CrisfieldTruss(e,a,n2b,n3b))
structure.add_element(CrisfieldTruss(e,a,n1b,n1a))
structure.add_element(CrisfieldTruss(e,a,n1,n3a))
structure.add_element(CrisfieldTruss(e,a,n1,n3b, prestress = 1e7))
structure.add_element(CrisfieldTruss(e,a,n2,n1a))
structure.add_element(CrisfieldTruss(e,a,n2,n1b))
structure.add_element(CrisfieldTruss(e,a,n3,n2a))
structure.add_element(CrisfieldTruss(e,a,n3,n3a))
structure.add_element(CrisfieldTruss(e,a,n4,n2b))
structure.add_element(CrisfieldTruss(e,a,n4,n1b))
structure.add_element(CrisfieldTruss(e,a,n5,n2b))
structure.add_element(CrisfieldTruss(e,a,n5,n3b))
structure.add_element(CrisfieldTruss(e,a,n6,n1a))
structure.add_element(CrisfieldTruss(e,a,n6,n2a))

u = structure.nonlinear_solve()
structure.print_structure()
structure.print_results()
v = MyViewer(structure)
v.set_displacement_scale()
v.set_force_scale()
v.draw_elements()
v.draw_constraints()
v.draw_nodal_forces()
v.draw_element_forces()
v.draw_displacements()
v.set_visible()