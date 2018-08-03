from Node import *
from Element import *

eqn = 0
n1 = Node(0, 0, 0)
n2 = Node(1, 1, 0)
e1 = Element(1, 1, n1, n2)
c1 = Constraint(False, True, True)
print(c1.constraint_print())
n1.set_constraint(c1)
eqn = n1.enumerate_dofs(eqn)
eqn = n2.enumerate_dofs(eqn)
print(e1.get_e1())
print(e1.compute_stiffness_matrix())
print(eqn)
print(n1.get_dof_numbers())
print(n2.get_dof_numbers())