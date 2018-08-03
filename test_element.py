from Structure import *
import math

n1 = Node(0, 0, 0)
n2 = Node(1, 1, 1)
e = Element(3, math.sqrt(3), n1, n2)
ke = e.compute_stiffness_matrix()
print(ke)
print(ke.shape)