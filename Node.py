from Force import *
from Constraint import *
from numpy import zeros, array


class Node:

    def __init__(self, x1, x2, x3):
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.__dofNumbers = []
        self.__Force = Force()
        self.__Constraint = Constraint(True, True, True)
        self.__displacement = zeros(3)
        self.__position = array([self.x1, self.x2, self.x3])

    def set_constraint(self, Constraint):
        self.__Constraint = Constraint

    def get_constraint(self):
        return self.__Constraint

    def set_force(self, force):
        self.__Force = force

    def get_force(self):
        return self.__Force

    def enumerate_dofs(self, start):
        for i in range(3):
            if self.__Constraint.is_free(i):
                self.__dofNumbers.append(start)
                start += 1
            else:
                self.__dofNumbers.append(-1)
        return start

    def get_dof_numbers(self):
        return self.__dofNumbers

    def get_position(self):
        return array([self.x1, self.x2, self.x3])

    def set_displacement(self, u):
        self.__displacement[0] = u[0]
        self.__displacement[1] = u[1]
        self.__displacement[2] = u[2]

    def get_displacement(self):
        return self.__displacement

    def node_print(self):
        print(self.__position)