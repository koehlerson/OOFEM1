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
        """set the constraints of the Node to a given constraint instance"""
        self.__Constraint = Constraint

    def get_constraint(self):
        """returns all constraints as a list"""
        return self.__Constraint

    def set_force(self, force):
        """sets the nodal force vector to a specific force class instance"""
        self.__Force = force

    def get_force(self):
        """gets the force instance which is related to the node"""
        return self.__Force

    def enumerate_dofs(self, start):
        """enumerate the dofs, input is the start of enumeration and output the last number of dofs"""
        for i in range(3):
            if self.__Constraint.is_free(i):
                self.__dofNumbers.append(start)
                start += 1
            else:
                self.__dofNumbers.append(-1)
        return start

    def get_dof_numbers(self):
        """gets the dof numbers of the node"""
        return self.__dofNumbers

    def get_position(self):
        """returns the spatial position of the node"""
        return array([self.x1, self.x2, self.x3])

    def set_displacement(self, u):
        """sets the displacement of the node to a given vector/list u"""
        self.__displacement[0] = u[0]
        self.__displacement[1] = u[1]
        self.__displacement[2] = u[2]

    def get_displacement(self):
        """returns the 3x1 displacement vector/list of the node"""
        return self.__displacement

    def node_print(self):
        """prints basic information about the node"""
        print(self.__position)