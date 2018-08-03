from numpy import outer, block, transpose, dot, array, sum
from numpy.linalg import norm
from math import sqrt
from decimal import Decimal


class Element:

    def __init__(self, e, a, n1, n2):
        self.n1 = n1
        self.n2 = n2
        self.__dofNumbers = array([self.n1.get_dof_numbers(), self.n2.get_dof_numbers()]).flatten()
        self.__eModulus = e
        self.__area = a

    def compute_stiffness_matrix(self):
        coeff = self.__eModulus * self.__area / self.get_length() ** 3
        a = outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        b = -outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        c = -outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        d = outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        return coeff*block([[a, b], [c, d]])

    def get_length(self):
        vector = self.n2.get_position() - self.n1.get_position()
        return sqrt(dot(vector, vector))

    def enumerate_dofs(self):
        self.__dofNumbers = array([self.n1.get_dof_numbers(), self.n2.get_dof_numbers()]).flatten()

    def get_dof_numbers(self):
        return self.__dofNumbers

    def compute_force(self):
        u1 = self.n1.get_displacement()
        u2 = self.n2.get_displacement()
        eps = dot(u2-u1, self.get_e1())/self.get_length()
        return self.__eModulus*self.__area*eps

    def get_e1(self):
        e1 = (self.n2.get_position() - self.n1.get_position()) / \
              self.get_length()
        return e1

    def get_node1(self):
        return self.n1

    def get_node2(self):
        return self.n2

    def get_area(self):
        return self.__area

    def get_e_modulus(self):
        return self.__eModulus

    def element_print(self):
        print('%.2E' % Decimal(self.__eModulus), "\t", '%.2E' % Decimal(self.__area), "\t", self.get_length())