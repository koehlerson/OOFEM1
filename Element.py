from abc import ABC, abstractmethod
from numpy import outer, block, dot, array, identity
from math import sqrt
from decimal import Decimal


class Element(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def compute_elastic_stiffness_matrix(self):
        pass

    @abstractmethod
    def enumerate_dofs(self):
        pass

    @abstractmethod
    def get_dof_numbers(self):
        pass

    @abstractmethod
    def element_print(self):
        pass

    @abstractmethod
    def update_displacements(self):
        pass


class TrussElementLinear(Element):

    def __init__(self, e, a, n1, n2):
        super().__init__()
        self.n1 = n1
        self.n2 = n2
        self.__dofNumbers = array([self.n1.get_dof_numbers(), self.n2.get_dof_numbers()]).flatten()
        self.__eModulus = e
        self.__area = a
        self.u = 0

    def compute_elastic_stiffness_matrix(self):
        """"Computes special case of tangential matrix"""
        coeff = self.__eModulus * self.__area / self.get_length() ** 3
        a = outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        b = -outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        c = -outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        d = outer(self.n2.get_position() - self.n1.get_position(), self.n2.get_position() - self.n1.get_position())
        return coeff*block([[a, b], [c, d]])

    def get_length(self):
        """Compute the length of the element"""
        vector = self.n2.get_position() - self.n1.get_position()
        return sqrt(dot(vector, vector))

    def enumerate_dofs(self):
        """Enumerates the dofs of the element"""
        self.__dofNumbers = array([self.n1.get_dof_numbers(), self.n2.get_dof_numbers()]).flatten()

    def get_dof_numbers(self):
        """Returns the dofs of the element"""
        return self.__dofNumbers

    def compute_force(self):
        """Computes inner force and returns it as a scalar"""
        u1 = self.n1.get_displacement()
        u2 = self.n2.get_displacement()
        eps = dot(u2-u1, self.get_e1())/self.get_length()
        return self.__eModulus*self.__area*eps

    def get_e1(self):
        """Computes eps scalar of the element"""
        e1 = (self.n2.get_position() - self.n1.get_position()) / \
              self.get_length()
        return e1

    def get_node1(self):
        """Returns first node of the element"""
        return self.n1

    def get_node2(self):
        """Returns second node of the element"""
        return self.n2

    def get_area(self):
        """Returns private property area"""
        return self.__area

    def get_e_modulus(self):
        """Returns private property E-Modulus"""
        return self.__eModulus

    def update_displacements(self):
        """Compute 6x1 displacement vector of the element and updates the zero vector"""
        X = block([self.n1.get_position(), self.n2.get_position()])
        x = X + block([self.n1.get_displacement(), self.n2.get_displacement()])
        self.u = x - X

    def element_print(self):
        """Print routine of the element"""
        print('%.2E' % Decimal(self.__eModulus), "\t", '%.2E' % Decimal(self.__area), "\t", self.get_length())


class CrisfieldTruss(Element):

    def __init__(self, e, a, n1, n2, prestress = 0):
        super().__init__()
        self.__eModulus = e
        self.__area = a
        self.n1 = n1
        self. n2 = n2
        self.u = array([0, 0, 0, 0, 0, 0])
        self.ri = 0
        self.X = block([self.n1.get_position(), self.n2.get_position()])
        self.A = block([[identity(3), -1*identity(3)], [-1*identity(3), identity(3)]])
        self.__dofNumbers = array([self.n1.get_dof_numbers(), self.n2.get_dof_numbers()]).flatten()
        self.s11 = prestress

    def compute_elastic_stiffness_matrix(self):
        """"Computes special case of tangential matrix"""
        coef = (self.__eModulus * self.__area) / self.get_length()
        a = outer(self.n1.get_position()-self.n2.get_position(), self.n1.get_position()-self.n2.get_position())
        b = outer(self.n1.get_position()-self.n2.get_position(), self.n2.get_position()-self.n1.get_position())
        c = outer(self.n2.get_position()-self.n1.get_position(), self.n1.get_position()-self.n2.get_position())
        d = outer(self.n2.get_position()-self.n1.get_position(), self.n2.get_position()-self.n2.get_position())
        return coef*block([[a, b], [c, d]])

    def compute_tangent_matrix(self):
        """Computes the tangent stiffness matrix for nonlinear analysis purposes"""
        coefm = (self.__eModulus * self.__area)/ self.get_length()**3
        coefg = self.__area/self.get_length()
        Km = outer(dot(self.A, (self.X + self.u)),dot((self.X + self.u), self.A))
        Kg = (self.compute_strains()*self.__eModulus + self.s11)*self.A
        return coefm*Km + coefg*Kg

    def get_length(self):
        """Compute the length of the element"""
        vector = self.n2.get_position() - self.n1.get_position()
        return sqrt(dot(vector, vector))

    def compute_strains(self):
        """Compute the Green-Lagrange strain. Scalar Output"""
        e = 1/self.get_length()*dot((self.X + 0.5*self.u), dot(self.A, self.u))
        return e

    def compute_force(self):
        """Computes inner force and returns it as a scalar"""
        u1 = self.n1.get_displacement()
        u2 = self.n2.get_displacement()
        eps = dot(u2-u1, self.get_e1())/self.get_length()
        return self.__eModulus*self.__area*eps

    def compute_internal_force(self):
        """Internal force vector. Outputs a 6x1 Vector"""
        u1 = self.n1.get_displacement()
        u2 = self.n2.get_displacement()
        x1 = self.n1.get_position()
        x2 = self.n2.get_position()
        self.ri = self.__area / self.get_length() * block([x1 + u1 - x2 - u2, -1 * (x1 + u1 - x2 - u2)])
        return self.ri

    def enumerate_dofs(self):
        """Enumerates the dofs of the element"""
        self.__dofNumbers = array([self.n1.get_dof_numbers(), self.n2.get_dof_numbers()]).flatten()

    def get_dof_numbers(self):
        """Returns the dofs of the element"""
        return self.__dofNumbers

    def get_node1(self):
        """Returns first node of the element"""
        return self.n1

    def get_node2(self):
        """Returns second node of the element"""
        return self.n2

    def get_area(self):
        """Returns private property area"""
        return self.__area

    def get_e_modulus(self):
        """Returns private property E-Modulus"""
        return self.__eModulus

    def get_e1(self):
        """Computes eps scalar of the element"""
        e1 = (self.n2.get_position() - self.n1.get_position()) / \
              self.get_length()
        return e1

    def update_displacements(self):
        """Compute 6x1 displacement vector of the element and updates the zero vector"""
        x = self.X + block([self.n1.get_displacement(), self.n2.get_displacement()])
        self.u = x - self.X

    def element_print(self):
        """Print routine of the element"""
        print('%.2E' % Decimal(self.__eModulus), "\t", '%.2E' % Decimal(self.__area), "\t", self.get_length())
