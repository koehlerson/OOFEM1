from Node import *
from Element import *
from numpy.linalg import solve, norm
from numpy import array, round, zeros
import numpy
from matplotlib import pyplot as plt
from decimal import Decimal


class Structure:

    def __init__(self, t = 0):
        """input of the time is a list of given time steps"""
        self.Elements = []
        self.Nodes = []
        self.u = []
        if type(t) == numpy.ndarray or type(t) == list:
            self.dt = t[1] - t[0]
            self.t_n = len(t)
        self.t = t
        self.t_cur = 0
        self.analysis_type = 'linear'

    def add_node(self, x1, x2, x3):
        """adds a node by 3 given coordinates x1,x2,x3"""
        self.Nodes.append(Node(x1, x2, x3))
        return self.Nodes[-1]

    def add_element(self, e):
        """add element e"""
        self.Elements.append(e)
        return self.Elements[-1]

    def get_number_of_nodes(self):
        """"returns the # of dofs"""
        return len(self.Nodes)

    def get_element(self, i):
        return self.Elements[i]

    def print_structure(self, dof = 1):
        """prints the information of the BVP as well as the analysis results"""
        counter = 0
        print("====================================================")
        print("====================================================")
        print("Problem Description: ")
        print("idx\t", "x1 \t\t\t", "x2 \t\t\t", "x3")
        for i in self.Nodes:
            print(counter, "\t",
                  '%.4E' % Decimal(i.get_position()[0]), "\t",
                  '%.4E' % Decimal(i.get_position()[1]), "\t",
                   '%.4E' % Decimal(i.get_position()[2]))
            counter += 1
        print("====================================================")
        print("idx\t", "u1   \t", "u2   \t", "u3   \t")
        for i in self.Nodes:
            print(counter, end='\t')
            i.get_constraint().constraint_print()
            counter += 1
        print("====================================================")
        print("idx\t", "   r1 \t\t", "    r2 \t\t", "r3")
        counter = 0
        for i in self.Nodes:
            print(counter, end='     \t')
            print('%.2E' % Decimal(i.get_force().get_component(0)), "\t",
                  '%.2E' % Decimal(i.get_force().get_component(1)), "\t",
                  '%.2E' % Decimal(i.get_force().get_component(2)))
        print("====================================================")
        counter = 0
        print("idx", "E\t\t\t" " A\t\t\t", "length")
        for i in self.Elements:
            print(counter, end='\t')
            i.element_print()
            counter += 1
        counter = 0
        print("====================================================")
        print("====================================================")
        print(" ")
        print("Analysis results")
        print("idx\t\t", "u1 \t\t\t", "u2 \t\t\t", "u3")
        for i in self.Nodes:
            print(counter,"\t\t",
                  '%.4E' % Decimal(i.get_displacement()[0]), "\t",
                  '%.4E' % Decimal(i.get_displacement()[1]), "\t"
                                                             '%.4E' % Decimal(i.get_displacement()[2]))
            counter += 1
        print("====================================================")

    def linear_solve(self):
        """When linear elements used, solves the system in a linear manner by calling numpy.solve()"""
        self.t_n = 1
        self.t_cur = 0
        eqn = self.__enumerate_dofs()
        K = zeros([eqn, eqn])
        f = zeros([eqn])
        K = self.__assemble_linear_stiffness_matrix(K)
        f = self.__assemble_external_load_vector(f)
        self.u.append((solve(K, f)))
        self.__select_displacements()
        return self.u[0]

    def nonlinear_solve(self):
        """When nonlinear elements used, method will use a Newton-Raphson scheme in order to find the equilibrium"""
        self.analysis_type = 'nonlinear'
        if all(isinstance(i, CrisfieldTruss) for i in self.Elements):
            eqn = self.__enumerate_dofs()
            K = zeros([eqn, eqn])
            f_ext = zeros([eqn])
            f_int = zeros([eqn])
            u = zeros([eqn])      #initial displacements
            du = zeros([eqn])     #initial increments
            G = zeros([eqn])
            for t in self.t:
                f_ext = self.__assemble_external_load_vector(f_ext)
                f_ext = f_ext*t
                niter = 0
                du = array([100,100])
                while norm(du, 2) > 1e-8 and niter < 200:
                    K = self.__assemble_nonlinear_stiffness_matrix(K)
                    f_int = self.__assemble_internal_load_vector(f_int)
                    G = f_ext - f_int
                    du = solve(K, G)
                    u += du
                    niter += 1
                self.u.append(u)
                self.__select_displacements()
                self.__update_time()
            return self.u[-1]
        else:
            print("Linear element chosen. Choose another one and start Program again")
            pass

    def __enumerate_dofs(self):
        """enumerates all dofs. blocked dofs are setted to -1"""
        start = 0
        for i in self.Nodes:
            start = i.enumerate_dofs(start)
        for i in self.Elements:
            i.enumerate_dofs()
        return start

    def __assemble_external_load_vector(self, rglob):
        """assembles the external load vector. overwrites the input vector, but returns it as well"""
        for i in self.Nodes:
            for j in range(3):
                if i.get_constraint().is_free(j):
                   rglob[i.get_dof_numbers()[j]] += i.get_force().get_component(j)
        return rglob

    def __assemble_linear_stiffness_matrix(self, kglob):
        """assemble the stiffness matrix in the linear case"""
        for i in self.Elements:
            ke = i.compute_elastic_stiffness_matrix()
            dofs = i.get_dof_numbers()
            for j in range(6):
                for k in range(6):
                    if (dofs[j] != -1) and (dofs[k] != -1):
                        kglob[dofs[j], dofs[k]] += ke[j, k]
        return kglob

    def __assemble_nonlinear_stiffness_matrix(self, kglob):
        """assembles the tangential stiffness matrix"""
        for i in self.Elements:
            ke = i.compute_tangent_matrix()
            dofs = i.get_dof_numbers()
            for j in range(6):
                for k in range(6):
                    if (dofs[j] != -1) and (dofs[k] != -1):
                        kglob[dofs[j], dofs[k]] += ke[j, k]
        return kglob

    def __assemble_internal_load_vector(self, rglob):
        """assembles the internal load vector"""
        for i in self.Elements:
            f_int = i.compute_internal_force()
            dofs = i.get_dof_numbers()
            for j in range(6):
                if (dofs[j] != -1):
                    rglob[dofs[j]] += f_int[j]
        return rglob

    def __select_displacements(self):
        """updates the displacement for every element and every node in the structure"""
        for i in self.Nodes:
            DOFs = i.get_dof_numbers()
            u = []
            for j in DOFs:
                if j == -1:
                    u.append(0)
                else:
                    u.append(self.u[int(self.t_cur * self.t_n)][j])
            i.set_displacement(u)
        for i in self.Elements:
            i.update_displacements()

    def __update_time(self):
        """update the current time"""
        self.t_cur += self.dt

    def print_results(self):
        """Prints only the results of the structure"""
        print("Global displacement vector u: ")
        print(self.u[0])
        print(" ")
        for i in self.Elements:
            print("Element force: ", round(i.compute_force(), 3))
