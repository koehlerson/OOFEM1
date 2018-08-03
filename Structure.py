from Node import *
from Element import *
from numpy.linalg import solve
from numpy import array, round
from decimal import Decimal


class Structure:

    def __init__(self):
        self.Elements = []
        self.Nodes = []
        self.u = array([])

    def add_node(self, x1, x2, x3):
        self.Nodes.append(Node(x1, x2, x3))
        return self.Nodes[-1]

    def add_element(self, e, a, n1, n2):
        self.Elements.append(Element(e, a, self.Nodes[n1], self.Nodes[n2]))
        return self.Elements[-1]

    def get_number_of_nodes(self):
        return len(self.Nodes)

    def get_element(self, i):
        return self.Elements[i]

    def print_structure(self):
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

    def solve(self):
        eqn = self.__enumerate_dofs()
        A = zeros([eqn, eqn])
        b = zeros([eqn])
        A = self.__assemble_stiffness_matrix(A)
        b = self.__assemble_load_vector(b)
        self.u = solve(A, b)
        self.__select_displacements()
        return self.u

    def __enumerate_dofs(self):
        start = 0
        for i in self.Nodes:
            start = i.enumerate_dofs(start)
        for i in self.Elements:
            i.enumerate_dofs()
        return start

    def __assemble_load_vector(self, rglob):
        for i in self.Nodes:
            for j in range(3):
                if i.get_constraint().is_free(j):
                   rglob[i.get_dof_numbers()[j]] += i.get_force().get_component(j)
        return rglob

    def __assemble_stiffness_matrix(self, kglob):
        for i in self.Elements:
            ke = i.compute_stiffness_matrix()
            dofs = i.get_dof_numbers()
            for j in range(6):
                for k in range(6):
                    if (dofs[j] != -1) and (dofs[k] != -1):
                        kglob[dofs[j], dofs[k]] += ke[j, k]
        return kglob

    def __select_displacements(self):
        for i in self.Nodes:
            DOFs = i.get_dof_numbers()
            u = []
            for j in DOFs:
                if j == -1:
                    u.append(0)
                else:
                    u.append(self.u[j])
            i.set_displacement(u)

    def print_results(self):
        print("Global displacement vector u: ")
        print(self.u)
        print(" ")
        for i in self.Elements:
            print("Element force: ", round(i.compute_force(),3))