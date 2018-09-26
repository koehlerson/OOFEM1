class Constraint:

    def __init__(self, dof1=True, dof2=True, dof3=True):
        self.dof1 = dof1
        self.dof2 = dof2
        self.dof3 = dof3
        self.__free = [self.dof1, self.dof2, self.dof3]

    def constraint(self, u1, u2, u3):
        """constraining the nodal constraint class with booleans"""
        self.dof1 = u1
        self.dof2 = u2
        self.dof3 = u3

    def is_free(self, i):
        """check if the i-th entry is free, if so returns true, else false"""
        return self.__free[i]

    def constraint_print(self):
        """prints basic information about the nodal constraints instance"""
        u1, u2, u3 = "fixed", "fixed", "fixed"
        if self.is_free(0):
            u1 = "free"
        if self.is_free(1):
            u2 = "free"
        if self.is_free(2):
            u3 = "free"
        return print(u1, "\t", u2, "\t", u3)
