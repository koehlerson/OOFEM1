class Force:

    def __init__(self, r1=0, r2=0, r3=0):
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.__components = [self.r1, self.r2, self.r3]

    def get_component(self, i):
        """returns the i-th component of the nodal force vector"""
        return self.__components[i]
    
    def force_print(self):
        """prints basic information about the force instance"""
        print(self.__components)