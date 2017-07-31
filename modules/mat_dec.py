import numpy as np
class matrix:
    def __init__(self,mat):
        self.mat=mat
        if self.mat.shape[0]!=self.mat.shape[1]:
            raise ValueError('不是方阵！请检查错误')
    def tri_decomposition(self):
        n=self.mat.shape[0]


