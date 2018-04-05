import numpy as np
from numba import jit
class matrix:
    def __init__(self,mat):
        self.mat=mat
        if self.mat.shape[0]!=self.mat.shape[1]:
            raise ValueError('不是方阵！请检查错误')
        elif np.abs(self.det)<1e-8:
            raise ValueError('奇异矩阵！抱歉我算不了^_^!!!')
    def det(self):
        return np.linalg.det(self.mat)
    @jit
    def _tri_decom(self):
        n = self.mat.shape[0]
        A = np.eye(n)
        D = np.copy(self.mat)
        for i in range(n - 1):
            for j in range(i+1,n):
                A[j,i]=D[j,i]/D[i,i]
            for j in range(i + 1, n):
                for jj in range(i + 1, n):
                    D[j, jj] = D[j, jj] - D[j, i] * D[i, jj] / D[i, i]
            for j in range(i+1,n):
                D[j,i]=0.0
                D[i,j]=0.0
        return A, D
    def tri_decomposition(self,origin='other'):
        if origin=='my':
            return self._tri_decom()
        else:
            P=np.linalg.cholesky(self.mat)
            D_tmp=np.diag(np.diag(P)**(-1))
            A=P@D_tmp
            D=np.diag(np.diag(P)**(2))
            return A,D
    def cholesky(self,origin='other'):
        if origin=='my':
            A,D=self.tri_decomposition(origin='my')
            return A@np.sqrt(D)
        else:
            return np.linalg.cholesky(self.mat)



