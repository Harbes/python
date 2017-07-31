import numpy as np
class matrix:
    def __init__(self,mat):
        self.mat=mat
        if self.mat.shape[0]!=self.mat.shape[1]:
            raise ValueError('不是方阵！请检查错误')
        elif np.abs(self.det)<1e-8:
            raise ValueError('奇异矩阵！抱歉我算不了^_^!!!')
    @property
    def det(self):
        return np.linalg.det(self.mat)
    def tri_decomposition(self,origin='other'):
        if origin=='my':
            n=self.mat.shape[0]
            A=np.eye(n)
            D=self.mat
            for i in range(n-1):
                E=np.eye(n)
                E[i+1:,i]=-D[i+1:,i]/D[i,i]
                D=E@D@E.T
                E[i + 1:, i] *=-1
                A=A@E
            return A,D
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



