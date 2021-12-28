import numpy as np
from scipy.sparse.linalg import svds

class ISTEMC:
    def __init__(self):
        pass

    def SoftTh(self,s, thld):
        z = np.multiply(np.sign(s), np.maximum(abs(s) - thld, np.zeros(s.shape)))
        return z

    def opRestriction(self,n, idx, matrix, mode, matrix1):
        # since python start array index from 0. We compare 0 and n-1
        if (np.min(idx) < 0 or np.max(idx) > n - 1):
            print("Index parameter must be integer and match dimensions of the operator")
        op = self.opRestriction_intrnl(n, idx, matrix, mode, matrix1)
        return op

    def opRestriction_intrnl(self,n, idx, x, mode, matrix1):
        m = len(idx)
        # self.checkDimensions(m,n,x,mode)
        if mode == 0:
            y = [m, n, [0, 1, 0, 1], {'Restriction'}]
        elif mode == 1:
            y = matrix1[idx]
        else:
            #         y = np.zeros((n,1))
            y = np.zeros(shape=(n))
            y[idx] = x
        return y

    def checkDimensions(self,m, n, x, mode):
        if (mode == 0):
            return 0
        if (x.shape[0] == 1 and x.shape[1] == 1):
            if (m != 1 or n != 1):
                print("Operator-scalar multiplication not yet supported")
        if (mode == 1):
            if (x.shape[0] != n):
                print("Incompatible dimensions")
            if (x.shape[1] != 1):
                print("Operator-matrix multiplication not yet supported")
        else:
            if (x.shape[1] != m):
                print("Incompatible dimensions")
            if (x.shape[1] != 1):
                print("Operator-matrix multiplication not yet supported")
        return 0

    def IST_eMC(self,y, sizeX, rankX, N1, IDX1):
        err = 1e-12
        # x_initial = np.zeros(shape=(np.prod(sizeX)))
        normfac = 1
        insweep = 20
        tol = 1e-4
        decfac = 0.7
        alpha = 1.1*normfac
        x = np.zeros(shape=(np.prod(sizeX)))
        lambdaInit = decfac*np.max(abs(self.opRestriction(N1,IDX1,y,2,np.ravel(y))))
        lambdaa = lambdaInit
        f_current = np.linalg.norm(y-self.opRestriction(N1,IDX1,x,1,np.ravel(x))) + lambdaa*np.linalg.norm(x,1)
        while(lambdaa >lambdaInit*tol):
            for ins in range(insweep):
                f_previous = f_current
                z = self.opRestriction(N1,IDX1,x,1,np.ravel(x))
                U,S,V = svds(np.reshape(x + (1/alpha)*self.opRestriction(N1,IDX1,y-z,2,np.ravel(y-z)),(sizeX),order='F'), k = rankX, which = 'LM')
                del z
                S = np.diag(self.SoftTh(S,lambdaa/(2*alpha)))
                X = np.matmul(np.matmul(U,S),V)
                del S
                del U
                del V
                X[X<0] = 0
                x = np.ravel(X)
                res = self.opRestriction(N1,IDX1,x,1,np.ravel(x))
                f_current = np.linalg.norm(y - res) + lambdaa*np.linalg.norm(x,1)
                if (np.linalg.norm(f_current - f_previous)/np.linalg.norm(f_current + f_previous)<tol):
                    break
            if (np.linalg.norm(y-res)< err):
                break
            lambdaa = decfac*lambdaa
        return X

    def istCalc(self,matrix, rank):
        # matrix1 = np.ravel(matrix, order='F')
        # N1 = len(matrix1)
        # t1 = list(range(N1))
        # IDX1 = np.where(matrix1 > 0)
        # y1 = self.opRestriction(N1, IDX1, matrix, 1, matrix1)
        Z = self.IST_eMC(self.opRestriction(len(np.ravel(matrix, order='F')), np.where(np.ravel(matrix, order='F') > 0) , matrix, 1, np.ravel(matrix, order='F')), matrix.shape,
                        rank,len(np.ravel(matrix, order='F')) , np.where(np.ravel(matrix, order='F') > 0))
        return Z
