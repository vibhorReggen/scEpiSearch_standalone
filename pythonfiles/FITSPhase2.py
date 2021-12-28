import numpy as np
import glob
from scipy.stats import spearmanr

class FITSPhase2:
    def __init__(self, dataX, topk, colWise, name2save):
        self.startProcess( dataX, topk, colWise, name2save)

    def startProcess(self, dataX, topk, colWise, name2save):
        dataX = np.log(dataX + 1.01)
        files = name2save+'_*.npy'
        trees = {}
        start = 1
        for t in glob.glob(files):
            try:
                trees[start] = {}
                final_imputed = np.load(t)
                trees[start]['val'] = final_imputed
                start = start + 1
            except:
                pass
                # print(name2save + "_" + str(t) + ".npy not exists")
        for t in range(3, 7):
            try:
                trees[start] = {}
                final_imputed = np.load(name2save + "_r_" + str(t) + ".npy")
                trees[start]['val'] = final_imputed
                start = start + 1
            except:
                # print(name2save + "_r_" + str(t) + ".npy not exists")
                pass
        if start < 3:
            topk = start
        self.maxCorrelated(dataX, trees, start - 1, topk, colWise, name2save)

    def maxCorrelated(self,mOriginal, mtree, count, topk, colWise, name2save):
        if colWise == 1:
            final_imputed = self.maxCorrelatedCol(mOriginal, mtree, count, topk)
            final_imputed = np.transpose(final_imputed)
            np.savetxt(name2save + "_feature.txt", final_imputed, delimiter=",")
        else:
            final_imputed = self.maxCorrelatedRow(mOriginal, mtree, count, topk)
            final_imputed = np.transpose(final_imputed)
            np.savetxt(name2save + "_sample.txt", final_imputed, delimiter=",")

    def maxCorrelatedCol(self,mOriginal, mtree, count, topk):
        res = np.zeros((mOriginal.shape[0], mOriginal.shape[1]))
        for i in range(mOriginal.shape[1]):
            corrAll = []
            for j in range(1, count + 1):
                cor, pval = spearmanr(mtree[j]['val'][:, i], mOriginal[:, i])
                corrAll = np.append(corrAll, cor)
            indices = np.argsort(corrAll)[::-1]
            corrAll = corrAll[indices]
            newMatrix = np.zeros((mOriginal.shape[0], topk))
            for j in range(topk):
                newMatrix[:, j] = mtree[indices[j] + 1]['val'][:, i]
            res[:, i] = np.max(newMatrix, axis=1).T
        return res

    def maxCorrelatedRow(self,mOriginal, mtree, count, topk):
        res = np.zeros((mOriginal.shape[0], mOriginal.shape[1]))
        for i in range(mOriginal.shape[0]):
            corrAll = []
            for j in range(1, count + 1):
                cor, pval = spearmanr(mtree[j]['val'][i, :], mOriginal[i, :])
                corrAll = np.append(corrAll, cor)
            indices = np.argsort(corrAll)[::-1]
            corrAll = corrAll[indices]
            newMatrix = np.zeros((topk, mOriginal.shape[1]))
            for j in range(topk):
                newMatrix[j, :] = mtree[indices[j] + 1]['val'][i, :]
            res[i, :] = np.max(newMatrix, axis=0)
        return res
