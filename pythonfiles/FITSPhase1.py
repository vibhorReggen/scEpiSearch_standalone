from ISTEMC import ISTEMC
import math
import numpy as np
# from sklearn.cluster import KMeans
from scipy.cluster.vq import kmeans2
# from sklearn.manifold import TSNE
from MulticoreTSNE import MulticoreTSNE as TSNE
import random
from scipy.sparse.linalg import svds

class FITSPhase1:
    def __init__(self,dataX,maxClusters,msc,aRFSR,maxAllowedLevel,name2save,num):
        self.istObj = ISTEMC()
        self.startProcess(dataX, maxClusters, msc, aRFSR, maxAllowedLevel, name2save, num)

    def startProcess(self,dataX,maxClusters,msc,aRFSR,maxAllowedLevel,name2save,num):
        dataX = np.log(dataX + 1.01)
        rank = random.randint(3, 6)
        impData = self.istObj.istCalc(dataX, rank)
        np.save(name2save + "_r_" + str(rank) + ".npy", impData)
        RFSR = aRFSR
        level = 1
        flag = 1
        clusters = {}
        records = {}
        cluster = maxClusters
        allowedLevel = random.randint(2, maxAllowedLevel)
        while (level <= allowedLevel and flag != 0):
            clusters[level] = {}
            records[level] = {}
            clusters[level]['val'] = self.clusterCount(cluster)
            if level == 1:
                while (1):
                    records[level]['predictedLabel'] = self.kclusters(impData, clusters[level]['val'], RFSR)
                    flag = self.checkClustersMinCount(records[level]['predictedLabel'], clusters[level]['val'], msc)
                    if flag != 0:
                        break
                    else:
                        clusters[level]['val'] = clusters[level]['val'] - 1

                records[level]['c'] = {}
                records[level]['unimpdata'] = {}
                for i in range(clusters[level]['val']):
                    records[level]['c'][i] = {}
                    records[level]['unimpdata'][i] = {}
                    records[level]['c'][i]['val'] = np.where(records[level]['predictedLabel'] == i)[0]
                    records[level]['unimpdata'][i]['val'] = dataX[records[level]['c'][i]['val'], :]

                records[level]['vec'] = {}
                records[level]['data'] = {}
                for i in range(clusters[level]['val']):
                    records[level]['vec'][i] = {}
                    records[level]['data'][i] = {}
                    records[level]['vec'][i]['val'] = np.where(~records[level]['unimpdata'][i]['val'].any(axis=0))[0]
                    records[level]['data'][i]['val'] = records[level]['unimpdata'][i]['val']
                    np.delete(records[level]['data'][i]['val'], records[level]['vec'][i]['val'], axis=1)

                for i in range(clusters[level]['val']):
                    records[level]['data'][i]['val'], records[level]['vec'][i]['val'] = self.operationOnMatrix(
                        records[level]['data'][i]['val'], rank, impData, records[level]['c'][i]['val'],
                        records[level]['vec'][i]['val'], records[level]['vec'][i]['val'])
            else:
                records[level]['predictedLabel'] = self.clusterPrediction(1, level, clusters, records[level - 1]['data'],
                                                                     RFSR)

                records[level]['c'] = self.clusterMaking(1, level, clusters, records[level]['predictedLabel'])
                records[level]['unimpdata'] = self.rawDataClusterMaking(1, level, clusters, records[level - 1]['unimpdata'],records[level]['c'])
                records[level]['vec'] = self.vecFormation(1, level, clusters, records[level]['unimpdata'])
                records[level]['data'] = self.dataFormation(1, level, clusters, records[level]['unimpdata'],records[level]['vec'])

                records[level]['data'], records[level]['vec'] = self.imputationOperation(rank, 1, level, clusters,records[level]['data'],
                                                                                    records[level - 1]['data'],records[level]['vec'],
                                                                                    records[level - 1]['vec'],records[level]['c'])

            cluster = math.ceil(clusters[level]['val'] / 2)
            if cluster < 2:
                cluster = 2
            if level-2 >0:
                records[level-2]['predictedLabel'] = {}
                records[level-2]['val'] = {}
                records['unimpdata'] = {}
            level = level + 1

        level = level - 1
        last = level
        prev = last
        while (level >= 0):
            if last == level:
                records[level]['newdata'] = self.newdataInitialization(1, level, clusters, records[level]['data'], dataX)
                records[level]['vec_negation'] = self.vec_negationInitialization(1, level, clusters,records[level]['vec'], dataX)
                records[level]['newdata'] = self.newDataMaking(1, level, clusters, records[level]['newdata'],records[level]['vec_negation'], records[level]['data'])

            elif level == 0:
                pre_final_imputed_new = np.zeros((dataX.shape))
                for i in range(clusters[1]['val']):
                    pre_final_imputed_new[records[prev]['c'][i]['val'], :] = records[prev]['newdata'][i]['val']
            else:
                records[level]['newdata'] = {}
                records[level]['newdata'] = self.newdataMerging(1, prev, clusters, records[prev]['newdata'],records[level]['data'], records[prev]['c'], dataX)
            prev = level
            level = level - 1
        final_imputed = np.zeros((dataX.shape))
        final_imputed = pre_final_imputed_new
        name2save = name2save + "_" + str(num) + ".npy"
        np.save(file=name2save, arr=final_imputed)

    def clusterCount(self,maxk):
        val = np.random.permutation(maxk)
        v = val[0]
        if v < 2:
            v = 2
        return v

    def kclusters(self,matrix, k, RFSR):
        random_vec = np.random.permutation(matrix.shape[1])
        r = random.randint(RFSR[0], RFSR[1]) / 100
        U, S, V = svds(matrix[:, random_vec[0:int(matrix.shape[1] * r)]], k=10)
        # matrix = TSNE(n_components=3).fit_transform(U)
        matrix = TSNE(n_jobs=2).fit_transform(U)
        loc = np.random.permutation(matrix.shape[0])[:k]
        init = matrix[loc, :]
        # Kmeans = KMeans(n_clusters=k, random_state=0, init=init).fit(matrix)
        Kmeans1, Kmeans2 = kmeans2(matrix, k= init)
        # return Kmeans.labels_
        return Kmeans2

    def operationOnMatrix(self,matrix, rank, impmatrix, rows, newVec, oldVec):
        matrix = np.array(matrix)
        if (matrix.shape[0] >= rank):
            toreturn = self.istObj.istCalc(matrix, rank)
            updatedvec = newVec
        else:
            toreturn = impmatrix[rows, :]
            updatedvec = oldVec
        return toreturn, updatedvec

    def clusterPrediction(self,startLevel, endLevel, clusters, data, RFSR):
        predictedLabel = {}
        for i in range(clusters[startLevel]['val']):
            predictedLabel[i] = {}
            if startLevel == endLevel - 1:
                predictedLabel[i]['val'] = self.kclusters(data[i]['val'], clusters[endLevel]['val'], RFSR)
            else:
                predictedLabel[i]['predictedLabel'] = {}
                predictedLabel[i]['predictedLabel'] = self.clusterPrediction(startLevel + 1, endLevel, clusters,
                                                                        data[i]['data'], RFSR)
        return predictedLabel

    def clusterMaking(self,startLevel, endLevel, clusters, predictedLabel):
        c = {}
        for i in range(clusters[startLevel]['val']):
            c[i] = {}
            c[i]['c'] = {}
            if startLevel == endLevel - 1:
                for m in range(clusters[endLevel]['val']):
                    c[i]['c'][m] = {}
                    c[i]['c'][m]['val'] = np.where(predictedLabel[i]['val'] == m)
            else:
                c[i]['c'] = self.clusterMaking(startLevel + 1, endLevel, clusters, predictedLabel[i]['predictedLabel'])
        return c

    def rawDataClusterMaking(self,startLevel, endLevel, clusters, unimpdata, c):
        unimpdata1 = {}
        for i in range(clusters[startLevel]['val']):
            unimpdata1[i] = {}
            if startLevel == endLevel - 1:
                unimpdata1[i]['unimpdata'] = {}
                for m in range(clusters[endLevel]['val']):
                    unimpdata1[i]['unimpdata'][m] = {}
                    unimpdata1[i]['unimpdata'][m]['val'] = {}
                    unimpdata1[i]['unimpdata'][m]['val'] = unimpdata[i]['val'][np.ravel(c[i]['c'][m]['val']), :]
            else:
                unimpdata1[i]['unimpdata'] = self.rawDataClusterMaking(startLevel + 1, endLevel, clusters,
                                                                  unimpdata[i]['unimpdata'], c[i]['c'])
        return unimpdata1

    def vecFormation(self,startLevel, endLevel, clusters, unimpdata):
        vec = {}
        for i in range(clusters[startLevel]['val']):
            vec[i] = {}
            vec[i]['vec'] = {}
            if startLevel == endLevel - 1:
                for m in range(clusters[endLevel]['val']):
                    vec[i]['vec'][m] = {}
                    vec[i]['vec'][m]['val'] = {}
                    vec[i]['vec'][m]['val'] = np.where(~unimpdata[i]['unimpdata'][m]['val'].any(axis=0))
            else:
                vec[i]['vec'] = self.vecFormation(startLevel + 1, endLevel, clusters, unimpdata[i]['unimpdata'])
        return vec

    def dataFormation(self,startLevel, endLevel, clusters, unimpdata, vec):
        data = {}
        for i in range(clusters[startLevel]['val']):
            data[i] = {}
            if startLevel == endLevel - 1:
                data[i]['data'] = {}
                for m in range(clusters[endLevel]['val']):
                    data[i]['data'][m] = {}
                    data[i]['data'][m]['val'] = unimpdata[i]['unimpdata'][m]['val']
                    if list(np.ravel(vec[i]['vec'][m]['val'])):
                        np.delete(data[i]['data'][m]['val'], list(np.ravel(vec[i]['vec'][m]['val'])), axis=1)
            else:
                data[i]['data'] = self.dataFormation(startLevel + 1, endLevel, clusters, unimpdata[i]['unimpdata'],
                                                vec[i]['vec'])
        return data

    def imputationOperation(self,rank, startLevel, endLevel, clusters, data, impdata, vec, oldvec, c):
        for i in range(clusters[startLevel]['val']):
            if startLevel == endLevel - 1:
                for m in range(clusters[endLevel]['val']):
                    data[i]['data'][m]['val'], vec[i]['vec'][m]['val'] = self.operationOnMatrix(data[i]['data'][m]['val'],
                                                                                           rank, impdata[i]['val'],
                                                                                           c[i]['c'][m]['val'],
                                                                                           vec[i]['vec'][m]['val'],
                                                                                           oldvec[i]['val'])
            else:
                data[i]['data'], vec[i]['vec'] = self.imputationOperation(rank, startLevel + 1, endLevel, clusters,
                                                                     data[i]['data'], impdata[i]['data'], vec[i]['vec'],
                                                                     oldvec[i]['vec'], c[i]['c'])
        return data, vec

    def newdataInitialization(self,startLevel, endLevel, clusters, data, dataX):
        newdata = {}
        for i in range(clusters[startLevel]['val']):
            newdata[i] = {}
            newdata[i]['newdata'] = {}
            if startLevel == endLevel - 1:
                for m in range(clusters[endLevel]['val']):
                    newdata[i]['newdata'][m] = {}
                    newdata[i]['newdata'][m]['val'] = {}
                    newdata[i]['newdata'][m]['val'] = np.zeros((data[i]['data'][m]['val'].shape[0], dataX.shape[1]))
            else:
                newdata[i]['newdata'] = self.newdataInitialization(startLevel + 1, endLevel, clusters, data[i]['data'],
                                                              dataX)
        return newdata

    def vec_negationInitialization(self,startLevel, endLevel, clusters, vec, dataX):
        vec_negation = {}
        for i in range(clusters[startLevel]['val']):
            vec_negation[i] = {}
            vec_negation[i]['vec_negation'] = {}
            if startLevel == endLevel - 1:
                for m in range(clusters[endLevel]['val']):
                    vec_negation[i]['vec_negation'][m] = {}
                    vec_negation[i]['vec_negation'][m]['val'] = {}
                    vec_negation[i]['vec_negation'][m]['val'] = list(range(dataX.shape[1]))
                    if np.ravel(vec[i]['vec'][m]['val']):
                        np.delete(vec_negation[i]['vec_negation'][m]['val'], np.ravel(vec[i]['vec'][m]['val']))
            else:
                vec_negation[i]['vec_negation'] = self.vec_negationInitialization(startLevel + 1, endLevel, clusters,
                                                                             vec[i]['vec'], dataX)
        return vec_negation

    def newDataMaking(self,startLevel, endLevel, clusters, newdata, vec_negation, data):
        for i in range(clusters[startLevel]['val']):
            if startLevel == endLevel - 1:
                for m in range(clusters[endLevel]['val']):
                    newdata[i]['newdata'][m]['val'][:, vec_negation[i]['vec_negation'][m]['val']] = data[i]['data'][m][
                        'val']
            else:
                newdata[i]['newdata'] = self.newDataMaking(startLevel + 1, endLevel, clusters, newdata[i]['newdata'],
                                                      vec_negation[i]['vec_negation'], data[i]['data'])
        return newdata

    def newdataMerging(self,startLevel, endLevel, clusters, newdata, data, c, dataX):
        newdata1 = {}
        for i in range(clusters[startLevel]['val']):
            newdata1[i] = {}
            if startLevel == endLevel - 1:
                newdata1[i]['val'] = np.zeros((data[i]['val'].shape[0], dataX.shape[1]))
                for m in range(clusters[endLevel]['val']):
                    newdata1[i]['val'][np.transpose(np.ravel(c[i]['c'][m]['val'])), :] = newdata[i]['newdata'][m]['val']
            else:
                newdata1[i]['newdata'] = {}
                newdata1[i]['newdata'] = self.newdataMerging(startLevel + 1, endLevel, clusters, newdata[i]['newdata'],
                                                        data[i]['data'], c[i]['c'], dataX)
        return newdata1

    def checkClustersMinCount(self,predicted_labels,clusters,msc):
        flag=1
        if clusters>1:
            for i in range(clusters):
                val=np.where(predicted_labels==i)
                if len(val)<msc:
                    flag=0
        return flag
