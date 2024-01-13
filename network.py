import numpy as np
from scipy import sparse,spatial
import os
import datetime
import time

class Network(object):
    def __init__(self):
        self.dateList = None
        self.baselineDict = {}
        self.pairsDates = []

    def computeBaseline(self, refDir, secDir):
        import isce
        import isceobj
        from mroipac.baseline.Baseline import Baseline
        import shelve

        #try:
        #    mdb = shelve.open( os.path.join(refDir, 'data'), flag='r')
        #except:
        with shelve.open( os.path.join(refDir, 'raw'), flag='r') as mdb:
            rFrame = mdb['frame']

        #try:
        #    sdb = shelve.open( os.path.join(secDir, 'data'), flag='r')
        #except:
        with shelve.open( os.path.join(secDir, 'raw'), flag='r') as sdb:
            sFrame = sdb['frame']


        bObj = Baseline()
        bObj.configure()
        bObj.wireInputPort(name='masterFrame', object=rFrame)
        bObj.wireInputPort(name='slaveFrame', object=sFrame)

        bObj.baseline()

        print('Baseline at top/bottom: %f %f'%(bObj.pBaselineTop,bObj.pBaselineBottom))

        avgBaseline = (bObj.pBaselineTop + bObj.pBaselineBottom) / 2.0

        return avgBaseline

    def get_baselineDict(self, dataDir):
        refDate = self.dateList[0]
        refDir = os.path.join(dataDir, refDate)
        self.baselineDict[refDate] = 0.0
        for d in self.dateList[1:]:
            secDir = os.path.join(dataDir, d)
            baseline = self.computeBaseline(refDir, secDir)
            self.baselineDict[d] = baseline

    def geometrical_coherence(self, rangeSpacing, theta, wvl, r):

        numDates = len(self.dateList)
        coh = np.zeros((numDates, numDates))

        for i in range(numDates - 1):
            for j in range(i+1,numDates):

                B = np.abs(self.baselineDict[self.dateList[j]] -
                            self.baselineDict[self.dateList[i]])

                corr = 1.0 - 2*B*rangeSpacing*((np.cos(theta))**2)/wvl/r
                if corr<=0.0:
                    corr = 0.001
                coh[i,j] = corr

        return coh

    def min_span_tree(self, coh, n = 1):
        self.pairsDates = []

        for i in range(n):
            weightMat = 1.0/coh
            mstMat = sparse.csgraph.minimum_spanning_tree(weightMat)

            # Convert MST index matrix into date12 list
            [s_idx_list, m_idx_list] = [date_idx_array.tolist()
                                for date_idx_array in sparse.find(mstMat)[0:2]]

            for i in range(len(m_idx_list)):
                idx = sorted([m_idx_list[i], s_idx_list[i]])
                #coh[idx[0],idx[1]] = 0.001
                date12 = self.dateList[idx[0]] + '_' + self.dateList[idx[1]]
                if date12 not in self.pairsDates:
                    #print(date12)
                    self.pairsDates.append(date12)

    def delaunay(self):

        idxs = []
        dec_year = list()
        for d in self.dateList:
            yr = d[0:4]
            mo = d[4:6]
            day = d[6:8]
            dt = datetime.date.toordinal(datetime.date(int(yr), int(mo), int(day)))
            d0 = datetime.date.toordinal(datetime.date(int(yr), 1, 1))
            doy = np.asarray(dt)-d0+1
            dec_year.append(float(yr) + (doy/365.25))

        def find_neighbors(pindex, triang):
            neighbors = list()
            for simplex in triang.simplices:
                if pindex in simplex:
                    neighbors.extend([simplex[i] for i in range(len(simplex)) if simplex[i] != pindex])
                    '''
                    this is a one liner for if a simplex contains the point we`re interested in,
                    extend the neighbors list by appending all the *other* point indices in the simplex
                    '''
            #now we just have to strip out all the dulicate indices and return the neighbors list:
            return list(set(neighbors))

        x_list = dec_year
        y_list = []
        for d in self.dateList:
            b = self.baselineDict[d]
            y_list.append(b/1000)

        tri = spatial.Delaunay(np.array([[x,y] for x,y in zip(x_list, y_list)]))

        self.pairsDates = []
        for ii in range(len(self.dateList)):
            neighbor_indices = find_neighbors(ii,tri)
            for jj in neighbor_indices:
                if int(self.dateList[ii]) < int(self.dateList[jj]):
                    self.pairsDates.append(self.dateList[ii] + '_' + self.dateList[jj])

        return None

    def maxbandwidth(self,bandwidth):

        dates = self.dateList
        self.pairsDates = []

        for ii,d in enumerate(dates[0:-1]):
            for jj in np.arange(1,bandwidth+1):
                if ii+jj < len(dates):
                    self.pairsDates.append(dates[ii] + '_' + dates[ii+jj])
        return None
                    
        return None
    def sequential1(self):
        dates = self.dateList
        self.pairsDates = []
        for ii,d in enumerate(dates[0:-1]):
            self.pairsDates.append(dates[ii] + '_' + dates[ii+1])
        return None

    def sequential2(self):
        skip=2
        dates = self.dateList
        self.pairsDates = []

        for ii,d in enumerate(dates[0:-1]):
            for jj in np.arange(1,skip+1):
                if ii+jj < len(dates):
                    self.pairsDates.append(dates[ii] + '_' + dates[ii+jj])
        return None


    def sequential3(self):
        skip=3
        dates = self.dateList
        self.pairsDates = []

        for ii,d in enumerate(dates[0:-1]):
            for jj in np.arange(1,skip+1):
                if ii+jj < len(dates):
                    self.pairsDates.append(dates[ii] + '_' + dates[ii+jj])
        return None

    def sequential4(self):
        skip=4
        dates = self.dateList
        self.pairsDates = []

        for ii,d in enumerate(dates[0:-1]):
            for jj in np.arange(1,skip+1):
                if ii+jj < len(dates):
                    self.pairsDates.append(dates[ii] + '_' + dates[ii+jj])
        return None

    def sequential5(self):
        skip=5
        dates = self.dateList
        self.pairsDates = []

        for ii,d in enumerate(dates[0:-1]):
            for jj in np.arange(1,skip+1):
                if ii+jj < len(dates):
                    self.pairsDates.append(dates[ii] + '_' + dates[ii+jj])
        return None

    def small_baseline(self, timeThreshold = 100, baselineThreshold = 100):

        return None

    def single_master(self):

        master_date = self.dateList[0]
        self.pairsDates = []

        for i in range(1, len(self.dateList)):
            date12 = master_date + '-' + self.dateList[i]
            self.pairsDates.append(date12)

        return None

    def plot_network(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt

        datefmt='%Y%m%d'
        fig1 = plt.figure(1)
        ax1=fig1.add_subplot(111)
        ax1.cla()

        for pair in self.pairsDates:
            print(pair)
            d0 = pair.split("-")[0]
            d1 = pair.split("-")[1]
            t0 = datetime.datetime.strptime(d0, datefmt)
            t1 = datetime.datetime.strptime(d1, datefmt)
            b0 = self.baselineDict[d0]
            b1 = self.baselineDict[d1]
            ax1.plot([t0,t1], [b0, b1],
                 '-ko',lw=1, ms=4, alpha=0.7, mfc='r')


        plt.title('Baseline plot')
        plt.xlabel('Time')
        plt.ylabel('Perp. Baseline')
        plt.tight_layout()

        plt.savefig("Pairs.png")
