from sklearn.cluster import KMeans
import numpy as np
import ROOT
import math
from scipy.stats import poisson
import matplotlib.pyplot as plt
from matplotlib import colors
from random import shuffle

class Cluster():
    def __init__( self, center, eventsTTH, eventsBkg):
        print 'initialising cluster with', len(eventsTTH), len(eventsBkg), 'events'
        self.center    = center
        self.eventsTTH    = eventsTTH
        self.eventsBkg    = eventsBkg
        self.subclusters = []
        self.pvalue = poisson.pmf( len(self.eventsTTH) + len(self.eventsTTH), len(self.eventsBkg))

    def recluster(self):
        global maxIndex
        print 'reclustering cluster with center', self.center,
        if hasattr(self, 'isClusterizable'): return
        kmeans = KMeans(n_clusters=2, random_state=0)
        kmeans.fit(np.concatenate( (self.eventsTTH, self.eventsBkg), axis =0))
        label = 0
        for center in kmeans.cluster_centers_:
            eventsTTH = np.array([point for point in self.eventsTTH if kmeans.predict(np.array([[point[0],point[1]]]))[0] == label])
            eventsBkg = np.array([point for point in self.eventsBkg if kmeans.predict(np.array([[point[0],point[1]]]))[0] == label])
            self.subclusters.append( Cluster(center, eventsTTH, eventsBkg))
            label = label+1

        if (self.pvalue < self.subclusters[0].pvalue*self.subclusters[1].pvalue):
            print '...new clustering doesnt add anything :('
            self.isClusterizable = False
            self.subclusters = []
        elif (len(self.subclusters[0].eventsTTH) < 10 or len(self.subclusters[1].eventsTTH) < 10):
            print '...too few signal events in sr'
            self.isClusterizable = False
            self.subclusters = []
        else:
            print '...is clusterizable'
            self.isClusterizable = True
            self.subclusters[0].recluster()
            self.subclusters[1].recluster()

        if not self.isClusterizable:
            self.index = maxIndex
            maxIndex = maxIndex + 1
        return

    def findSubCluster(self, point):
        if not hasattr(self, 'isClusterizable'):
            self.recluster()

        if not self.isClusterizable: return self.index

        minDist = 10000
        cluster_counter = -1
        counter = 0
        for cluster in self.subclusters:
            center = cluster.center
            dist = (point[0] - center[0])**2 + (point[1] - center[1])**2
            if dist < minDist:
                minDist = dist
                cluster_counter = counter
            counter = counter +1
        return self.subclusters[cluster_counter].findSubCluster(point)

class recursiveClustering():
    def __init__(self):
        self.getDataFromFile()
        self.clusterize()

    def getDataFromFile(self):

        TTbar = []
        TTV   = []
        TTH   = []
        TTbar2 = []
        TTV2  = []
        TTH2   = []

        ttbar = open('data/ttbar.txt','r')
        ttv   = open('data/ttv.txt' ,'r')
        tth   = open('data/tth.txt' ,'r')

        counter = 0
        for line in ttbar:
            x = float( line.split(' ')[0] )
            y = float( line.split(' ')[5] )
            if counter%2: TTbar.append([x,y])
            else: TTbar2.append([x,y])
            counter = counter +1
            if counter > 200: break

        counter = 0
        for line in ttv:
            x = float( line.split(' ')[0] )
            y = float( line.split(' ')[5] )
            if counter%2: TTV.append([x,y])
            else:  TTV2.append([x,y])
            counter = counter +1
            if counter > 200: break

        counter = 0
        for line in tth:
            x = float( line.split(' ')[0] )
            y = float( line.split(' ')[5] )
            if counter%2: TTH.append([x,y])
            else: TTH2.append([x,y])
            counter = counter +1
            if counter > 200: break

        self.eventsBkg = np.array( TTV + TTbar)
        self.eventsTTH = np.array( TTH )
        self.eventsTTV = np.array( TTV )
        self.eventsTTbar = np.array( TTbar )

        self.eventsTTH2 = np.array( TTH2 )
        self.eventsTTV2 = np.array( TTV2 )
        self.eventsTTbar2 = np.array( TTbar2 )

    def clusterize(self):
        global maxIndex
        print 'clusterizing the phase space'
        self.mainCluster = Cluster( np.array([-999,-999]), self.eventsTTH, self.eventsBkg)
        self.mainCluster.recluster()
        print 'obtained', maxIndex, 'clusters'

    def makeHistos(self):
        global maxIndex
        Httbar = ROOT.TH1F('ttbar', '', maxIndex, -0.5, maxIndex-0.5)
        HttV  = ROOT.TH1F('ttV' , '', maxIndex, -0.5, maxIndex-0.5)
        HttH  = ROOT.TH1F('ttH' , '', maxIndex, -0.5, maxIndex-0.5)
        for point in self.eventsTTH2:
            HttH.Fill( self.mainCluster.findSubCluster(point) )
        for point in self.eventsTTbar2:
            Httbar.Fill( self.mainCluster.findSubCluster(point) )
        for point in self.eventsTTV2:
            HttV.Fill( self.mainCluster.findSubCluster(point) )
        mc = ROOT.THStack('mc','mc')
        Httbar.SetLineColor(ROOT.kRed)    ;
        HttV  .SetLineColor(ROOT.kBlue)   ;
        HttH  .SetLineColor(ROOT.kMagenta);
        Httbar.SetFillColor(ROOT.kRed)    ;
        HttV  .SetFillColor(ROOT.kBlue)   ;
        HttH  .SetFillColor(ROOT.kMagenta);
        mc.Add( HttH )
        mc.Add( Httbar )
        mc.Add( HttV )
        c = ROOT.TCanvas()
        mc.Draw('hist')
        c.SaveAs('recursiveClustering.pdf')
    def makeScatterPlot(self):
        print 'making scatter plot to see how the clustering went'
        global maxIndex
        clustersX = []
        clustersY = []
        for i in range(maxIndex):
            clustersX.append([])
            clustersY.append([])

        for point in self.eventsBkg:
            subCluster = self.mainCluster.findSubCluster(point)
            clustersX[subCluster].append(point[0])
            clustersY[subCluster].append(point[1])
        for point in self.eventsTTH:
            subCluster = self.mainCluster.findSubCluster(point)
            clustersX[subCluster].append(point[0])
            clustersY[subCluster].append(point[1])

        mycolors = []
        for key in colors.cnames:
            mycolors.append( colors.cnames[key] )
        shuffle(mycolors)
        for i in range(maxIndex):
            plt.scatter(clustersX[i], clustersY[i], c = mycolors[i])
        plt.show()

if __name__ == "__main__":
    global maxIndex
    maxIndex = 0
    a = recursiveClustering()
    a.makeHistos()
#    a.makeScatterPlot()
