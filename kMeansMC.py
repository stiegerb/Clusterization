from sklearn.cluster import KMeans
import numpy as np
import ROOT

def predictByHand(point, centers):
    minDist = 10000
    cluster = -1
    counter = 0
    for center in centers:
        dist = (point[0] - center[0])**2 + (point[1] - center[1])**2
        if dist < minDist:
            minDist = dist
            cluster = counter
        counter = counter +1

    return cluster

def classicalBinning(point):
    if (-1. < point[0] <= 1)   and (-1.0 < point[1] <= -0.2) : return 0
    if (-1. < point[0] <= 1)   and (-0.2 < point[1] <=  0.1) : return 1
    if (-1. < point[0] <= 0.3) and (0.1  < point[1] <=  0.4) : return 2
    if (0.3 < point[0] <= 1.)  and (0.1  < point[1] <=  0.4) : return 3
    if (-1. < point[0] <= 0.1) and (0.4  < point[1] <=  1.0) : return 4
    if (0.1 < point[0] <= 0.4) and (0.4  < point[1] <=  1.0) : return 5
    if (0.4 < point[0] <= 1.)  and (0.4  < point[1] <=  1.0) : return 6
    else: print 'one bin is missing', point[0], point[1]

Data  = []
TTbar = []
TTV   = []
TTH   = []

ttbar = open('data/ttbar.txt','r')
ttv   = open('data/ttv.txt' ,'r')
tth   = open('data/tth.txt' ,'r')

counter = 0
for line in ttbar:
    x = float( line.split(' ')[0] )
    y = float( line.split(' ')[5] )
    if counter%2 == 0: Data.append([x,y])
    else: TTbar.append([x,y])
    counter = counter +1
    if counter > 200: break

counter = 0
for line in ttv:
    x = float( line.split(' ')[0] )
    y = float( line.split(' ')[5] )
    if counter%2 == 0: Data.append([x,y])
    else: TTV.append([x,y])
    counter = counter +1
    if counter > 200: break

counter = 0
for line in tth:
    x = float( line.split(' ')[0] )
    y = float( line.split(' ')[5] )
    if counter%2 == 0: Data.append([x,y])
    else: TTH.append([x,y])
    counter = counter +1
    if counter > 200: break


X = np.array(Data)
n_clusters = 4
kmeans = KMeans(n_clusters=n_clusters, random_state=0)
kmeans.fit(X)


Hdata  = ROOT.TH1F('data' , '', n_clusters, -0.5, n_clusters-0.5)
Httbar = ROOT.TH1F('ttbar', '', n_clusters, -0.5, n_clusters-0.5)
HttV  = ROOT.TH1F('ttV' , '', n_clusters, -0.5, n_clusters-0.5)
HttH  = ROOT.TH1F('ttH' , '', n_clusters, -0.5, n_clusters-0.5)

Hdataold  = ROOT.TH1F('dataold' , '', n_clusters, -0.5, n_clusters-0.5)
Httbarold = ROOT.TH1F('ttbarold', '', n_clusters, -0.5, n_clusters-0.5)
HttVold   = ROOT.TH1F('ttVold'  , '', n_clusters, -0.5, n_clusters-0.5)
HttHold   = ROOT.TH1F('ttHold'  , '', n_clusters, -0.5, n_clusters-0.5)

for datapoint in Data:
    Hdata.Fill(predictByHand(datapoint, kmeans.cluster_centers_))
    Hdataold.Fill(classicalBinning(datapoint))
for datapoint in TTbar:
    Httbar.Fill(predictByHand(datapoint, kmeans.cluster_centers_))
    Httbarold.Fill(classicalBinning(datapoint))
for datapoint in TTV:
    HttV.Fill(predictByHand(datapoint  , kmeans.cluster_centers_))
    HttVold.Fill(classicalBinning(datapoint))
for datapoint in TTH:
    HttH.Fill(predictByHand(datapoint  , kmeans.cluster_centers_))
    HttHold.Fill(classicalBinning(datapoint))

Httbar.SetLineColor(ROOT.kRed)    ;  Httbarold.SetLineColor(ROOT.kRed)
HttV  .SetLineColor(ROOT.kBlue)   ;  HttVold.SetLineColor(ROOT.kBlue)
HttH  .SetLineColor(ROOT.kMagenta);  HttHold.SetLineColor(ROOT.kMagenta)
Httbar.SetFillColor(ROOT.kRed)    ;  Httbarold.SetFillColor(ROOT.kRed)
HttV  .SetFillColor(ROOT.kBlue)   ;  HttVold.SetFillColor(ROOT.kBlue)
HttH  .SetFillColor(ROOT.kMagenta);  HttHold.SetFillColor(ROOT.kMagenta)
Hdata.SetMarkerStyle(ROOT.kFullCircle)

mc = ROOT.THStack('mc','mc')
mcold = ROOT.THStack('mcold','mcold')
mc.Add(Httbar); mc.Add(HttV); mc.Add(HttH)
mcold.Add(Httbarold); mcold.Add(HttVold); mcold.Add(HttHold)

c = ROOT.TCanvas()
mc.Draw("HIST")
Hdata.Draw("P,SAME")
c.SaveAs('stack.pdf')
mcold.Draw("HIST")
Hdataold.Draw("P,SAME")
c.SaveAs('stackold.pdf')

StoBkg = []
StoBkgOld = []
for i in range(1,HttH.GetNbinsX()+1):
    StoBkg   .append( HttH.GetBinContent(i)    / ( Httbar.GetBinContent(i)    + HttV.GetBinContent(i)))
    StoBkgOld.append( HttHold.GetBinContent(i) / ( Httbarold.GetBinContent(i) + HttVold.GetBinContent(i)))

StoBkg.sort()
StoBkgOld.sort()

for i in range(0,HttH.GetNbinsX()):
    print 'cluster', i, StoBkg[i], StoBkgOld[i]
