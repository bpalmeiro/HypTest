from TestStatClass import HypTest
from ROOT import *
import numpy as np


background = [9.,18.,28.,30.,16.,34.,18.,16.,18.,8.]

signal = [0.00036, 0.01330, 0.21367, 1.35614, 3.41925, 3.42042, 1.34877, 0.21450, 0.01333, 0.00026]

measurement = [9.,18.,28.,31.,19.,37.,21.,16.,18.,8.]

npoints = 30

test = HypTest(background,signal,False)

aux = test.muScan(measurement,npoints,3)


nbin = 10
npoints = 10
Emin = 2.4057
Emax = 2.5103
sh = TH1D('signal','signal',nbin,Emin,Emax)
bh = TH1D('bb','bg',nbin,Emin,Emax)
dh = TH1D('data','data',nbin,Emin,Emax)
dh0 = TH1D('data0','data0',nbin,Emin,Emax)

background = [9.,18.,28.,30.,16.,34.,18.,16.,18.,8.]
signal = [0.00036, 0.01330, 0.21367, 1.35614, 3.41925, 3.42042, 1.34877, 0.21450, 0.01333, 0.00026]
measurement = [9.,18.,28.,31.,19.,37.,21.,16.,18.,8.]

for i in range(1,11):
    bh.SetBinContent(i,background[i-1])
    dh.SetBinContent(i,measurement[i-1])
    dh0.SetBinContent(i,background[i-1])

median = [1]
data = [1]
mulist = np.linspace(0,1,11)

for mu in mulist[1:]:
    for i in range(1,11):
        sh.SetBinContent(i,mu*signal[i-1])

    mydatasource = TLimitDataSource(sh,bh,dh);
    myconfidence = TLimit.ComputeLimit(mydatasource,50000);
    data.append(myconfidence.CLs())
    mydatasource0 = TLimitDataSource(sh,bh,dh0);
    myconfidence0 = TLimit.ComputeLimit(mydatasource0,50000);
    median.append(myconfidence0.CLs())

gd = TGraph()
gm = TGraph()

for i in range(npoints+1):
    gd.SetPoint(i,mulist[i],data[i])
    gm.SetPoint(i,mulist[i],median[i])

aux['graphs'][0].Draw()
raw_input('cacolas')
gd.Draw('same')
gm.Draw('same')
