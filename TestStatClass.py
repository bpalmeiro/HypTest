# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 13:46:32 2015

@author: brais
"""
from __future__ import division
import ROOT as rt
from math import *
from copy import *
import numpy as np

def prob(x,lam):
    return exp(-lam)*(lam)**x/gamma(x+1.)

rd = rt.TRandom3(0)


class HypTest:

    def __init__(self,H0,H1,BFit = True):

        if not len(H1) == len(H0):
            print 'Ambas hipotesis tienen que tener la misma dimension'
            return

        self.h0 = np.array(H0)
        self.h1 = np.array(H1)
        self.dim = len(H0)
        if BFit:
            self.Calc()


    def LogLikeVal(self,data):
        '''
        Calcula el valor del LogLike para un dato dado
        '''
        if not self.CheckData(data):
            return
        else:
            n = self.dim

        mu0 = self.h0
        mu1 = self.h1
        L = 1.

        for j in range(n):
            L *= prob(data[j],mu1[j])/prob(data[j],mu0[j])
        L = 2*log(L)
        return L

    def Calc(self,xinf = 'nada',xsup = 'nada',Nbin = 100):
        '''
        Crea los histogramas de LogLike para las dos hipotesis
        '''

        try:
            self.L0list
            self.L1list
        except:
            self.LogLikeList()

        L0 = self.L0list
        L1 = self.L1list
        N = len(self.L0list)
        rt.gStyle.SetOptStat(0)
        if xinf == 'nada': xinf = min(min(L0),min(L1))
        if xsup == 'nada': xsup = max(max(L0),max(L1))
        hist0 = rt.TH1D("","",Nbin,xinf,xsup)
        hist1 = rt.TH1D("","",Nbin,xinf,xsup)
        for i in range(N):
            hist0.Fill(L0[i])
            hist1.Fill(L1[i])
        self.hist0 = hist0
        self.hist1 = hist1


    def LogLikeList(self,Ndata = 100000):
        '''
        Crea una lista con los valores que se introducirían en los
        histogramas
        '''
        mu0 = self.h0
        mu1 = self.h1
        n = self.dim
        L0list = []
        L1list = []
        for i in range(Ndata):

            data0 = []
            data1 = []
            for j in range(n):
                data0.append(rd.Poisson(mu0[j]))
                data1.append(rd.Poisson(mu1[j]))
            L0 = self.LogLikeVal(data0)
            L1 = self.LogLikeVal(data1)

            L0list.append(L0)
            L1list.append(L1)

        self.L0list = sorted(L0list)
        self.L1list = sorted(L1list)

    def Plot(self,rebin = 0,Refresh = False):
        '''
        Grafica los histogramas de LogLike
        '''
        try:
            self.hist0
            self.hist1
        except:
            self.Calc()

        if Refresh: self.Calc()

        c = rt.TCanvas()
        h0 = deepcopy(self.hist0)
        h0.Rebin(rebin)
        h0.Draw()
        h0.SetFillColor(2)
        h0.SetFillStyle(3001)
        h0.SetLineColor(2)

        h1 = deepcopy(self.hist1)
        h1.SetFillColor(rt.kBlue)
        h1.Rebin(rebin)
        h1.SetFillStyle(3051)
        h1.Draw("same")


        #leg = rt.TLegend(.65,.70,.90,.90)
        #leg.AddEntry(h0,"H0")
        #leg.AddEntry(h1,"H1")
        #leg.Draw("same")
        return c,h0,h1#,leg


    def Integral(self,Value,Hlist,rev = False):
        '''
        Integra un histograma en forma de lista, para no depender del
        bineado, desde un punto dado hacia la direccion
        indicada por direc, si es positivo hacia la izquierda y hacia la
        derecha en caso negativo.
        '''

        if not rev:
            Llist = sorted(Hlist[:])
        else:
            Llist = sorted(Hlist[:])[::-1]

        N = len(Hlist)
        aux = True
        i = 0

        while aux:
            if not rev:
                if Value <= Llist[i] or i == N-1:
                    alpha = float(i+1)/N
                    aux = False
            else:
                if Value >= Llist[i] or i == N-1:
                    alpha = float(i+1)/N
                    aux = False
            i += 1

        return alpha

    def CheckData(self,data):
        '''
        Comprueba que el dato introducido en las funciones correspondientes
        tiene la dimensión adecuada
        '''
        ok = True

        if not len(data) == self.dim:
            ok = False
            print 'Las dimensiones del dato tienen que coincidir!'
        return ok

    def PValues(self,data,CLout = False,logopt=False):
        '''
        Calcula el alfa y el beta de las distribuciones con respecto a un
        valor dado
        '''
        if not logopt:
            if not self.CheckData(data):
                return

        try:
            self.L0list
            self.L1list
        except:
            self.LogLikeList()
        if not logopt:
            LLdata = self.LogLikeVal(data)
        else:
            LLdata = data
        List0 = self.L0list
        List1 = self.L1list
        beta = self.Integral(LLdata,List1,+1)
        alpha = self.Integral(LLdata,List0,-1)

        if CLout:
            return alpha,(1-beta)*100
        else:
            return alpha,1-beta

    def CLs(self,data,logopt=False):
        '''
        Calcula el cls de una medida
        '''
        if not logopt:
            if not self.CheckData(data):
                return

        alpha,beta = self.PValues(data,False,logopt)
        return beta/(1-alpha)

    def median(self,lista):
        '''
        Calcula la mediana de un set de datos en una lista
        '''
        Llist = sorted(lista[:])
        N = len(Llist)
        if len(Llist)%2:
            median = (Llist[N/2]+Llist[N/2+1])/2.
        else:
            median = Llist[int(ceil(N/2))]
        return median

    def findSigmaData(self,sigma,lista):
        '''
        Calcula el dato para el sigma dada
        '''

        i = int(sigma*len(lista)-1)

        return lista[i]


    def ROCcurve(self,N=100):
        '''
        Representa la curva de ROC con N puntos
        '''

        Xmin = max(min(self.L0list),min(self.L1list))
        Xmax = min(max(self.L0list),max(self.L1list))

        Graph = rt.TGraph()
        i = 0

        for data in np.linspace(Xmin,Xmax,N):
            a,b = self.PValues(data,logopt=True)

            Graph.SetPoint(i,a,1-b)
            i += 1

        Graph.Draw('alp')
        Graph.SetLineColor(3)
        Graph.GetXaxis().SetTitle("#alpha")
        Graph.GetYaxis().SetTitle("1-#beta")
        return Graph

    def muScan(self,data,npuntos,maxmu=1.):
        '''
        Hace un barrido para un CL deseado y un numero de puntos dado
        de como se comporta la señal frente al ruido a medida que se
        aumenta esta en factor mu
        '''

        print 'Starting ...'
        shift = 1./npuntos
        npuntos += 1
        mus = np.linspace(0,maxmu,npuntos,endpoint=True)
        H0backup = list(self.h0[:])
        H1backup = list(self.h1[:])

        CLsMedian = []
        CLsData = []
        CLsSigmaU = []
        CLsSigmaD = []
        CLsSigma2U = []
        CLsSigma2D = []
        for k in mus:

            self.h1 = list(np.array(H0backup)+k*np.array(H1backup)) #
            self.LogLikeList()

            median = self.median(self.L0list)

            sigmaU = self.findSigmaData( 0.8413,self.L0list)
            sigmaD = self.findSigmaData( 0.1587,self.L0list)
            sigma2U = self.findSigmaData( 0.9772,self.L0list)
            sigma2D = self.findSigmaData( 0.0228,self.L0list)

            CLsData.append(self.CLs(data))
            CLsMedian.append(self.CLs(median,logopt=True))
            CLsSigmaU.append(self.CLs(sigmaU,logopt=True))
            CLsSigmaD.append(self.CLs(sigmaD,logopt=True))
            CLsSigma2U.append(self.CLs(sigma2U,logopt=True))
            CLsSigma2D.append(self.CLs(sigma2D,logopt=True))

        self.h1 = H1backup

        GraphMedian = rt.TGraph()
        GraphData = rt.TGraph()
        GraphSigma = rt.TGraph()
        GraphSigma2 = rt.TGraph()
        for i in range(npuntos):
                GraphData.SetPoint(i,mus[i],CLsData[i])
                GraphMedian.SetPoint(i,mus[i],CLsMedian[i])
                GraphSigma.SetPoint(i,mus[i],CLsSigmaU[i])
                GraphSigma.SetPoint(npuntos+i,mus[-i-1],CLsSigmaD[-i-1])
                GraphSigma2.SetPoint(i,mus[i],CLsSigma2U[i])
                GraphSigma2.SetPoint(npuntos+i,mus[-i-1],CLsSigma2D[-i-1])

        c = rt.TCanvas()
        GraphSigma2.Draw("af")
        GraphSigma2.GetXaxis().SetTitle("#mu")
        GraphSigma2.GetYaxis().SetTitle("CLs")
        GraphSigma2.SetMinimum(min(CLsData+CLsMedian+CLsSigmaU+CLsSigmaD))
        GraphSigma2.SetFillStyle(3051)
        GraphSigma2.SetFillColor(4)

        GraphMedian.Draw("l*")
        GraphMedian.SetFillColor(0)


        GraphData.Draw("*l")
        GraphData.SetFillColor(0)
        GraphData.SetLineColor(3)

        GraphSigma.Draw("f")
        GraphSigma.SetFillStyle(3051)
        GraphSigma.SetFillColor(2)


        leg = rt.TLegend(.65,.70,.90,.90)
        leg.AddEntry(GraphMedian,"Median")
        leg.AddEntry(GraphData,"Data")
        leg.AddEntry(GraphSigma,"1 #sigma")
        leg.AddEntry(GraphSigma2,"2 #sigma")
        leg.Draw("same")

        output = {}
        output['CLsData'] = CLsData
        output['CLsMedian'] = CLsMedian
        output['CLsSigmaU'] = CLsSigmaU
        output['CLsSigmaD'] = CLsSigmaD
        output['CLsSigma2U'] = CLsSigma2U
        output['CLsSigma2D'] = CLsSigma2D
        output['graphs'] = [c,GraphMedian,GraphData,GraphSigma,GraphSigma2,leg]
        return output


if __name__ == '__main__':

    a = [9.,18.,28.,30.,16.,34.,18.,16.,18.,8.]#input('Introduzca la hipotesis 0 como una lista: ')
    b = [0.0003600000054575503, 0.013299999758601189, 0.21367000043392181, 1.3561400175094604, 3.419250011444092, 3.420419931411743, 1.348770022392273, 0.21449999511241913, 0.013330000452697277, 0.00026000000070780516]

    k = [9.,18.,28.,31.,19.,37.,21.,16.,18.,8.]#input('Introduzca la hipotesis 1 como una lista: ')
    c = HypTest(a,b,False)

    text = '''
            Metodos posibles:
                #1. Graficar los histogramas de las dos hipotesis
                #2. Graficar la curva de ROC
                #3. Calcular los P-Value de un dato
                #4. Calcular el CLs de un dato
                #5. Hacer un mu-scan
           '''
    print text
    aux = True
    while aux:
        d = input('Seleccione el metodo deseado: ')
        if d == 1:
            e = c.Plot()
        if d == 2:
            e = c.ROCcurve()
        if d == 3:
            e = input('Introduzca la medida: ')
            print c.PValues(e)
        if d == 4:
            e = input('Introduzca la medida: ')
            print c.CLs(e)
        if d == 5:
            e = k#input('Introduzca la medida: ')
            f = 50 #input('Introduzca el numero de puntos: ')
            g = c.muScan(e,f)

        print 'Desea hacer otro cálculo, si es asi escriba True, en caso contrario escriba False'
        aux = input()
