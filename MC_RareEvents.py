# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 16:22:39 2016
@author: Tarek Frahi
"""
from openturns import *
from math import *
from numpy import *
from scipy import *
import openturns as ot
from openturns.viewer import View
from time import *
#%%
#Choix de la situation
k = 3
print("\n Situation : " + str(k))
#%%
####Partie 1
#Creation de la copule
R = CorrelationMatrix(2)
R[0, 1] = 2*sin(pi/6*0.17)
copula1 = NormalCopula(R)
copula2 = IndependentCopula(2)
collection = [copula1, copula2]
copula = ComposedCopula(collection)
#Creation des distributions
norm = TruncatedDistribution(ot.Normal(30., 7.5), 0., math.inf)
tri1 = Triangular(49., 50., 51.)
tri2 = Triangular(54., 55., 56.)
gamm = Gamma(3.6239, 1/134.827)
gumb = TruncatedDistribution(Gumbel(0.00068951, 663.75), 0., math.inf)
collDist = [gamm, gumb]
weight = [0.85, 0.15]
mixt = Mixture(collDist, weight)
#Situation
if k==1 :
    marginals = [tri1, tri2, norm, gamm]
else:
    if k==2:
        marginals = [tri1, tri2, norm, gumb]
    else:
        marginals = [tri1, tri2, norm, mixt]
#Creation de la copules      
X = ot.ComposedDistribution(marginals, copula)
#%%
#Lois marginales
graph1 = X.drawMarginal1DPDF(0, 40, 60, 1000)
graph1.setLegends(["dP Zm"])
graph1.setTitle("Densité de la loi marginale de Zm")
View(graph1).show()
graph2 = X.drawMarginal1DPDF(1, 45, 65, 100)
graph2.setLegends(["dP Zv"])
graph2.setTitle("Densité de la loi marginale de Zv")
View(graph2).show()
graph3 = X.drawMarginal1DPDF(2, -10, 70, 100)
graph3.setLegends(["dP Ks"])
graph3.setTitle("Densité de la loi marginale de Ks")
View(graph3).show()
graph4 = X.drawMarginal1DPDF(3, -1000, 14000, 100)
graph4.setLegends(["dP Q"])
graph4.setTitle("Densité de la loi marginale de Q dans la situation : "+str(k))
View(graph4).show()
#%%
#Indicateurs numeriques
print("\n Tableau des moyennes")
print("\n",X.getMean())
print("\n Tableau des ecarts-types")
print("\n",X.getStandardDeviation())
print("\n Matrice de covariance")
print("\n",X.getCovariance())
print("\n Matrice de Spearman")
print("\n",X.getSpearmanCorrelation())
print("\n Matrice de tau de Kendall")
print("\n",X.getKendallTau())
#%%
####Partie 2
Zv = X.getMarginal(0)
Zm = X.getMarginal(1)
Ks = X.getMarginal(2)
Q = X.getMarginal(3)
B = 10
L = 100
Zd = 58
#Fonction f
f = NumericalMathFunction(['Zv','Zm','Ks','Q'],['Zv+(Q/(Ks*10*sqrt((Zm-Zv)/100)))^(3/5)'])
#Vecteur Y
Y = RandomVector(f,RandomVector(X))
#Echantillon de Y
N = 10000
sample = Y.getSample(N)
print("\n Situation : " + str(k))
print("\n Moyenne =", mean(sample))
print("\n Ecart-type =", std(sample))
#%%
#Noyau
kernel = KernelSmoothing()
fitDist = kernel.build(sample)
graph5 = fitDist.drawPDF()
View(graph5).show()
graph5.setLegends(["dP Y"])
graph5.setTitle("Estimation par noyau de la densité de Y dans la situation : "+str(k))
#%%
#Evenement rare
ev = Event(Y, Greater(), 58)
#%%
####Partie 3
#Monte-Carlo
montecarlo = MonteCarlo(ev)
montecarlo.setMaximumCoefficientOfVariation(0.1)
montecarlo.setMaximumOuterSampling(10000)
montecarlo.setBlockSize(1000)
Log.Show(Log.INFO)
f.enableHistory()
t0 = time()
montecarlo.run()
t1 = time()
result = montecarlo.getResult()
p = result.getProbabilityEstimate()
q = result.getCoefficientOfVariation()
v = result.getVarianceEstimate()
l = result.getConfidenceLength(0.95)
b = result.getOuterSampling()
print("\n Situation : " + str(k))
print("\n ***Methode de Monte-Carlo***")
print ("\n Probabilite estimee :", p)
print("\n Coefficient de variation :", q)
print("\n Variance :", v)
print ("\n Intervalle de confiance à 95% :", Interval(p - 0.5 * l, p + 0.5 * l))
print("\n Budget : ", b)
print("\n Temps :", 1000*(t1-t0),"ms")
#%%
graph = Graph("Monte-Carlo")
graph.add(Cloud(result.getSample(10000)))
view = View(graph)
#%%
#Echantillon
xsample = f.getHistoryInput().getSample()
ysample = f.getHistoryOutput().getSample()
#%%
#Conditionnement
def CondX(Y,X,Zd):
    dim = Y.getSize();
    ind = []
    y = array(Y)
    x = array(X)
    for i in range(dim):
        if (y[i]>Zd):
            ind.append(x[i][0:4]);
    return ind
Xcond = CondX(ysample,xsample,Zd)
#Estimation par noyau
est = kernel.build(Xcond)
#%%
#Tirage preferentiel
impsamp = ImportanceSampling(ev, est)
impsamp.setMaximumCoefficientOfVariation(0.1)
impsamp.setMaximumOuterSampling(10000)
impsamp.setBlockSize(1000)
Log.Show(Log.INFO)
f.enableHistory()
t0 = time()
impsamp.run()
t1 = time()
result = impsamp.getResult()
p = result.getProbabilityEstimate()
q = result.getCoefficientOfVariation()
l = result.getConfidenceLength(0.95)
v = result.getVarianceEstimate()
b = result.getOuterSampling()
print("\n Situation : " + str(k))
print("\n ***Methode de tirage préférentiel***")
print ("\n Probabilite estimee :", p)
print("\n Coefficient de variation :", q)
print("\n Variance :", v)
print ("\n Intervalle de confiance à 95% :", Interval(p - 0.5 * l, p + 0.5 * l))
print("\n Budget : ", b)
print("\n Temps :", 1000*(t1-t0),"ms")
#%%
####Partie 4
#Copule independante
DX = ComposedDistribution(marginals, IndependentCopula(4))
#Construction de la base de polynômes
d= 4
p = 3
enumerateFunction = LinearEnumerateFunction(d)
H = [StandardDistributionPolynomialFactory(AdaptiveStieltjesAlgorithm(DX.getMarginal(i))) for i in range(4)]
productBasis = OrthogonalProductPolynomialFactory(H,LinearEnumerateFunction(d))
m = enumerateFunction.getStrataCumulatedCardinal(p)
print("\n Cardinal de la base : ",m)
#%%
#Echantillonage de la fonction f
n = 100
insample = X.getSample(n)
outsample = f(insample)
#%%
#Methode LASSO, algo LARS, selection leave-one out
algoSelection = LeastSquaresMetaModelSelectionFactory(LAR(),CorrectedLeaveOneOut())
algometa = FunctionalChaosAlgorithm(insample,outsample,DX,FixedStrategy(productBasis,m),LeastSquaresStrategy(algoSelection))
algometa.run()
result = algometa.getResult()
metamodel = result.getMetaModel()
vect = FunctionalChaosRandomVector(result)
#%%
#Evenement rare par l'approximation modele
Y = RandomVector(metamodel,RandomVector(X))
Ev = Event(Y, Greater(), 58)
#%%
#Simulation Monte-Carlo
metamodel.enableHistory()
montecarlo = MonteCarlo(Ev)
montecarlo.setMaximumCoefficientOfVariation(0.01)
montecarlo.setMaximumOuterSampling(10000)
montecarlo.setBlockSize(1000)
Log.Show(Log.INFO)
t0 = time()
montecarlo.run()
t1 = time()
result = montecarlo.getResult()
result = montecarlo.getResult()
p = result.getProbabilityEstimate()
q = result.getCoefficientOfVariation()
v = result.getVarianceEstimate()
l = result.getConfidenceLength(0.95)
b = result.getOuterSampling()
print("\n Situation : " + str(k))
print("\n ***Methode de Monte-Carlo sur l'approximation model*** :")
print ("\n Probabilite estimee :", p)
print("\n Coefficient de variation :", q)
print("\n Variance :", v)
print("\n Intervalle de confiance à 95% :", Interval(p - 0.5 * l, p + 0.5 * l))
print("\n Budget : ", b)
print("\n Temps :", 1000*(t1-t0),"ms")
#%%
#Methode de reduction de variance
#%%
#Partie 
#Indices de Sobol
for i in range(d):
    print("\n Sobol", i, "=", vect.getSobolIndex(i))
    print("\n SobolT", i, "=", vect.getSobolTotalIndex(i))
#Moyenne et covariance
print("\n Moyenne :",vect.getMean()[0])
print("\n Variance :", vect.getCovariance()[0,0])
#%%
#Estimation de la moyenne
xsample = metamodel.getHistoryInput().getSample()
ysample = metamodel.getHistoryOutput().getSample()
Xcond = CondX(ysample, xsample, Zd)
#%%