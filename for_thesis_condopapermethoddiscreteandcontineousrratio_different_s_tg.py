# Date: 2019
#
# Description: The purpose of this file is to plot Polystyrene (PS) Thermodynamics Properties
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
import winsound    # Play a beep sound
#from matplotlib.ticker import AutoMinorLocator
# from all_p_params import *
# from loadSpecificHeatExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from findVectors import findVectors
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
# from wrapperFlexibilityFunctions import calculateSpecificHeat
from isListOrNpyArray import *
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
from optimizeResidualFunctions import pureEOSResidual,pureChemicalPotentialResidual
# from loadSpecificHeatExperimentalData import *
from sympy import Symbol, nsolve
import sympy
import mpmath
from Parameters_of_Different_Polymers import *


def TgCriteriaForRratio(Rratio,P,T,M,Tratio,S_Tg,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Tstarstar=Tratio*Tstar

	#Condo Theory:
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))

	res=S_condo

	return res

def Rrat(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Rratio = bisect(TgCriteriaForRratio,0.4,100,args=(P,T,M,Tratio,S_Tg,Pstar,Tstar,Rstar))
	return Rratio


def ResidualArrayContineous(params,P,T):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	M = params['M'].value
	epsilon_2 = params['epsilon_2'].value
	S_Tg = params['S_Tg'].value
	Tg_atm = params['Tg_atm'].value
	P_atm = params['P_atm'].value

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	Rratio=Rrat(P_atm,Tg_atm,M,Tratio=Tratio,S_Tg=S_Tg,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	print epsilon_2
	print Rratio
	# kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'S_Tg':S_Tg}
	# Rratio=Rrat(P_atm,Tg_atm,M,**kwargs)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Rratio':Rratio,'S_Tg':S_Tg}
	
	# print Rratio
	# print epsilon_2
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],M,**kwargs)
		residual[j] = abs((T[j]-Tg_calculated))/T[j]
	
	return residual


def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R


def TgCriteriaForEpsilon(epsilon_2,P,T,M,Rratio,S_Tg,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	#Condo Theory:
	S_condo=-S_Tg+(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))

	res=S_condo

	return res

def eps(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	# epsilon_2 = bisect(TgCriteriaForEpsilon,1500,20000,args=(P,T,M,Rratio,S_Tg,Pstar,Tstar,Rstar))
		
	for i in range(0,50000,500):
		
		epsilon_2=0.0
		try:
			epsilon_2 = bisect(TgCriteriaForEpsilon,i,i+500,args=(P,T,M,Rratio,S_Tg,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get epsilon_2")
			pass
		if epsilon_2!=0.0:
			print 'Hurry! epsilon_2_dependent is:', epsilon_2, 'Rratio_independent is:', Rratio
			break
	
	if epsilon_2==0.0:
		print 'Program Failed to get value of epsilon_2 against Rratio_independent=', Rratio
		# epsilon_2=5000

	return epsilon_2


def glassTransitionCriteria(T,P,M,Rratio,Tratio,S_Tg,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Tstarstar=Tratio*Tstar

	#Condo Theory:
	S_condo=-S_Tg+(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))

	res=S_condo

	return res

def glassTemp(P,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,Rratio,Tratio,S_Tg,Pstar,Tstar,Rstar))
	
	return Tg

def ResidualArray(P,T,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	# kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'S_Tg':S_Tg}
	# Rratio=Rrat(P_atm,Tg_atm,M,**kwargs)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'Rratio':Rratio,'S_Tg':S_Tg}
	
	# print Rratio
	# print epsilon_2
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],M,**kwargs)
		residual[j] = abs((T[j]-Tg_calculated))/T[j]
	
	return residual


Program_Running_For=['PMMA Olabisi 44kilo']

Pick_List_Element = Program_Running_For[0]
Divide_List_Picked_Element = Pick_List_Element.split()

print(Divide_List_Picked_Element)

Polymer_Type=Divide_List_Picked_Element[0]
Reference=Divide_List_Picked_Element[1]
Polymer_Weight=Divide_List_Picked_Element[2]
# class Polymer_Type

kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

M=M_infinity
dP_dT=1/dTg_dP_atm

P_line = npy.linspace(0.101325,P_upper,10)
Tg_line = npy.zeros(len(P_line))

#Ideal Experimental Straight Line Data for Fitting:
for i in range(0,len(P_line)):
	Tg_line[i]=((P_line[i]-P_atm)/dP_dT)+Tg_atm				#Condo Ref[53] Slope Line

'''
S_Tg=[-0.2,-0.1,0.0,0.1,0.2]
Rratio_GD=npy.zeros(len(S_Tg))
epsilon_2_GD=npy.zeros(len(S_Tg))

for k in range(0,len(S_Tg)):

	Rratio=[2.0,3.0,4.0,5.0]
	epsilon_2_list = npy.zeros(len(Rratio))
	residual_list=npy.zeros(len(Rratio))

	for j in range(0,len(Rratio)):
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Rratio':Rratio[j],'S_Tg':S_Tg[k]}
		epsilon_2_list[j]=eps(P_atm,Tg_atm,M,**kwargs)

	for j in range(0,len(Rratio)):
		kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Rratio':Rratio[j],'S_Tg':S_Tg[k],'epsilon_2':epsilon_2_list[j],'Tg_atm':Tg_atm,'P_atm':P_atm}
		residual_array=ResidualArray(P_line,Tg_line,**kwargs)
		SumSquareError=0

		for i in residual_array:
			SumSquareError+=i*i

		residual_list[j]=SumSquareError

	print Rratio
	print epsilon_2_list
	print residual_list

	residual_min=min(residual_list)
	index_min=npy.argmin(residual_list)
	epsilon_2_min=epsilon_2_list[index_min]
	Rratio_min=Rratio[index_min]

	print Rratio_min
	print epsilon_2_min

	Rratio_GD[k]=Rratio_min
	epsilon_2_GD[k]=epsilon_2_min

print S_Tg
print Rratio_GD
print epsilon_2_GD

########################################################################


########################################################################
#For Continuous Fitting:
#Fitting Idealized Experimental Straight Line Data:
params = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#						(Name,			Value,		        Vary?,	Min,		Max,	Expr)
params.add_many((		'epsilon_2',	7000.0,		    	True,	3000.0,		20000,	None),				
				(		'S_Tg',		S_Tg,				False,	0.0,		1.0,	None),				
				(		'M',			M,					False,	0.0,		None,	None),				
				(		'P_atm',		P_atm,				False,	0.0,		None,	None),
				(		'Tg_atm',		Tg_atm,				False,	0.0,		None,	None),
				(		'Pstar',		Pstar,				False,	0.0,		None,	None),
				(		'Tstar',		Tstar,				False,	0.0,		None,	None),
				(		'Rstar',		Rstar,				False,	0.0,		None,	None))

#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
#These need to be added to the experimental datapints to find the fitted specific heat.
fit = minimize(ResidualArrayContineous,params,args=(P_line,Tg_line))

#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
report_fit(fit.params)

if 'epsilon_2' in fit.params:
	epsilon_2 = fit.params['epsilon_2'].value

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio,'S_Tg':S_Tg}
Rratio=Rrat(P_atm,Tg_atm,M,**kwargs)

print 'Rratio is', Rratio
print 'epsilon_2 is', epsilon_2

#Play a Beep Sound
duration = 1000  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)

##########################################################################################
'''
##########################################################################################

S_Tg=[-0.2, -0.1, 0.0, 0.1, 0.2]
Rratio_GD=[2.0, 2.0, 3.0, 3.0, 4.0]
epsilon_2_GD=[8807.33003575, 6971.27016426, 7428.00611103, 6084.32031276, 6388.57482328]

Tratio_GD=npy.zeros(len(S_Tg))
for i in range(0,len(S_Tg)):
	Tratio_GD[i]=epsilon_2_GD[i]/(kB*Tstar)

#Initializing the array of densities.
P=npy.linspace(0.101325,P_upper,20)
R=npy.zeros(len(P))
Tg_calculated1=npy.zeros(len(P))
Tg_calculated2=npy.zeros(len(P))
Tg_calculated3=npy.zeros(len(P))
Tg_calculated4=npy.zeros(len(P))
Tg_calculated5=npy.zeros(len(P))

for i in range(0,len(P)):
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio_GD[0],'Rratio':Rratio_GD[0],'S_Tg':S_Tg[0]}
	Tg_calculated1[i] = glassTemp(P[i],M,**kwargs)
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio_GD[1],'Rratio':Rratio_GD[1],'S_Tg':S_Tg[1]}
	Tg_calculated2[i] = glassTemp(P[i],M,**kwargs)
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio_GD[2],'Rratio':Rratio_GD[2],'S_Tg':S_Tg[2]}
	Tg_calculated3[i] = glassTemp(P[i],M,**kwargs)
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio_GD[3],'Rratio':Rratio_GD[3],'S_Tg':S_Tg[3]}
	Tg_calculated4[i] = glassTemp(P[i],M,**kwargs)
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Tratio':Tratio_GD[4],'Rratio':Rratio_GD[4],'S_Tg':S_Tg[4]}
	Tg_calculated5[i] = glassTemp(P[i],M,**kwargs)
	print P[i]

#Setting font size
axis_size = 20
title_size = 20
size = 14
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'Plot_Thesis'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#General line properties.
linewidth = 1
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
fig=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

plt.plot(P,Tg_calculated1,'k',color='r',lw=linewidth,ls='-',label='Condo et al. for S(Tg)={} J/g.K'.format(S_Tg[0]))
# plt.plot(P,Tg_calculated2,'k',color='b',lw=linewidth,ls='-',label='S(Tg)={} J/g.K'.format(S_Tg[1]))
plt.plot(P,Tg_calculated3,'k',color='g',lw=linewidth,ls='-',label='S(Tg)={} J/g.K'.format(S_Tg[2]))
# plt.plot(P,Tg_calculated4,'k',color='k',lw=linewidth,ls='-',label='S(Tg)={} J/g.K'.format(S_Tg[3]))
plt.plot(P,Tg_calculated5,'k',color='m',lw=linewidth,ls='-',label='S(Tg)={} J/g.K'.format(S_Tg[4]))

# plt.plot(P_line,Tg_line,'k',color='r',lw=linewidth,ls='-',label='Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Experiment')

# plt.axhline(y=243.5,lw=0.5,color='k', linestyle='-.')
# plt.axvline(x=0.101325,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

fig.savefig('./'+output_folder+r'\PMMA_Olabisi_44kilo_Tg-P Condo Theory GD Challenging'+img_extension,dpi=240)

plt.show()
