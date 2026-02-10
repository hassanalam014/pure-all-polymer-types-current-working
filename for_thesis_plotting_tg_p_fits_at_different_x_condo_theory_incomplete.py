# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
# from loadSpecificHeatExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
from isListOrNpyArray import *
from Parameters_of_Different_Polymers import *
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
from optimizeResidualFunctions import pureEOSResidual,pureChemicalPotentialResidual
from loadSpecificHeatExperimentalData import *
from sympy import Symbol, nsolve
import sympy
import mpmath
import winsound
# import pyttsx		#Text to speech
# from lazyme.string import palette, highlighter, formatter, color_print		#For command text color and underline etc.
# color_print('abc', color='red', underline=True, bold=True, highlight='white', faint=True, reverse=True)

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def glassTransition_for_Rratio(Rratio,P,T,R,M,S_Tg,epsilon_2,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	S_condo=-S_Tg+(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))
	res=S_condo

	return res

def Rrat(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	for i in range(0,100,1):
		
		Rratio=0.0
		try:
			Rratio = bisect(glassTransition_for_Rratio,i,i+1,args=(P,T,R,M,S_Tg,epsilon_2,Vratio,Pstar,Tstar,Rstar))
		except:
			# print("Failure to get Rratio")
			pass
		if Rratio!=0.0:
			print 'Hurry! Rratio_dependent is:', Rratio, 'epsilon_2_independent is:', epsilon_2
			break
	
	if Rratio==0.0:
		print 'Program Failed to get value of Rratio at epsilon_2=',epsilon_2
		# Rratio=50000

	return Rratio

def glassTransitionCriteria(T,P,M,S_Tg,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	# Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	# Pstarstar=Pratio*Pstar
	# Rstarstar=Rratio*Rstar
	epsilon_2=kB*Tstarstar

	S_condo=-S_Tg+(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T)))))-((((Rratio)*exp(-Tstarstar/(T)))/(1+(Rratio)*exp(-Tstarstar/(T))))*Tstarstar/(T))))
	res=S_condo
	
	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,S_Tg,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg

def ResidualArray(params,P,Tg):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	M = params['M'].value
	epsilon_2 = params['epsilon_2'].value
	Vratio = params['Vratio'].value
	x = params['x'].value
	Tg_atm = params['Tg_atm'].value

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	Rg_atm=density(P_atm,Tg_atm,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Rratio=Rrat(P_atm,Tg_atm,Rg_atm,M,S_Tg=S_Tg,epsilon_2=epsilon_2,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio,'x':S_Tg,'Tg_atm':Tg_atm}
	
	print 'epsilon_2=', epsilon_2
	print 'Rratio	=', Rratio
	print 'x	=', x
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],**kwargs)
		residual[j] = abs((Tg[j]-Tg_calculated))
	
	return residual

# engine = pyttsx.init()	## Text to Speech
# engine.say('The program has started')
# engine.runAndWait()


# Program_Running_For=['PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PS Quach 36kilo','PVAc Sandberg 189kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo','LPP Passaglia 15kilo','LPP Hollander 15kilo','BPP Passaglia 15kilo','BPP Hollander 15kilo']
# Program_Running_For=['PMMA Grassia 140kilo','PMMA Olabisi 44kilo','PS Quach 36kilo','PVAc Roland 189kilo','PVME Casalini 60kilo','PC Zoller 00kilo']

Polymer_Type='PMMA'
Reference='Olabisi'
Polymer_Weight='44kilo'
# class Polymer_Type

kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

Abelow,Bbelow,Aabove,Babove,A,B,deltaCp,T0_excluding_Tg,M0_excluding_Tg,C0_excluding_Tg,P0_excluding_Tg,I0_excluding_Tg,Tg0_excluding_Tg,T0_above_Tg,M0_above_Tg,C0_above_Tg,P0_above_Tg,I0_above_Tg,Tg0_above_Tg,T0_at_Tg,M0_at_Tg,C0_at_Tg,P0_at_Tg,I0_at_Tg,Tg0_at_Tg,T0_below_Tg,M0_below_Tg,C0_below_Tg,P0_below_Tg,I0_below_Tg,Tg0_below_Tg,T0_complete_Tg,M0_complete_Tg,C0_complete_Tg,P0_complete_Tg,I0_complete_Tg,Tg0_complete_Tg=loadSpecificHeatExperimentalData(**kwargs)
Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

P = P_atm
T=Tg_atm
M=M_infinity
# R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Rg_atm=density(P_atm,Tg_atm,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
Vratio=1.0
r = (Pstar*M)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

# Ptilde=P/Pstar
# Ttilde=T/Tstar
# Rtilde=R/Rstar
# dPtilde_dT=dP_dT_atm/Pstar
# dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar

######################################################################################################

######################################################################################################
#Ideal Experimental Straight Line Data
P_line = npy.linspace(0.101325,P_upper,20)
T_line = npy.zeros(len(P_line))
R_line=npy.zeros(len(P_line))

for i in range(0,len(P_line)):
	T_line[i]=((P_line[i]-P)/dP_dT_atm)+T
	#R_line[i]=density(P_line[i],T_line[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

# x_min_before=0.30
# x_min_after=0.40

# S_Tg= npy.linspace(x_min_before,x_min_after,10)
S_Tg=[0.3,0.6,0.8]

epsilon_2=npy.zeros(len(S_Tg))
Tratio=npy.zeros(len(S_Tg))
Rratio=npy.zeros(len(S_Tg))
#######################################################################################################

#######################################################################################################
#Fitting on Straight Line Linear Data for Own_Criteria_1:
for i in range(0,len(S_Tg)):
	print '---------------------------------'
	print 'Program is iterating for S_Tg= ', S_Tg[i]
	print '---------------------------------'
	#Fitting Data to the base curve below glass transition:
	params = Parameters()
	#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
	#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
	#						(Name,			Value,		        Vary?,	Min,	Max,	Expr)
	params.add_many((		'S_Tg',			S_Tg[i],			False,	0.0,	1.0,	None),
					(		'Vratio',		Vratio,				False,	0.0,	None,	None),
					(		'epsilon_2',	7000.0,				True,	3000.0,	None,	None),
					(		'Pstar',		Pstar,				False,	0.0,	None,	None),
					(		'Tstar',		Tstar,				False,	0.0,	None,	None),
					(		'Rstar',		Rstar,				False,	0.0,	None,	None),
					(		'M',			M,					False,	0.0,	None,	None),
					(		'Tg_atm',		Tg_atm,		        False,	0.0,	None,	None))

	#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
	#These need to be added to the experimental datapoints to find the fitted specific heat.
	fit = minimize(ResidualArray,params,args=(P_line,T_line))

	#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
	report_fit(fit.params)

	if 'epsilon_2' in fit.params:
		epsilon_2[i] = fit.params['epsilon_2'].value
		S_Tg[i] = fit.params['x'].value
		Rratio[i]=Rrat(P_atm,Tg_atm,Rg_atm,M,S_Tg=S_Tg[i],epsilon_2=epsilon_2[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
		Tratio[i]=(epsilon_2[i]/(kB*Tstar))
		# Vratio = fit.params['Vratio'].value
		#kwargs = {'A':A,'B':B}

duration = 2000  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)

#######################################################################################################

# Tstarstar_min=epsilon_2_min/kB
# Tratio_min=Tstarstar_min/Tstar
Vratio=1.0

#Condo Theory
S_Tg=[0.3,0.6,0.8]
Rratio=[2.1557238,5.73362622,63.68224322]
epsilon_2=[9182.03041965,5858.84806115,8063.93644207]

Tratio=npy.zeros(len(S_Tg))

for i in range(0,len(S_Tg)):
	Tratio[i]=epsilon_2[i]/(kB*Tstar)

#Initializing the array of densities.
P = npy.linspace(0.101325,P_upper,20)
R=npy.zeros(len(P))
Tg_calculated1=npy.zeros(len(P))
Tg_calculated2=npy.zeros(len(P))
Tg_calculated3=npy.zeros(len(P))

for i in range(0,len(P)):

	Tg_calculated1[i]=glassTemp(P[i],M=M,S_Tg=S_Tg[0],Rratio=Rratio[0],Tratio=Tratio[0],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Tg_calculated2[i]=glassTemp(P[i],M=M,S_Tg=S_Tg[1],Rratio=Rratio[1],Tratio=Tratio[1],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Tg_calculated3[i]=glassTemp(P[i],M=M,S_Tg=S_Tg[2],Rratio=Rratio[2],Tratio=Tratio[2],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	# print P[i]

# for i in range(0,len(P0)):
# 	R0[i]=density(P0[i],T0[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)	

print S_Tg
print Rratio
print epsilon_2

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

###################################################################################
###    Tg VERSUS P Plots  #### 
###################################################################################
plt.plot(P,Tg_calculated1,'k',color='r',lw=linewidth,ls='-',label='Present Theory for S_Tg={}'.format(S_Tg[0]))
plt.plot(P,Tg_calculated2,'k',color='b',lw=linewidth,ls='-',label='S_Tg={}'.format(S_Tg[1]))
plt.plot(P,Tg_calculated3,'k',color='g',lw=linewidth,ls='-',label='S_Tg={}'.format(S_Tg[2]))
# plt.plot(P_line,T_line,'k',color='r',lw=linewidth,ls='-',label='Pure PMMA Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Experiment')
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)
#############################################################################################################

#############################################################################################################
# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)
fig.savefig('./'+output_folder+r'\PMMA_Olabisi_Tg-P Straight Line Fit'+img_extension,dpi=240)

plt.show()
