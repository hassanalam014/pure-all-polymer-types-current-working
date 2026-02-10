# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
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
from POST_THESIS_A_and_B_of_Cp_line_without_Kier_Base import*

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

### Programs My Theory

def glassTransitionCriteria(T,P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	# Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	# Pstarstar=Pratio*Pstar
	# Rstarstar=Rratio*Rstar

	epsilon_2=Tstarstar*kB

	# MY Theory:
	# Own_Criteria_1_incorrect=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
	# Own_Criteria_1_incorrect=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x)-((x/Vratio)*(ln(1+Rratio))))

	res=Own_Criteria_1

	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg
'''
def ResidualArray(params,P,Tg):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	M = params['M'].value
	epsilon_2 = params['epsilon_2'].value
	Rratio = params['Rratio'].value
	Vratio = params['Vratio'].value
	x = params['x'].value
	Tg_atm = params['Tg_atm'].value

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio,'x':x,'Tg_atm':Tg_atm}
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],**kwargs)
		residual[j] = abs((Tg[j]-Tg_calculated))
	
	return residual
'''

Program_Running_For=['PVME Casalini 60kilo']		#'PS Self_Grassia 02kilo_POST_THESIS'

Pick_List_Element = Program_Running_For[0]
Divide_List_Picked_Element = Pick_List_Element.split()

print(Divide_List_Picked_Element)

Polymer_Type=Divide_List_Picked_Element[0]
Reference=Divide_List_Picked_Element[1]
Polymer_Weight=Divide_List_Picked_Element[2]
# class Polymer_Type

kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

Abelow,Bbelow,Aabove,Babove,A,B,deltaCp,T0_excluding_Tg,M0_excluding_Tg,C0_excluding_Tg,P0_excluding_Tg,I0_excluding_Tg,Tg0_excluding_Tg,T0_above_Tg,M0_above_Tg,C0_above_Tg,P0_above_Tg,I0_above_Tg,Tg0_above_Tg,T0_at_Tg,M0_at_Tg,C0_at_Tg,P0_at_Tg,I0_at_Tg,Tg0_at_Tg,T0_below_Tg,M0_below_Tg,C0_below_Tg,P0_below_Tg,I0_below_Tg,Tg0_below_Tg,T0_complete_Tg,M0_complete_Tg,C0_complete_Tg,P0_complete_Tg,I0_complete_Tg,Tg0_complete_Tg=loadSpecificHeatExperimentalData(**kwargs)
Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

#######################################################
#POST THESIS LINE:
Abelow,Bbelow,Aabove,Babove,deltaCp,Tg_used_for_deltaCp = load_A_and_B_of_Cp_line_without_Kier_Base(**kwargs)
#######################################################
# deltaCp=0.31
# deltaCp=(Aabove-Abelow)+(Babove-Bbelow)*Tg_atm
print 'deltaCp is:', deltaCp
P = P_atm
T=Tg_atm
M=M_infinity
R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

Ptilde=P/Pstar
Ttilde=T/Tstar
Rtilde=R/Rstar
dPtilde_dT=dP_dT_atm/Pstar
dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar
Vratio=1.0

#####################################################################################################

#####################################################################################################
#Simultaneous Equation Solver, Finding g_min by changing 'x':

x= npy.linspace(0.28,0.293,100)
epsilon_2_array=npy.zeros(len(x))
Rratio_array=npy.zeros(len(x))
mpmath.mp.dps = 15

for i in range(0,len(x)):
	
	print 'Program is iterating for the cycle number = ',i+1,' with x= ', x[i]

	Rratio = Symbol('Rratio')
	epsilon_2 = Symbol('epsilon_2')
	# Incorrect Equation below:
	# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))
	#All correct Equations below:
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x[i])-((x[i]/Vratio)*(ln(1+Rratio))))
	Cp_GlassRHS=((Pstar/(Rstar*Tstar))*(Vratio)*((epsilon_2/(kB*T))**2)*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))*(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))))-deltaCp
	# F=((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))
	# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((epsilon_2/(kB*T))*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x[i])-((x[i]/Vratio)*ln(1+Rratio)))
	# Cp_GlassRHS=((Pstar/(Rstar*Tstar))*((((((Vratio*((epsilon_2)**2))*Rratio)/((kB*T)**2))*(exp(-(((Vratio*epsilon_2))/(kB*T)))))/((1+(Rratio*(exp(-(((Vratio*epsilon_2))/(kB*T))))))**2))))-deltaCp

	answer=nsolve((Own_Criteria_1, Cp_GlassRHS), (Rratio, epsilon_2), (1.50, 8000.0),verify=True)

	Rratio_array[i]=answer[0]
	epsilon_2_array[i]=answer[1]

	print Rratio_array[i]
	print epsilon_2_array[i]

Rratio=Rratio_array
epsilon_2=epsilon_2_array

Rratio_min=min(Rratio)
index_min=npy.argmin(Rratio)
epsilon_2_min=epsilon_2[index_min]
x_min=x[index_min]
x_min_upper=x[index_min+1]
x_min_lower=x[index_min-1]

print 'Thus, the answers are:'
print kwargs
print 'Rratio_min is:',Rratio_min
print 'epsilon_2_min is:',epsilon_2_min
print 'x_min_lower is:',x_min_lower
print 'x_min is:',x_min
print 'x_min_upper is:',x_min_upper
print '--------'
print Rratio[0],Rratio_min,Rratio[75],Rratio[-1]
print epsilon_2[0],epsilon_2_min,epsilon_2[75],epsilon_2[-1]
print x[0],x_min,x[75],x[-1]

Rratio=[Rratio[0],Rratio_min,Rratio[75],Rratio[-1]]
epsilon_2=[epsilon_2[0],epsilon_2_min,epsilon_2[75],epsilon_2[-1]]
x=[x[0],x_min,x[75],x[-1]]

# print '--------'
# print Rratio[15]
# print epsilon_2[15]
# print x[15]
# print '--------'
# print Rratio[-1]
# print epsilon_2[-1]
# print x[-1]
# print '--------'
# print Rratio[75]
# print epsilon_2[75]
# print x[75]
# print '--------'

######################################################################################################

P_line = npy.linspace(0.101325,200,15)
T_line = npy.zeros(len(P_line))
R_line=npy.zeros(len(P_line))

#Ideal Experimental Straight Line Data
for i in range(0,len(P_line)):
	T_line[i]=((P_line[i]-P)/dP_dT_atm)+T
	#R_line[i]=density(P_line[i],T_line[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

###############################################################################

###############################################################################
# Rratio=[1.94, 1.66, 7.27, 46.19]
# epsilon_2=[10717, 8094, 5177, 6907]
# x=[0.24, 0.32, 0.66, 0.8]

# Rratio=[11.079429620806625,1.9446027654187243,1.6569083635911956,7.2728170505998175,46.194161858927906]
# epsilon_2=[20579.413734465055,10717.354884201131,8045.45591899061,5177.636343694266,6907.8981251276855]
# x=[0.1,0.24,0.32484848484848483,0.6642424242424243,0.8]

Tratio=npy.zeros(len(x))

for i in range(0,len(x)):
	Tratio[i]=epsilon_2[i]/(kB*Tstar)

#Initializing the array of densities.
P = npy.linspace(0.101325,300,20)
R=npy.zeros(len(P))
Tg_calculated1=npy.zeros(len(P))
Tg_calculated2=npy.zeros(len(P))
Tg_calculated3=npy.zeros(len(P))
Tg_calculated4=npy.zeros(len(P))
Tg_calculated5=npy.zeros(len(P))

for i in range(0,len(P)):

	Tg_calculated1[i]=glassTemp(P[i],M=M,x=x[0],Rratio=Rratio[0],Tratio=Tratio[0],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Tg_calculated2[i]=glassTemp(P[i],M=M,x=x[1],Rratio=Rratio[1],Tratio=Tratio[1],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Tg_calculated3[i]=glassTemp(P[i],M=M,x=x[2],Rratio=Rratio[2],Tratio=Tratio[2],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Tg_calculated4[i]=glassTemp(P[i],M=M,x=x[3],Rratio=Rratio[3],Tratio=Tratio[3],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	# Tg_calculated5[i]=glassTemp(P[i],M=M,x=x[4],Rratio=Rratio[4],Tratio=Tratio[4],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	# print P[i]

# for i in range(0,len(P0)):
# 	R0[i]=density(P0[i],T0[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)	

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
output_folder = 'Post_Thesis_Iteration_Results'

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
figPURE_POLYMER=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

plt.plot(P,Tg_calculated1,'k',color='r',lw=linewidth,ls='-',label='Present Theory for x={}'.format(round(x[0],2)))
plt.plot(P,Tg_calculated2,'k',color='y',lw=linewidth,ls='-',label='x={}'.format(round(x[1],2)))
plt.plot(P,Tg_calculated3,'k',color='b',lw=linewidth,ls='-',label='x={}'.format(round(x[2],2)))
plt.plot(P,Tg_calculated4,'k',color='m',lw=linewidth,ls='-',label='x={}'.format(round(x[3],2)))
# plt.plot(P,Tg_calculated5,'k',color='g',lw=linewidth,ls='-',label='x={}'.format(round(x[4],2)))
# plt.plot(P_line,T_line,'k',color='r',lw=linewidth,ls='-',label='Pure PMMA Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Experiment')
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
# plt.axis([0,4,370,380])
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.title(kwargs, fontdict=None, loc='center', pad=None)
# plt.subplots_adjust(bottom=0.3)

figPURE_POLYMER.savefig('./'+output_folder+r'\Tg(P)_pure_'+Polymer_Type+'_'+Reference+'_'+Polymer_Weight+img_extension,dpi=240)

plt.show()
