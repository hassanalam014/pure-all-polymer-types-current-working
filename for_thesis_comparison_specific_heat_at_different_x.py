# Date: April 2017
#
# Description: The purpose of this file is to plot Polystyrene (PS) density information based on experiment and theory for comparison.
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
# from all_p_params import *
# from loadSpecificHeatExperimentalData import *
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
from loadSpecificHeatExperimentalData import *
from Parameters_of_Different_Polymers import *
from matplotlib.ticker import FormatStrFormatter #To fix axis decimal places

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.0000000000000001,0.9999999999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def specificHeat_myTheory(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters and "three" flexibility parameters.
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	

	Rratio=Rstarstar/Rstar

	# Pratio=Pstarstar/Pstar
	# Tratio=Tstarstar/Tstar
	# C_myTheory_this_is_also_correct=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	
	epsilon_2=kB*Tstarstar
	Vratio=1.0

	C_myTheory=((Pstar/(Rstar*Tstar))*(Vratio)*((epsilon_2/(kB*T))**2)*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))*(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T)))))))

	return C_myTheory

def specificHeat_kier(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	C_kier=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
	
	return C_kier

def specificHeat_line_below_Tg(P,T,R,M,**kwargs):     #Cp/mass  and A, B
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	C_line_below_Tg=A+B*T
	
	return C_line_below_Tg

def specificHeat_line_above_Tg(P,T,R,M,**kwargs):     #Cp/mass  and A, B
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	C_line_above_Tg=A+B*T
	
	return C_line_above_Tg


Program_Running_For=['PMMA Olabisi 44kilo']
Pick_List_Element = Program_Running_For[0]
Divide_List_Picked_Element = Pick_List_Element.split()
print(Divide_List_Picked_Element)
Polymer_Type=Divide_List_Picked_Element[0]
Reference=Divide_List_Picked_Element[1]
Polymer_Weight=Divide_List_Picked_Element[2]
kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}
Abelow,Bbelow,Aabove,Babove,A,B,deltaCp,T0_excluding_Tg,M0_excluding_Tg,C0_excluding_Tg,P0_excluding_Tg,I0_excluding_Tg,Tg0_excluding_Tg,T0_above_Tg,M0_above_Tg,C0_above_Tg,P0_above_Tg,I0_above_Tg,Tg0_above_Tg,T0_at_Tg,M0_at_Tg,C0_at_Tg,P0_at_Tg,I0_at_Tg,Tg0_at_Tg,T0_below_Tg,M0_below_Tg,C0_below_Tg,P0_below_Tg,I0_below_Tg,Tg0_below_Tg,T0_complete_Tg,M0_complete_Tg,C0_complete_Tg,P0_complete_Tg,I0_complete_Tg,Tg0_complete_Tg=loadSpecificHeatExperimentalData(**kwargs)
Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper=Parameters_of_Different_Polymers(**kwargs)

P0 = P_atm
M0=M_infinity
r = (Pstar*M0)/(kB*Tstar*Rstar)
Vratio=1.0
##########################################################################################################

##########################################################################################################

# x=			[0.24,			0.24842105,		0.25684211,		0.26526316,		0.27368421,		0.28210526,		0.29052632,		0.29894737,		0.30736842,		0.31578947,		0.32421053,		0.33263158,		0.34105263,		0.34947368,		0.35789474,		0.36631579,		0.37473684,		0.38315789,		0.39157895,		0.4]
# Rratio=		[2.34727994,	2.28477248,		2.23297183,		2.19054442,		2.15654731,		2.13023182,		2.11066765,		2.09756025,		2.0902686,		2.08839558,		2.09165171,		2.10004339,		2.11287613,		2.1305158,		2.15249371,		2.17893098,		2.20970812,		2.24489437,		2.28449351,		2.32856757]
# epsilon_2=	[11018.26544103,10707.69623088,	10413.43213637,	10134.40522688,	9869.76716689,	9618.76944159,	9380.44041676,	9154.32356951,	8939.61178405,	8735.68433104,	8541.99727317,	8358.28210935,	8183.54789234,	8017.85457869,	7860.40305671,	7710.98299848,	7569.12672176,	7434.55012821,	7306.89133645,	7185.83798169]

# x=[0.24,0.28210526,0.31578947, 0.6, 0.8]
# Rratio=[2.34727994,2.13023182,2.08839558, 5.73362622, 63.68224322]
# epsilon_2=[11018.26544103,9618.76944159,8735.68433104, 5858.84806115, 8063.93644207]

# x=[0.3,0.6,0.8]
# Rratio=[2.1557238,5.73362622,63.68224322]
# epsilon_2=[9182.03041965,5858.84806115,8063.93644207]

#Olabisi:
Rratio=[11.079429620806625,1.9446027654187243,1.6569083635911956,7.2728170505998175,46.194161858927906]
epsilon_2=[20579.413734465055,10717.354884201131,8045.45591899061,5177.636343694266,6907.8981251276855]
x=[0.1,0.24,0.32484848484848483,0.6642424242424243,0.8]

#Grassia:
# Rratio=[1.53388852,1.23147997,1.14731937,1.73562219,5.04649589]
# epsilon_2=[9300.86567,7651.159035,6533.291519,4615.37796,4531.818697]
# x=[0.2,0.25263158,0.30526316,0.51578947,0.7]

Tstarstar=npy.zeros(len(x))		
Rstarstar=npy.zeros(len(x))		
for i in range(0,len(x)):
	Tstarstar[i]=epsilon_2[i]/kB
	Rstarstar[i]=Rratio[i]*Rstar

#Initializing the array of densities.
T_max=T0_above_Tg[-1]
T_min=T0_above_Tg[0]
T0 = npy.linspace(T_min,T_max,100)
# T0 = npy.linspace(200,600,100)

C_line_below_Tg=npy.zeros(len(T0))				#Cp = A + B*T
C_line_above_Tg=npy.zeros(len(T0))				#Cp = A + B*T
C_kier=npy.zeros(len(T0))				#Cp Kier
C_baseFit_below_Tg=npy.zeros(len(T0))			#Cp_kier + C_line  i.e. Cp Base Curve Fit
C_baseFit_above_Tg=npy.zeros(len(T0))			#Cp_kier + C_line  i.e. Cp Base Curve Fit
C1=npy.zeros(len(T0))			#Cp My Theory Only
C2=npy.zeros(len(T0))			#Cp My Theory Only
C3=npy.zeros(len(T0))			#Cp My Theory Only
C4=npy.zeros(len(T0))			#Cp My Theory Only
C5=npy.zeros(len(T0))			#Cp My Theory Only
C_total1=npy.zeros(len(T0))		#Total Cp from My Theory
C_total2=npy.zeros(len(T0))		#Total Cp from My Theory
C_total3=npy.zeros(len(T0))		#Total Cp from My Theory
C_total4=npy.zeros(len(T0))		#Total Cp from My Theory
C_total5=npy.zeros(len(T0))		#Total Cp from My Theory
R0=npy.zeros(len(T0))

for i in range(0,len(T0)):
	R0[i]=density(P0,T0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

for i in range(0,len(T0)):
	C_line_below_Tg[i] = specificHeat_line_below_Tg(P0,T0[i],R0[i],M0,A=A,B=B)
	C_line_above_Tg[i] = specificHeat_line_above_Tg(P0,T0[i],R0[i],M0,A=Aabove,B=Babove)
	C_kier[i] = specificHeat_kier(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	C_baseFit_below_Tg[i] = C_line_below_Tg[i]+C_kier[i]
	C_baseFit_above_Tg[i] = C_line_above_Tg[i]+C_kier[i]
	C1[i] = specificHeat_myTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Tstarstar=Tstarstar[0],Rstarstar=Rstarstar[0])
	C2[i] = specificHeat_myTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Tstarstar=Tstarstar[1],Rstarstar=Rstarstar[1])
	C3[i] = specificHeat_myTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Tstarstar=Tstarstar[2],Rstarstar=Rstarstar[2])
	C4[i] = specificHeat_myTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Tstarstar=Tstarstar[3],Rstarstar=Rstarstar[3])
	C5[i] = specificHeat_myTheory(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Tstarstar=Tstarstar[4],Rstarstar=Rstarstar[4])
	C_total1[i] = C1[i] + C_kier[i] + C_line_below_Tg[i]
	C_total2[i] = C2[i] + C_kier[i] + C_line_below_Tg[i]
	C_total3[i] = C3[i] + C_kier[i] + C_line_below_Tg[i]
	C_total4[i] = C4[i] + C_kier[i] + C_line_below_Tg[i]
	C_total5[i] = C5[i] + C_kier[i] + C_line_below_Tg[i]

# delete_first_n_elements=23
# T0 = T0[delete_first_n_elements:]
# C_total = C_total[delete_first_n_elements:]

# delete_first_n_experiment_data=12
# T0_complete_Tg=T0_complete_Tg[delete_first_n_experiment_data:]
# C0_complete_Tg=C0_complete_Tg[delete_first_n_experiment_data:]
# Plotting Theoretical Curve For Difference Values of Pressure:

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

# plt.axvline(x=Tg_atm,lw=0.5,color='k', linestyle='--',label='Glass Temperature')
# plt.plot(T0,C_line_below_Tg,'k',color='r',lw=linewidth,ls='-',label='C_line_below_Tg')
# plt.plot(T0,C_line_above_Tg,'k',color='r',lw=linewidth,ls='-',label='C_line_above_Tg')
# plt.plot(T0,C_kier,'k',color='b',lw=linewidth,ls='--',label='C_Kier theory')
# plt.plot(T0,C_baseFit_below_Tg,'k',color='g',lw=linewidth,ls='-.',label='Linear-Fit below Tg')
# plt.plot(T0,C_baseFit_above_Tg,'k',color='g',lw=linewidth,ls='-.',label='Linear-Fit above Tg')
# plt.plot(T0,C,'k',color='c',lw=linewidth,ls=':',label='C Only')
plt.plot(T0,C_total1,'k',color='r',lw=linewidth,ls='-',label='Present Theory for x={}'.format(round(x[0],2)))
# plt.plot(T0,C_total2,'k',color='y',lw=linewidth,ls='-',label='x={}'.format(round(x[1],2)))
plt.plot(T0,C_total3,'k',color='b',lw=linewidth,ls='-',label='x={}'.format(round(x[2],2)))
# plt.plot(T0,C_total4,'k',color='m',lw=linewidth,ls='-',label='x={}'.format(round(x[3],2)))
plt.plot(T0,C_total5,'k',color='g',lw=linewidth,ls='-',label='x={}'.format(round(x[4],2)))

plt.plot(T0_above_Tg,C0_above_Tg,'sk',ms=markersize,label='Experiment')
# plt.plot(T0_complete_Tg,C0_complete_Tg,'sk',ms=markersize,label='Experiment')

# plt.axis([325,342,1.80,2.00])
plt.xlabel('T(K)',fontsize=axis_size)
plt.ylabel(r'$\mathrm{C_P}$( $\mathrm{J/g.K}$ )',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
# plt.subplots_adjust(left=0.15,right=0.95,top=0.95,bottom=0.10,wspace=0.30,hspace=0.25)
fig.savefig('./'+output_folder+r'PMMA_Olabisi_Cp_different x Method 2_new'+img_extension,dpi=240)
plt.show()
