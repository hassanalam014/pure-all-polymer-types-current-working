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
from scipy.interpolate import spline

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def glassTransition_for_Rratio(Rratio,P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Vratio*epsilon_2))/(kB*T)))/(1+Rratio*exp(-((Vratio*epsilon_2))/(kB*T))))+((1/Vratio)*ln(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))-(x)-((x/Vratio)*ln(1+Rratio)))
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x)-((x/Vratio)*(ln(1+Rratio))))
	res=Own_Criteria_1

	return res

def Rrat(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	for i in range(0,1000,1):
		
		Rratio=0.0
		try:
			Rratio = bisect(glassTransition_for_Rratio,i,i+1,args=(P,T,R,M,x,epsilon_2,Vratio,Pstar,Tstar,Rstar))
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
	epsilon_2=kB*Tstarstar

	# MY Theory:
	# Own_Criteria_1_incorrect_perhaps=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1/Rtilde)*(1-Rtilde)*(ln(1-Rtilde)))-((1/r)*(ln(Rtilde)))+((epsilon_2/(kB*T))*((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))-((1/Vratio)*(ln(1-((Rratio*exp(-Vratio*epsilon_2/(kB*T)))/(1+(Rratio*exp(-Vratio*epsilon_2/(kB*T))))))))-(x)-((x/Vratio)*(ln(1+Rratio))))

	res=Own_Criteria_1

	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
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
	Rratio=Rrat(P_atm,Tg_atm,Rg_atm,M,x=x,epsilon_2=epsilon_2,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'M':M,'Tratio':Tratio,'Rratio':Rratio,'Vratio':Vratio,'x':x,'Tg_atm':Tg_atm}
	
	print 'epsilon_2=', epsilon_2
	print 'Rratio	=', Rratio
	print 'x	=', x
	
	residual=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],**kwargs)
		residual[j] = abs((Tg[j]-Tg_calculated))
	
	return residual

Program_Running_For=['PMMA Olabisi 44kilo']

Pick_List_Element = Program_Running_For[0]
Divide_List_Picked_Element = Pick_List_Element.split()

print(Divide_List_Picked_Element)

Polymer_Type=Divide_List_Picked_Element[0]
Reference=Divide_List_Picked_Element[1]
Polymer_Weight=Divide_List_Picked_Element[2]

kwargs = {'Polymer_Type':Polymer_Type,'Reference':Reference,'Polymer_Weight':Polymer_Weight}

# Abelow,Bbelow,Aabove,Babove,A,B,deltaCp,T0_excluding_Tg,M0_excluding_Tg,C0_excluding_Tg,P0_excluding_Tg,I0_excluding_Tg,Tg0_excluding_Tg,T0_above_Tg,M0_above_Tg,C0_above_Tg,P0_above_Tg,I0_above_Tg,Tg0_above_Tg,T0_at_Tg,M0_at_Tg,C0_at_Tg,P0_at_Tg,I0_at_Tg,Tg0_at_Tg,T0_below_Tg,M0_below_Tg,C0_below_Tg,P0_below_Tg,I0_below_Tg,Tg0_below_Tg,T0_complete_Tg,M0_complete_Tg,C0_complete_Tg,P0_complete_Tg,I0_complete_Tg,Tg0_complete_Tg=loadSpecificHeatExperimentalData(**kwargs)
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
P_line = npy.linspace(0.101325,P_upper,5)
T_line = npy.zeros(len(P_line))
R_line=npy.zeros(len(P_line))

for i in range(0,len(P_line)):
	T_line[i]=((P_line[i]-P)/dP_dT_atm)+T
	#R_line[i]=density(P_line[i],T_line[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

# x_min_before=0.30
# x_min_after=0.40
'''
# x= npy.linspace(x_min_before,x_min_after,10)
x = npy.linspace(0.24,0.40,20)

epsilon_2=npy.zeros(len(x))
Tratio=npy.zeros(len(x))
Rratio=npy.zeros(len(x))
#######################################################################################################

#######################################################################################################
#Fitting on Straight Line Linear Data for Own_Criteria_1:
for i in range(0,len(x)):
	print '---------------------------------'
	print 'Program is iterating for x= ', x[i]
	print '---------------------------------'
	#Fitting Data to the base curve below glass transition:
	params = Parameters()
	#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
	#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
	#						(Name,			Value,		        Vary?,	Min,	Max,	Expr)
	params.add_many((		'x',			x[i],				False,	0.0,	1.0,	None),
					(		'Vratio',		Vratio,				False,	0.0,	None,	None),
					(		'epsilon_2',	5000.0,				True,	0.0,	None,	None),
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
		x[i] = fit.params['x'].value
		Rratio[i]=Rrat(P_atm,Tg_atm,Rg_atm,M,x=x[i],epsilon_2=epsilon_2[i],Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
		Tratio[i]=(epsilon_2[i]/(kB*Tstar))

duration = 500  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)
'''
#######################################################################################################

# Tstarstar_min=epsilon_2_min/kB
# Tratio_min=Tstarstar_min/Tstar
# Vratio=1.0

#PMMA Olabisi For Thesis:
x=			[0.24,				0.24842105,		0.25684211,		0.26526316,		0.27368421,		0.28210526,		0.29052632,		0.29894737,		0.30736842,		0.31578947,0.32421053,0.33263158,0.34105263,0.34947368,0.35789474,0.36631579,0.37473684,0.38315789,0.39157895,0.4]
Rratio=		[2.34727994,		2.28477248,		2.23297183,		2.19054442,		2.15654731,		2.13023182,		2.11066765,		2.09756025,		2.0902686,		2.08839558,2.09165171,2.10004339,2.11287613,2.1305158,2.15249371,2.17893098,2.20970812,2.24489437,2.28449351,2.32856757]
epsilon_2=	[11018.26544103,	10707.69623088,	10413.43213637,	10134.40522688,	9869.76716689,	9618.76944159,	9380.44041676,	9154.32356951,	8939.61178405,	8735.68433104,8541.99727317,8358.28210935,8183.54789234,8017.85457869,7860.40305671,7710.98299848,7569.12672176,7434.55012821,7306.89133645,7185.83798169]

#PMMA Grassia for APS Meeting Full Data:
# Rratio = [1.53388852,1.34254729,1.23147997,1.17135786,1.14731937,1.15068658,1.17575506,1.21947618,1.28169671,1.36194033,1.46221882,1.58537485,1.73562219,1.91957776,2.14763818,2.43322031,2.80703228,3.30298643,3.99760896,5.04649589]
# x = [0.2,0.22631579,0.25263158,0.27894737,0.30526316,0.33157895,0.35789474,0.38421053,0.41052632,0.43684211,0.46315789,0.48947368,0.51578947,0.54210526,0.56842105,0.59473684,0.62105263,0.64736842,0.67368421,0.7]
# epsilon_2 = [9300.86567,8394.329994,7651.159035,7038.787026,6533.291519,6115.156293,5767.539993,5477.224129,5235.386477,5033.218525,4865.365107,4727.275844,4615.37796,4527.404744,4462.435904,4419.041613,4402.037587,4409.848973,4449.253426,4531.818697]

#PMMA Grassia for APS Meeting Reduced Data:
# Rratio = npy.array([1.53388852,1.34254729,1.23147997,1.17135786,1.14731937,1.15068658,1.17575506,1.21947618,1.28169671,1.36194033,1.46221882])
# x = npy.array([0.2,0.22631579,0.25263158,0.27894737,0.30526316,0.33157895,0.35789474,0.38421053,0.41052632,0.43684211,0.46315789])
# epsilon_2 = npy.array([9300.86567,8394.329994,7651.159035,7038.787026,6533.291519,6115.156293,5767.539993,5477.224129,5235.386477,5033.218525,4865.365107])

# # TO smoothen plot for APS Meeting: 300 represents number of points to make between x.min and x.max
# xnew = npy.linspace(x.min(), x.max(), 300)  
# Rratio_smooth = spline(x, Rratio, xnew)

# Tratio=npy.zeros(len(x))

# for i in range(0,len(x)):
# 	Tratio[i]=epsilon_2[i]/(kB*Tstar)

#Initializing the array of densities.
# for i in range(0,len(P0)):
# 	R0[i]=density(P0[i],T0[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)	

print x
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
###    Degeneracy VERSUS x Plot  #### 
###################################################################################
# plt.plot(x,epsilon_2,'k',color='b',lw=linewidth,ls='-',label='Present Theory')
plt.plot(x,Rratio,'k',color='b',lw=linewidth,ls='-',label='Present Theory')
# plt.plot(xnew,Rratio_smooth,'k',color='b',lw=linewidth,ls='-',label='Present Theory')
plt.xlabel(r'$x_p$',fontsize=axis_size)
plt.ylabel(r'Degeneracy $g_p$',fontsize=axis_size)
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)
# fig.savefig('./'+output_folder+r'\PMMA_Olabisi_Rratio and epsilon versus x Straight Line Fit'+img_extension,dpi=240)
plt.show()
