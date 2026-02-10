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

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.0000000000000001,0.9999999999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

'''
def OwnCriteria1_Deriv_Equation(epsilon_2,Rratio,Vratio,Tg_atm,dP_dT_atm,M,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	T=Tg_atm
	dP_dT=dP_dT_atm
	P=0.101325
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar
	epsilon_2=kB*Tstarstar

	dPtilde_dTtilde=dP_dT*(Tstar/Pstar)

	F=(Rratio*exp(-Tratio**2/(Ttilde*Pratio)))/(1+Rratio*exp(-Tratio**2/(Ttilde*Pratio)))
	
	res=((((Tratio*Ttilde)/Pratio)*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((dPtilde_dTtilde)/(Rtilde))))/(Rtilde*(((Ttilde/Rtilde)*(((Rtilde)/(1-Rtilde))+(1/r)))-2)))+((((Tratio*Tstarstar)/(Pratio*T))**2)*F*(1-F))

	return res

def Rratio_From_Own_Criteria_1(M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	# Tstarstar=epsilon_2/kB
	# Tratio=Tstarstar/Tstar
	division_gap=100
	for i in range(10,20000,division_gap):
		# print i
		try:
   			# print epsilon_2
			epsilon_2 = bisect(OwnCriteria1_Deriv_Equation,i,i+division_gap,args=(Rratio,Vratio,Tg_atm,dP_dT_atm,M,Pstar,Tstar,Rstar))
			print 'Your supposed range is sussessful'
			print epsilon_2
			# print 'Rratio is'
			# print Rratio
			pass
		
		except:
   			print 'Your supposed range does not work'
   			pass

	return epsilon_2

'''

def specificHeat(P,T,R,M,**kwargs):     #Cp/mass  and **kwargs must contain "three" characteristic parameters and "three" flexibility parameters.
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	Pratio=Pstarstar/Pstar
	Tratio=Tstarstar/Tstar
	Rratio=Rstarstar/Rstar


	C1=(Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))))
	#C2_with v equal v0=(Pstar/(Rstar*Tstar))*((((((Tratio**2)*Rratio)/(Ttilde**2))*(exp(-(Tratio/Ttilde))))/((1+(Rratio*(exp(-(Tratio/Ttilde)))))**2)))
	C2=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	# C2=(Pstar/(Rstar*Tstar))*(((((   (Tratio**2)*r*Rratio/N   )/(Ttilde**2))*(exp(-(((r*Tratio)/N)/Ttilde))))/((1+(Rratio*(exp(-(((r*Tratio)/N)/Ttilde)))))**2)))

	# A=3.687e-11#0.4747231	#3.363225e-10	#3.363225e-10
	# B=0.00358531	#0.003023843	#0.004156399	#0.003992431# 0.00394548
	C3=A+B*T
	#Crotation=0.0								#Pstar/(Rstar*Tstar)
	#Tstarstarstar=700.0
	#Cvibration=2.0*(((Tstarstarstar/T)**2)*(exp(-(Tstarstarstar/T)))/((1-(exp(-(Tstarstarstar/T))))**2))#*(Pstar/(Rstar*Tstar))
	C=C1+C2+C3#+Crotation+Cvibration
	#CexcludingBending=C1+Crotation+Cvibration
	return C,C1,C2,C3

P0 = P_atm
# T0=Tg_atm
M0=M_infinity
# R0=density(P0,T0,M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M0)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

# Ptilde=P0/Pstar
# Ttilde=T0/Tstar
# Rtilde=R0/Rstar
# vtilde=1/Rtilde	
dPtilde_dT=dP_dT_atm/Pstar
dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar
##########################################################################################################

##########################################################################################################

Rratio=	1.427345691
epsilon_2= 3897.431932

# x=0.32

# Rratio=1.6425014628866696	#1.75451544356446			#PMMA gmin Grassia Data = 1.1473193			#PMMA gmin Condo Data = 2.5883053		#PS gmin = 1.642864790		#PMMA 44kilo = 63.5105956	#PMMA 140kilo = 1.20630412	#PS=1.85399320
# epsilon_2=6638.92051810446	#PMMA gmin Grassia Data = 6533.2915188		#PMMA gmin Condo Data = 7775.6803   	#PS gmin = 7123.66472	    #PMMA 44kilo = 7550.38882	#PMMA 140kilo = 6251.31324	#PS=6480.99556
Vratio=1.0
# A= 0.34997622				#PVAc 00kilo=0.34997622		#PMMA 44kilo = 3.6870e-11	#PMMA 140kilo = 3.6870e-11
# B= 0.00231786				#PVAc 00kilo=0.00231786 	#PMMA 44kilo = 0.00551656	#PMMA 140kilo = 0.00390689	
#A:	#PVAc 01kilo=2.6799e-09	#PVAc 189kilo=-0.51533657	#PC 00kilo=0.11396855		#PC 01kilo=8.4781e-11
#B:	#PVAc 01kilo=0.00287418	#PVAc 189kilo=0.00524197	#PC 00kilo=0.00306293		#PC 01kilo=0.00336962
#A:	#LPP 15kilo=0.28422215	#PS= 3.6870e-11				#PVME=0.35758515
#B:	#LPP 15kilo=0.00381332	#PS= 0.00372975				#PVME=0.00280566

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Rstarstar=Rratio*Rstar
Pratio=Tratio/Vratio
Pstarstar=Pratio*Pstar


# kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'Rratio':Rratio,'Vratio':Vratio,'A':A,'B':B,'Tg_atm':Tg_atm,'dP_dT_atm':dP_dT_atm}
# epsilon_2=Rratio_From_Own_Criteria_1(M0,**kwargs)
# print 'this my epsilon_2 is'
# print epsilon_2
# Rratio=3.79311229297
# epsilon_2=6817.9

'''
#Condo Flexibiity Parameters
Vratio=1.0    #Vratio = v/v_0
epsilon_2=5257.282680671158
# z=5.0
# Rratio=z-2.0
Rratio=7.340435509599036
Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Pratio=Tratio/Vratio
Pstarstar=Pratio*Pstar
Rstarstar=Rratio*Rstar

P0=0.101325
M0=3600000000.0
'''
# Pstarstar=0.00249935#1.2855e-05	#1.7274e-04	#5.1759e-04
# Tstarstar=11.9757060#1.58456913	#4.36981412	#6.518#1817.02246
# Rstarstar=1.1013e+40#3.999e+141	#1.3155e+80	#6.5113e+59#23.7286863
# Tratio=Tstarstar/Tstar

# N=(36000*Pstarstar)/(8.313*Rstar*Tstarstar)

# Pratio = 0.00009			#0.00009		#0.00008					 															#0.00339
# Tratio = 0			#0.07		#0.03																					#0.75
# Rratio = 5			#150000000000000000000000000000000000000000000		#150000000																				#90

# Pstarstar= Pratio*Pstar		#4236.99551			#8737.82469							#
# Tstarstar= Tratio*Tstar		#3259.82974			#5924.59463							#
# Rstarstar= Rratio*Rstar    	#6.39600739			#5.40346613							#

# P0=0.101325
# M0=36000

#Initializing the array of densities.
T0 = npy.linspace(100,600,100)
C0=npy.zeros(len(T0))		#Total Cp
C1=npy.zeros(len(T0))		#Cp Kier
C2=npy.zeros(len(T0))		#Cp at Glass Transition
C3=npy.zeros(len(T0))		#Cp = A + B*T
C4=npy.zeros(len(T0))		#Cp_kier + A + B*T  i.e. Cp excluding Glass Transition
R0=npy.zeros(len(T0))

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
output_folder = 'plot_specificHeat'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Defining linetype
# M179_line = '-'

for i in range(0,len(T0)):
	R0[i]=density(P0,T0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

for i in range(0,len(T0)):
	C0[i],C1[i],C2[i],C3[i] = specificHeat(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)


for i in range(0,len(T0)):
	C4[i]=C1[i]+C3[i]

# Plotting Theoretical Curve For Difference Values of Pressure:

'''
#==============================================================================================================
#Calculating Isobars.
#==============================================================================================================

press = ['0MPa']#,'15MPa','30MPa','45MPa','60MPa','125MPa','150MPa','175MPa','200MPa']
isobar= [0.101325]#,15,30,45,60,125,150,175,200]
#print press[3]
for i in range(0,len(press)):
	exec "P0_%s = npy.full((1,len(T0)),%s)" %(press[i],isobar[i])
	exec "C0_%s = npy.zeros(len(T0))" %(press[i])
	exec "C1_%s = npy.zeros(len(T0))" %(press[i])
	exec "C2_%s = npy.zeros(len(T0))" %(press[i])
	exec "C3_%s = npy.zeros(len(T0))" %(press[i])

	exec "P0_%s =P0_%s[0] " %(press[i],press[i])
	exec "R0_%s = calculateDensity(P0_%s[0],T0,M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (press[i],press[i])
	exec "R0_%s = R0_%s[1]" %(press[i],press[i])
	for j in range(0,len(T0)):
		exec "C0_%s[j],C1_%s[j],C2_%s[j],C3_%s[j] = specificHeat(P0_%s[j],T0[j],R0_%s[j],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)"  % (press[i],press[i],press[i],press[i],press[i],press[i])
	#exec "P%s_T = result[0]" % (press[i])
	#exec "vector_%s = findVectors(P%s_T,P0,T0_%s,R0_%s)" % (press[i],press[i],press[i],press[i])



#==============================================================================================================
#Calculating isomasss.
#==============================================================================================================

mass = ['179kgpermol','208kgpermol','248kgpermol','36kgpermol','268kgpermol']
isomass= [179000,208000,248000,36000,268000]
P0=0.0034
#print mass[3]
for i in range(0,len(mass)):
	exec "M0_%s = npy.full((1,len(T0)),%s)" %(mass[i],isomass[i])
	exec "C0_%s = npy.zeros(len(T0))" %(mass[i])
	exec "C1_%s = npy.zeros(len(T0))" %(mass[i])
	exec "C2_%s = npy.zeros(len(T0))" %(mass[i])
	exec "C3_%s = npy.zeros(len(T0))" %(mass[i])

	exec "M0_%s =M0_%s[0] " %(mass[i],mass[i])
	exec "R0_%s = calculateDensity(P0,T0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (mass[i],mass[i])
	exec "R0_%s = R0_%s[1]" %(mass[i],mass[i])
	for j in range(0,len(T0)):
		exec "C0_%s[j],C1_%s[j],C2_%s[j],C3_%s[j] = specificHeat(P0,T0[j],R0_%s[j],M0_%s[j],N=N,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)"  % (mass[i],mass[i],mass[i],mass[i],mass[i],mass[i])
	#exec "P%s_T = result[0]" % (mass[i])
	#exec "vector_%s = findVectors(P%s_T,P0,T0_%s,R0_%s)" % (mass[i],mass[i],mass[i],mass[i])
'''


# C0tilde=C0*Rstar*Tstar/Pstar
# T0tilde=T0/Tstar


# molarmass = ['179K','240K','36K']							#K=kilo; not Kelvin

# for i in range(0,len(molarmass)):
# 	exec "R0=calculateDensity(P0_%s[0],T0,M0_%s[0],Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)" % (molarmass[i],molarmass[i])
# 	exec "R0=R0[1]"	#Because caclulateDensity returns nested list whose 2nd row is list of R0 values
# 	exec "result = calculateSpecificHeat(P0_%s[0],T0,R0,M0_%s[0],N=N,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)" % (molarmass[i],molarmass[i])
# 	exec "M%s_C = result[0]" % (molarmass[i])
# 	#exec "vector_%s = findVectors(T%s_P,R0,P0_%s,R0_%s)" % (temp[i],temp[i],temp[i],temp[i])

#General line properties.
linewidth = 1
markersize = 1

arrow_ls = 'dashdot'
show_arrows = True
#print M179K_C
#==================================================================================
#P versus R plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

plt.axvline(x=Tg_atm,lw=0.5,color='k', linestyle='-.')

plt.plot(T0,C0,'k',color='r',lw=linewidth,ls='-',label='C_Total')
plt.plot(T0,C1,'k',color='b',lw=linewidth,ls='--',label='C_Kier theory')
plt.plot(T0,C2,'k',color='g',lw=linewidth,ls='-.',label='C_Glass theory')
plt.plot(T0,C3,'k',color='c',lw=linewidth,ls=':',label='C_A+BT theory')
plt.plot(T0,C4,'k',color='m',lw=linewidth,ls=':',label='C_excluding_Tg theory')

# plt.plot(T0,C0_0MPa,'k',color='r',lw=linewidth,ls='-',label='C_0MPa')
# plt.plot(T0,C0_15MPa,'k',color='b',lw=linewidth,ls='-',label='C_15MPa')
# plt.plot(T0,C0_30MPa,'k',color='g',lw=linewidth,ls='-',label='C_30MPa')
# plt.plot(T0,C0_45MPa,'k',color='c',lw=linewidth,ls='-',label='C_45MPa')
# plt.plot(T0,C0_60MPa,'k',color='m',lw=linewidth,ls='-',label='C_60MPa')
'''
plt.plot(T0,C0_179kgpermol,'k',color='r',lw=linewidth,ls='-',label='C_179kgpermol')
plt.plot(T0,C0_208kgpermol,'k',color='b',lw=linewidth,ls='-',label='C_208kgpermol')
plt.plot(T0,C0_248kgpermol,'k',color='g',lw=linewidth,ls='-',label='C_248kgpermol')
plt.plot(T0,C0_36kgpermol,'k',color='c',lw=linewidth,ls='-',label='C_36kgpermol')
plt.plot(T0,C0_268kgpermol,'k',color='m',lw=linewidth,ls='-',label='C_268kgpermol')
'''
#plt.plot(T524K_P,R0,'y')
#plt.plot(T0_179K,C0_179K,'ok',ms=markersize)#,label='179kilo experiment')
#plt.plot(T0_240K,C0_240K,'^k',ms=markersize)#,label='240kilo experiment')
plt.plot(T0_complete_Tg,C0_complete_Tg,'sk',ms=markersize)#,label='36kilo experiment')
#plt.plot(P0_524K,R0_524K,'*y')
plt.xlabel('Temperature T (K)',fontsize=axis_size)
plt.ylabel(r'Specific Heat $c_P$ ($J/g.K$)',fontsize=axis_size)
# plt.axis([250,450,1.00,2.25])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)
#minorLocator = AutoMinorLocator()
#ax.xaxis.set_minor_locator(minorLocator)
#plt.tick_params(which='both', width=1)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4)
#minorLocator = AutoMinorLocator()
#ax.yaxis.set_minor_locator(minorLocator)
#plt.tick_params(which='both', width=1)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4)
# figPUREPS.savefig('./'+output_folder+r'\pure_PS_specificHeatFitted_36kilo'+img_extension,dpi=img_dpi)

plt.show()
