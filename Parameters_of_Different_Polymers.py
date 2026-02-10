#Sanchez-Lacombe parameters for the pure POLYMERS.

#======================================================
#Characteristic Parameters
#======================================================
# Pstar = [MPa]
# Tstar = [K]
# Rstar = [g/cm3]

#############################################################################################
# Polymer_Type='PVME' #PS or PMMA or DME or LPP or BPP or PLA or LDPE
# Reference='Roland'  #Condo or Kier
#############################################################################################

P_atm = 0.101325
M_infinity = 9E9
# global Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper

def Parameters_of_Different_Polymers(**kwargs):
	
    for key,value in kwargs.items():
		exec "%s='%s'" % (key,value)
    
    print 'Polymer Type', Polymer_Type, 'Referenced from', Reference

    #################################################################################
    #PMMA Values Regressed By Myself:
    #################################################################################

    if Polymer_Type=='PMMA' and Reference=='Self_Schmidt':
        # Data 1 Very Good Matching 1
        # Paper: Good 2000 PVT of PMMA Schmidt ma991722h
        Pstar=577.933210    #+/- 22.7932121 (3.94%)
        Tstar=677.583131    #+/- 7.01760682 (1.04%)
        Rstar=1.27392163    #+/- 0.00421348 (0.33%)
        M_infinity=75000
        #Tg Data Given By Schimidt
        Tg_atm=378.15
        dTg_dP_atm = 0.24            #Self-Regressed From Internet
        Pg_exp=[0.101325,50.0,100.0,150.0,200.0]
        Tg_exp=[378.15,390.82,405.33,416.39,424.52]
        P_upper=150.0
        T_upper=416.39

    if Polymer_Type=='PMMA' and Reference=='Self_Grassia':
        # Data 2 Very Good Matching 2
        # Paper: Best 2011 PVT and Tg of PMMA Grassia 1-s2.0-S0022309310005338-main
        Pstar=562.397384    #+/- 7.52128331 (1.34%)
        Tstar=654.343418    #+/- 2.72548523 (0.42%)
        Rstar=1.28146781    #+/- 0.00153252 (0.12%)
        M_infinity=120000
        #Tg Data Given By Grassia
        Tg_atm = 352.0        #Given in Paper. However, Seems Abnormally Low Tg             
        dTg_dP_atm = 0.3      #Given Value in Paper is: 0.3        #Unit: K/MPa, "It is straight line fit" upto 150MPa
        Pg_exp=[0.101325,10.0,30.0,60.0,80.0,100.0,120.0,150.0]
        Tg_exp=[352.0,356.0,363.5,373.0,379.5,385.0,390.0,397.5]
        P_upper=150.0
        T_upper=397.5

    if Polymer_Type=='PMMA' and Reference=='Self_Walsh':
        # Data 3 Overlap and Significant Matching but Not Complete Matching
        # Paper: Old 1992 PVT of PMMA PII_ 0032-3861(92)90694-R
        Pstar=564.823646    #+/- 17.9901923 (3.19%)
        Tstar=667.892308    #+/- 5.07983058 (0.76%)
        Rstar=1.28420005    #+/- 0.00323659 (0.25%)
        #Tg Data Given By Walsh
        Tg_atm=376.9
        dTg_dP_atm = 0.2174   #Self-Calculated Since Data is Already Straight Line
        Pg_exp=[0.101325,39.23,78.45,117.68,156.91,196.13]
        Tg_exp=[376.9,386.24,396.15,405.18,413.01,419.54]
        P_upper=196.13
        T_upper=419.54

    if Polymer_Type=='PMMA' and Reference=='Self_Wen':
        # Data 4 Neither Matching Nor Overlap
        # Paper: Good 2001 PVT and Tg of PMMA Wen ma010023d
        Pstar=479.377365    #+/- 31.3221097 (6.53%)
        Tstar=709.595728    #+/- 7.99245188 (1.13%)
        Rstar=1.28295203    #+/- 0.00451412 (0.35%)
        M_infinity=387000
        #Tg Data Given By Wen
        Tg_atm=382.20
        dTg_dP_atm = 0.1515    #Self-Calculated Since Data is Already Straight Line
        Pg_exp=[0.101325,60.0,120.0,180.0]
        Tg_exp=[382.20,389.84,400.68,409.46]
        P_upper=180.0
        T_upper=409.46

    #################################################################################
    #PC Values Regressed By Myself:
    #################################################################################

    if Polymer_Type=='PC' and Reference=='Self_Aravind':
        # Data 1 Matching 1
        # Paper: Best 2012 PVT and Tg of PC Aravind 1-s2.0-S014294181100153X-main
        Pstar=477.201600      #+/- 6.50120144 (1.36%) (init = 60)
        Tstar=745.977589      #+/- 2.32835065 (0.31%) (init = 50)
        Rstar=1.28039144      #+/- 0.00124800 (0.10%) (init = 2)
        #Tg Data Given By Aravind
        Tg_atm =413.19
        dTg_dP_atm =0.3088      #Self-Calculated. Since Data Is Already Straight Line
        Pg_exp=[0.101325,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0]
        Tg_exp=[413.19,417.10,421.11,424.52,427.93,431.64,435.40,438.12,441.14,444.65,447.81,450.88,453.45,456.17,459.04,462.06,464.83,467.15,469.82,472.30,474.92]
        P_upper=200.0
        T_upper=474.92

    if Polymer_Type=='PC' and Reference=='Self_Sato':
        # Data 2 Matching 2
        # Paper: Good 1997 PVT of PC (SICI)1097-4628(19971003)66_1_141__AID-APP17_3.0.CO;2-4
        Pstar=500.834279      #+/- 11.7676771 (2.35%) (init = 60)
        Tstar=751.670242      #+/- 4.73421735 (0.63%) (init = 50)
        Rstar=1.28175322      #+/- 0.00261206 (0.20%) (init = 2)
        M_infinity=60000
        #Tg Data Given By Aravind  (No Tg data given by "Sato")
        # Tg_atm =413.19
        # dTg_dP_atm =0.3088      #Self-Calculated. Since Data Is Already Straight Line
        # Pg_exp=[0.101325,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0]
        # Tg_exp=[413.19,417.10,421.11,424.52,427.93,431.64,435.40,438.12,441.14,444.65,447.81,450.88,453.45,456.17,459.04,462.06,464.83,467.15,469.82,472.30,474.92]
        # P_upper=200.0
        # T_upper=474.92
        #Experiment Data PC; Zoller Paper
        Tg_atm = 423.4           #Zoller Paper: Unit K
        dTg_dP_atm = 0.530        #Zoller Paper, Unit: K/MPa , Linear fit value upto 60MPa given by Zoller = 0.530  
        Pg_exp=[0.101325,9.7,19.9,29.7,39.5,49.1,59.3,69.3,79.1,89.1,98.9,108.7,118.7,128.5,137.9,147.9,157.9,167.9,177.3]
        Tg_exp=[423.4,429.3,433.4,436.8,444.5,450.8,455.5,459.6,464.8,469.5,472.8,475.8,481.0,484.3,486.9,490.0,490.4,494.9,496.3]
        #Do not take more data linear fit is until 60MPa
        P_upper=59.3
        T_upper=455.5

    if Polymer_Type=='PC' and Reference=='Self_Rudolph':
        # Data 3 Overlap But Not Matching
        # Paper: Best 2016 PVT and Tg of PC Rudolph Rudolph2016_Article_WLFModelForThePressureDependen
        Pstar=378.30332       #+/- 11.5504425 (3.05%) (init = 60)
        Tstar=625.52729       #+/- 3.99028295 (0.64%) (init = 50)
        Rstar=1.33244470      #+/- 0.00373602 (0.28%) (init = 2)
        M_infinity=30500
        #Tg Data Given By Rudolph
        Tg_atm = 409.78
        dTg_dP_atm = 0.4212
        Pg_exp=[0.101325,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0]
        Tg_exp=[409.78,418.81,429.11,437.25,446.85,454.87,463.20,470.19,477.95,485.20,493.98]
        P_upper= 200.0
        T_upper= 493.98

    if Polymer_Type=='PC' and Reference=='Self_Kikuchi':
        # Data 4 Neither Matching Nor Overlap
        # Paper: Good 2003 PVT and Tg of PS and PC Kikuchi Thermal Conductivity
        Pstar=899.363789      #+/- 14.6137301 (1.62%) (init = 60)
        Tstar=659.224862      #+/- 2.94947589 (0.45%) (init = 50)
        Rstar=1.22050074      #+/- 0.00213321 (0.17%) (init = 2)
        M_infinity=73700
        #Tg Data Given By Kikuchi
        Tg_atm = 419.47
        dTg_dP_atm = 0.2505
        Pg_exp=[0.101325,10.0,50.0,100.0,150.0,200.0]
        Tg_exp=[419.47,422.96,433.14,443.01,452.89,473.11]
        P_upper= 200.0
        T_upper= 473.11

    #################################################################################
    #PS Values Regressed By Myself:
    #################################################################################

    if Polymer_Type=='PS' and Reference=='Self_Kier':
        # Data # 1: Matching But Not Exactly Kier Zoller Wash Data
        Pstar=  430.896437      #+/- 9.16502716 (2.13%) (init = 60)
        Tstar=  690.581024      #+/- 2.83119250 (0.41%) (init = 50)
        Rstar=  1.11801998      #+/- 0.00128308 (0.11%) (init = 2)
        M_infinity=110000
        #Tg Data Given By Kikuchi 
        Tg_atm = 380.00
        dTg_dP_atm = 0.208
        Pg_exp=[0.101325,10.0,50.0,100.0,150.0,200.0]
        Tg_exp=[380.00,383.05,392.72,402.75,412.42,423.19]
        P_upper= 200.0
        T_upper= 423.19

    if Polymer_Type=='PS' and Reference=='Self_Kikuchi':
        # Data#2: Overlap But Not Matching
        # Paper: Good 2003 PVT and Tg of PS and PC Kikuchi Thermal Conductivity
        Pstar=  624.048947      #+/- 16.3388656 (2.62%) (init = 60)
        Tstar=  623.692994      #+/- 4.00763066 (0.64%) (init = 50)
        Rstar=  1.13205301      #+/- 0.00257298 (0.23%) (init = 2)
        M_infinity=25000
        #Tg Data Given By Kikuchi 
        Tg_atm = 380.00
        dTg_dP_atm = 0.208
        Pg_exp=[0.101325,10.0,50.0,100.0,150.0,200.0]
        Tg_exp=[380.00,383.05,392.72,402.75,412.42,423.19]
        P_upper= 200.0
        T_upper= 423.19

    if Polymer_Type=='PS' and Reference=='Self_Park':
        # Data#3: Matching But Not Exactly 2004 Very few points
        # Paper: Good 2004 PVT of PS and PP adv.20020
        Pstar=  265.473727      #+/- 24.9282071 (9.39%) (init = 275.8)
        Tstar=  798.859724      #+/- 14.5314703 (1.82%) (init = 803.1)
        Rstar=  1.07238654      #+/- 0.00508400 (0.47%) (init = 1.072)
        M_infinity=238000
        #Tg Data Given By 
        Tg_atm = 372.0
        #Tg Data Given By Grassia
        Tg_atm = 355.0
        dTg_dP_atm = 0.3554
        Pg_exp=[0.101325,10.0,30.0,60.0,80.0,100.0,120.0,150.0]
        Tg_exp=[355.0,360.26,370.87,383.95,391.13,397.45,403.02,410.13]
        P_upper= 150.0
        T_upper= 410.13

    if Polymer_Type=='PS' and Reference=='Self_Grassia':
        # Data#4: Matching But Not Exactly 2011 Best
        # Paper: Best 2011 PVT Tg of PS app.34789
        Pstar=  462.343780      #+/- 6.99594722 (1.51%) (init = 60)
        Tstar=  647.334764      #+/- 2.98020363 (0.46%) (init = 50)
        Rstar=  1.14149853      #+/- 0.00157953 (0.14%) (init = 2)
        M_infinity=290000
        #Tg Data Given By Grassia
        Tg_atm = 355.0
        dTg_dP_atm = 0.3554
        Pg_exp=[0.101325,10.0,30.0,60.0,80.0,100.0,120.0,150.0]
        Tg_exp=[355.0,360.26,370.87,383.95,391.13,397.45,403.02,410.13]
        P_upper= 150.0
        T_upper= 410.13

    #################################################################################
    #Values Directly Taken From Other's Papers:
    #################################################################################

    if Polymer_Type=='PS' and Reference=='Kier':
        # For PS Ref: Kier 
        Pstar = 421.8 
        Tstar = 687.8 
        Rstar = 1.118
        Tg_atm = 374.0         #Ref[54] of Condo Paper: Unit K
        dTg_dP_atm = 0.316     #Ref[54] of Condo Paper: Unit K/MPa
        
        #Experiment Data PS; Condo Ref [54]
        Pg_exp=[0.101325,40,60,80,120,160]
        Tg_exp=[373.0,388.4,392.7,402.8,413.2,428.8]
        P_upper=160.0

    if Polymer_Type=='PS' and Reference=='Quach':
        # For PS Ref: Quach Paper
        Pstar = 357.0
        Tstar = 735.0
        Rstar = 1.105
        Tg_atm = 374.0          #Ref[54] of Quach Paper: Unit K
        dTg_dP_atm = 0.316       #Value in Paper: 0.316      #Ref[54] of Quach Paper, Unit: K/MPa   

        #Experiment Data PS; Quach Ref [54]
        Pg_exp=[0.101325,40,60,80,120,160]
        Tg_exp=[374.0,388.4,392.7,402.8,413.2,428.8]
        P_upper=160.0
        T_upper=428.8

    if Polymer_Type=='PC' and Reference=='Zoller':
        # For PC Ref: Zoller Paper, A Studey of PVT Relationships of Four Related Amorphous Polymers
        Pstar = 574.4             #From huge List of SL EOS Parameters
        Tstar = 728.0             #From huge List of SL EOS Parameters
        Rstar = 1.2925            #From huge List of SL EOS Parameters

        #Experiment Data PC; Zoller Paper
        Tg_atm = 423.4           #Zoller Paper: Unit K
        dTg_dP_atm = 0.530        #Zoller Paper, Unit: K/MPa , Linear fit value upto 60MPa given by Zoller = 0.530  
        Pg_exp=[0.101325,9.7,19.9,29.7,39.5,49.1,59.3,69.3,79.1,89.1,98.9,108.7,118.7,128.5,137.9,147.9,157.9,167.9,177.3]
        Tg_exp=[423.4,429.3,433.4,436.8,444.5,450.8,455.5,459.6,464.8,469.5,472.8,475.8,481.0,484.3,486.9,490.0,490.4,494.9,496.3]
        #Do not take more data linear fit is until 60MPa
        P_upper=59.3
        T_upper=455.5

    if Polymer_Type=='PMMA' and Reference=='Grassia':
        # For P*T*R* are from Condo Paper
        Pstar = 503.0
        Tstar = 696.0
        Rstar = 1.269
        Tg_atm = 352.0             
        dTg_dP_atm = 0.3      #Given Value in Paper is: 0.3        #Unit: K/MPa, "It is straight line fit" upto 150MPa

        #Experiment Data PMMA; Luigi Grassia: Isobaric and isothermal glass transition of PMMA
        Pg_exp=[0.101325,10.0,30.0,60.0,80.0,100.0,120.0,150.0]
        Tg_exp=[352.0,356.0,363.5,373.0,379.5,385.0,390.0,397.5]
        P_upper=150.0
        T_upper=397.5

    if Polymer_Type=='PMMA' and Reference=='Olabisi':
        # For PMMA Ref: Olabisi Paper
        Pstar = 503.0
        Tstar = 696.0
        Rstar = 1.269
        Tg_atm = 378.0              #Ref[53] of Olabisi Paper
        dTg_dP_atm = 0.236          #Ref[53] of Olabisi Paper, Unit: K/MPa

        #Experiment Data PMMA; Olabisi Ref [53]
        Pg_exp=[0.101325,30,40,80,120,140,180]
        Tg_exp=[378.0,386.5,386.5,397.5,408.1,408.1,419.7]
        P_upper=180.0
        T_upper=419.7

    if Polymer_Type=='PVAc' and Reference=='Sandberg':
        # For PVAc Ref: Roland Paper, Dynamic properties of polyvinylmethylether near the glass transition
        Pstar = 504.2             #From huge List of SL EOS Parameters
        Tstar = 592.0             #From huge List of SL EOS Parameters
        Rstar = 1.2822            #From huge List of SL EOS Parameters

        #Experiment Data PVAc; Sandberg Paper #This data does not seems right to me. Shifted by 10K.
        Tg_atm = 319.0           #Sandberg Paper: Unit K
        dTg_dP_atm = 0.264        #Sandberg Paper, Unit: K/MPa, in low P limit. So take only low pressure values, linear fit value upto 80 MPa.   
        Pg_exp=[0.101325,26.142,79.299,162.381,241.831,353.225,493.251]
        Tg_exp=[319.0,325.53,338.30,355.91,370.70,388.63,407.90]
        #Take only low pressure values. Data has curvature.   
        P_upper=162.381
        T_upper=355.91

    if Polymer_Type=='PVAc' and Reference=='Roland':
        # For PVAc Ref: Roland Paper, Dynamic properties of polyvinylmethylether near the glass transition
        Pstar = 504.2             #From huge List of SL EOS Parameters
        Tstar = 592.0             #From huge List of SL EOS Parameters
        Rstar = 1.2822            #From huge List of SL EOS Parameters
        #Use Experiment Data PVAc; Roland Paper
        Tg_atm = 311.0           #Roland Paper: Unit K
        dTg_dP_atm = 0.216     #My value=0.216 upto 150MPa,     #Roland Paper value is 0.25 in limit of P=0, Unit: K/MPa
        #Experiment Data PVAc; Roland Paper
        Pg_exp=[0.101325,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0]
        Tg_exp=[311.0,323.0,333.0,343.5,351.5,359.5,366.5,373.5,380.0]
        #My slope value=0.216 upto 150MPa
        P_upper=150.0
        T_upper=343.5

    if Polymer_Type=='PVME' and Reference=='Casalini':
        # For PVME Ref: Casalini Paper, Dynamic properties of polyvinylmethylether near the glass transition
        Pstar = 463.0             #From huge List of SL EOS Parameters
        Tstar = 567.0             #From huge List of SL EOS Parameters
        Rstar = 1.1198            #From huge List of SL EOS Parameters
        Tg_atm = 247.60           #Casalini Paper: Unit K
        dTg_dP_atm = 0.149 #My Value=0.149 upto 180MPa      #Casalini Paper=0.177 in P=0 limit, Unit: K/MPa

        #Experiment Data PVME; Casalini Paper
        #Curve has significant curvature even at low values of pressure
        Pg_exp=[0.101325,50.16,111.50,177.70,249.75,309.13,375.32,441.51,497.96,556.36,622.55,657.59,687.76]
        Tg_exp=[247.60,256.22,265.21,274.08,282.69,289.25,296.19,302.48,307.49,312.12,317.64,320.34,322.40]

        # Pg_exp=[0.101325,250,375,690]
        # Tg_exp=[247.6,282.5,296.0,322.5]
        #Curve has significant curvature even at low values of pressure
        P_upper=177.70
        T_upper=274.08

    if Polymer_Type=='BPP' and Reference=='Hollander':
        # For Brached PP Ref: Kier 
        Pstar = 356.4 
        Tstar = 656.0 
        Rstar = 0.8950 
        Tg_atm = 251.2          #Ref: Actactic PP at High Pressure. Deturon NMR of Glass Temp; Hollander 2001 
        dTg_dP_atm = 0.158      #My Value=0.158 upto 100.4MPa   #In P=0 limit value is given to be 0.158   #Ref: Actactic PP at High Pressure. Deturon NMR of Glass Temp; Hollander 2001 
    
        #Experiment Data PP; Ref: Actactic PP at High Pressure. Deturon NMR of Glass Temp; Hollander 2001
        Pg_exp=[0.101325,0.4,50.5,100.4,199.5,500.3]
        Tg_exp=[251.2,250.0,262.0,267.0,279.0,311.0]
        #My slope value=0.158 upto 100.4MPa 
        P_upper=100.4
        T_upper=267.0

    if Polymer_Type=='BPP' and Reference=='Passaglia':
        # For Brached PP Ref: Kier 
        Pstar = 356.4 
        Tstar = 656.0 
        Rstar = 0.8950 
        Tg_atm = 243.5        #Ref: "Variation of Glass Temperature With Pressure in Polypropylene,Passaglia,1962"
        dTg_dP_atm = 0.204    #Note:It is straight line fit slope on whole data    #Ref: "Variation of Glass Temperature With Pressure in Polypropylene,Passaglia,1962"
        #Experiment Data Linear Polypropylene; Ref: "Variation of Glass Temperature With Pressure in Polypropylene,Passaglia,1962"
        Pg_exp=[0.101325,15,30,40,50,70]
        Tg_exp=[243.5,249.0,251.3,252.3,254.0,258.3]
        P_upper=70.0
        T_upper=258.3

    if Polymer_Type=='LPP' and Reference=='Hollander':
        # For Linear PP Ref: Kier 
        Pstar = 316.2 
        Tstar = 662.8 
        Rstar = 0.8685
        Tg_atm = 251.2          #Ref: Actactic PP at High Pressure. Deturon NMR of Glass Temp; Hollander 2001 
        dTg_dP_atm = 0.158      #My Value=0.158 upto 100.4MPa   #In P=0 limit value is given to be 0.158   #Ref: Actactic PP at High Pressure. Deturon NMR of Glass Temp; Hollander 2001 
    
        #Experiment Data PP; Ref: Actactic PP at High Pressure. Deturon NMR of Glass Temp; Hollander 2001
        Pg_exp=[0.101325,0.4,50.5,100.4,199.5,500.3]
        Tg_exp=[251.2,250.0,262.0,267.0,279.0,311.0]
        #My slope value=0.158 upto 100.4MPa 
        P_upper=100.4
        T_upper=267.0

    if Polymer_Type=='LPP' and Reference=='Passaglia':
        # For Linear PP Ref: Kier 
        Pstar = 316.2 
        Tstar = 662.8 
        Rstar = 0.8685
        Tg_atm = 243.5        #Ref: "Variation of Glass Temperature With Pressure in Polypropylene,Passaglia,1962"
        dTg_dP_atm = 0.204    #Note:It is straight line fit slope on whole data    #Ref: "Variation of Glass Temperature With Pressure in Polypropylene,Passaglia,1962"
        #Experiment Data Linear Polypropylene; Ref: "Variation of Glass Temperature With Pressure in Polypropylene,Passaglia,1962"
        Pg_exp=[0.101325,15,30,40,50,70]
        Tg_exp=[243.5,249.0,251.3,252.3,254.0,258.3]
        P_upper=70.0
        T_upper=258.3

    if Polymer_Type=='DME' and Reference=='Kier':
        # For DME Ref: Kier 
        Pstar = 313.8 
        Tstar = 450.0 
        Rstar = 0.8146 
        Tg_atm = 1.0
        dTg_dP_atm = 0.0

    if Polymer_Type=='LDPE' and Reference=='Kier':
        # For Low Density PE Ref: Kier 
        Pstar = 407.5 
        Tstar = 586.6 
        Rstar = 0.9271
        Tg_atm = 1.0
        dTg_dP_atm = 0.0

    if Polymer_Type=='PLA' and Reference=='Kier':
        # For PLA Ref: Kier 
        Pstar = 598.4 
        Tstar = 617.3 
        Rstar = 1.347
        Tg_atm = 1.0
        dTg_dP_atm = 0.0

    return (Pstar,Tstar,Rstar,Tg_atm,dTg_dP_atm,Pg_exp,Tg_exp,P_upper,T_upper)
