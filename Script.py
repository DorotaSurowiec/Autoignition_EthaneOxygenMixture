# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
import numpy as np
from cantera import *
from matplotlib.pylab import *
import csv

#Mechanism used for the process 
gas = Solution('gri30.xml')

#Keybord input
Tmin = 1000
Tmax = 2000
Pmin = 101325
Pmax = 405300
fimin = 0.2 
fimax = 3.8

#Number of iterations
npoints = 11
fipoints = 10

s = 0
nt = 100000 #number of time steps
dt = 10**(-6) #time step size

#Creating lists for data storage
Ti = np.zeros(npoints, 'd')
Pi = np.zeros(npoints, 'd')
fi = np.zeros(fipoints, 'd')
tim = np.zeros(nt, 'd')
temp_cas = np.zeros(nt, 'd')
dtemp_cas = np.zeros(nt-1, 'd')
Autoignition_cas = np.zeros(npoints**2 * fipoints, 'd')
FinalTemp_cas = np.zeros(npoints**2 *fipoints, 'd')

#Lists for plots
Autoignition_casTemp = np.zeros(npoints, 'd')
FinalTemp_casTemp = np.zeros(npoints, 'd')
Autoignition_casPressure = np.zeros(npoints, 'd')
FinalTemp_casPressure = np.zeros(npoints, 'd')
Autoignition_casFi = np.zeros(fipoints, 'd')
FinalTemp_casFi = np.zeros(fipoints, 'd')



for j in range(npoints):
    Ti[j]=Tmin + (Tmax-Tmin)*j/(npoints-1)
    
    for p in range(npoints):
        Pi[p] = Pmin + (Pmax-Pmin)*p/(npoints-1)
        
        for f in range(fipoints):
            fi[f] = fimin + (fimax-fimin)*f/(fipoints-1)
            no = float(1/fi[f])
            X='C2H6:0.2857 O2:'+str(no)
            gas.TPX = Ti[j], Pi[p], X #initial temperature, pressure and stoichiometry
            r = IdealGasReactor(gas) # creating the batch reactor
            sim = ReactorNet([r]) #creating a reactor network consisting of single batch reactor
            time = 0.0 #initial simulation time
            
            #Running the simulation
            for n in range(nt): #loop for nt times steps of dt seconds
                time += dt
                sim.advance(time)
                tim[n] = time
                temp_cas[n] = r.T
                     
            #catching the autoignition timing
            Dtmax=[0,0.0]
            for n in range(nt-1):
                dtemp_cas[n] = (temp_cas[n+1] - temp_cas[n])/dt
                if (dtemp_cas[n] > Dtmax[1]):
                    Dtmax[0] = n
                    Dtmax[1] = dtemp_cas[n]
            Autoignition = (tim[Dtmax[0]] + tim[Dtmax[0] + 1])/2.
           # print 'For T=' +str(Ti[j]) +'K P=' +str(Pi[p])+'Pa and fi='+str(fi[f]) +', Autoignition time=' +str(Autoignition) + '(s)'
            Autoignition_cas[s] = Autoignition*1000 #converting [s] to [ms]
            FinalTemp_cas[s] = temp_cas[nt-1]
            s += 1

#Saving data for plots
            if Pi[p] == 101325 and fi[f] == 1:
                FinalTemp_casTemp[j] = temp_cas[nt-1]
                Autoignition_casTemp[j] = Autoignition*1000
            if Ti[j] == 1100 and fi[f] == 1:
                FinalTemp_casPressure[p] = temp_cas[nt-1]
                Autoignition_casPressure[p] = Autoignition*1000
            if Pi[p] == 101325 and Ti[j] == 1100:
                FinalTemp_casFi[f] = temp_cas[nt-1]
                Autoignition_casFi[f] = Autoignition*1000



#Data import to csv file        
s = 0
csv_file = 'DS_Autoignition_EthaneOxygen.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Initial temperature', 'Pressure', 'Fi', 'Autoignition time', 'Final Temperature'])
    for i in range(npoints):
        writer.writerow([Ti[i]])
        for n in range(npoints):
            writer.writerow(['', Pi[n]])
            for k in range(fipoints):
                writer.writerow(['', '', fi[k], Autoignition_cas[s], FinalTemp_cas[s]])
                s += 1

print 'output written to ' +csv_file

#PLOTS

#Autoign_time(InitialTemp)

plot(1000/Ti,Autoignition_casTemp,'-',color='blue')
xlabel(r'Temp [1000/K]',fontsize=20)
ylabel("Autoignition [ms]")
title(r'Autoignition of $C_{2}H_{6}$ + P=1atm + Oxygen at $\Phi$=1', fontsize=22,horizontalalignment='center')
axis([0.5,1.1,0.0,100.0])
grid()
savefig('Autoign_inittemp.png',bbox_inches='tight')

#FinalTemp(InitialTemp)

plot(1000/Ti,FinalTemp_casTemp,'-',color='blue')
xlabel(r'Temp [1000/K]',fontsize=20)
ylabel("FinalTemp [K]")
title(r'Autoignition of $C_{2}H_{6}$ + P=1atm + Oxygen at $\Phi$=1', fontsize=22,horizontalalignment='center')
axis([0.5,1.1,1000,4000])
grid()
savefig('Finaltemp_temp.png',bbox_inches='tight')

#Autoign_time(pressure)

plot(Pi,Autoignition_casPressure,'-',color='blue')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("Autoignition [ms]")
title(r'Autoignition of $C_{2}H_{6}$ + T=1100K + Oxygen at $\Phi$=1', fontsize=22,horizontalalignment='center')
axis([100000,410000,0.5,40.0])
grid()
savefig('Autoign_pressure.png',bbox_inches='tight')

#FinalTemp(pressure)

plot(Pi,FinalTemp_casPressure,'-',color='blue')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("Final Temp[K]")
title(r'Autoignition of $C_{2}H_{6}$ + T=1100K + Oxygen at $\Phi$=1', fontsize=22,horizontalalignment='center')
axis([100000,410000,3000,4000])
grid()
savefig('Finaltemp_pressure.png',bbox_inches='tight')

#Autoign_time(Fi)

plot(fi,Autoignition_casFi,'-',color='blue')
xlabel(r'Fi',fontsize=20)
ylabel("Autoignition [ms]")
title(r'Autoignition of $C_{2}H_{6}$ + P=1atm + T=1100K', fontsize=22,horizontalalignment='center')
axis([0.2,4.0,1.0,50.0])
grid()
savefig('Autoign_fi.png',bbox_inches='tight')

#FinalTemp(Fi)

plot(fi,FinalTemp_casFi,'-',color='blue')
xlabel(r'Fi',fontsize=20)
ylabel("Final Temp[K]")
title(r'Autoignition of $C_{2}H_{6}$ + P=1atm + T=1100K', fontsize=22,horizontalalignment='center')
axis([0.2,4.0,2000,4000])
grid()
savefig('Finaltemp_fi.png',bbox_inches='tight')
