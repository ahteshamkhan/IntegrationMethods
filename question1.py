# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 20:57:18 2017

@author: leoke_000
"""

# -*- coding: utf-8 -*-
"""
Numerical Methods For Physics Section 2.2 Problem 26

Ash and Leo
"""

import math
import numpy as np
import matplotlib.pyplot as plt

#Set the physical constants and other variables as global variables
g_over_l = 1.0 #constant
time = 0
irrev = 0
nstep = 30000 #constant
omega = 0

def theoretical(theta, tau, time, nstep):
    """Determines thetas using small-angle approximation"""
    timesT = [time]
    thetasT = [theta * 180/math.pi]
    
    for istep in range(1, nstep):
        time = time + tau
        period = 2*math.pi*math.sqrt(1/g_over_l)
        theta_n = theta * math.cos((2*math.pi*time)/period)
        timesT.append(time)
        thetasT.append(theta_n * 180/math.pi)
    
    plt.figure(1)
    plt.plot(timesT, thetasT)
    plt.show()
    print(timesT[-1], thetasT[-1])
    return [timesT, thetasT]


def euler_cromer(theta, tau, time, omega, nstep):
    """Determines thetas using the Euler-Cromer method"""
    timesE = [time]
    thetasE = [theta * 180/math.pi]
    ##theoryThetas = []
    for istep in range(1, nstep):
        accel = -g_over_l * math.sin(theta)
        omega_n = omega + tau*accel
        theta_n = theta + tau*omega_n
        
        omega = omega_n
        theta = theta_n
      
        time = time + tau
        
        timesE.append(time)
        thetasE.append(theta * 180/math.pi)
        
    plt.figure(2)
    plt.plot(timesE, thetasE)
    plt.show()
    print(timesE[-1], thetasE[-1])
    return [timesE, thetasE]

def leapfrog(theta, tau, time, omega, nstep):
    """Determines thetas using the Leap-Frog method"""
    timesL = [time]
    thetasL = [theta * 180/math.pi]
    
    #Uses the Euler-Cromer method to compute omega[n-1] and theta[n]
    accel = -g_over_l * math.sin(theta)
    omega = omega + accel*tau
    theta = theta + omega*tau
    time = time + tau*2
    
    timesL.append(time)
    thetasL.append(theta * 180/math.pi)

    for istep in range(2, nstep):
        
        #Computes omega[n+1]
        if istep % 2 == 0:
            accel = -g_over_l * math.sin(theta)
            omega = omega + 2*tau*accel
        
        #Computes omega[r+2]    
        else:
            theta = theta + 2*tau*omega
            time = time + 2*tau
            
            timesL.append(time)
            thetasL.append(theta * 180/math.pi)
    
    plt.figure(3)
    plt.plot(timesL, thetasL)
    plt.show()
    
    print(timesL[-1], thetasL[-1])
    return [timesL, thetasL]
    
            

    
def leapfrog_old(theta, tau, time, omega, nstep):
    """Determines thetas using the Leapfrog method"""
    timesL = []
    thetasL = []
    
    #Start the Leap-Frog method
    accel = -g_over_l * math.sin(theta)
    omega_old = omega - accel*tau
    theta_old = theta - omega*tau
    
    for istep in range(1, nstep):
        accel = -g_over_l * math.sin(theta)
        omega_n = omega_old + 2*tau*accel
        theta_n = theta_old + 2*tau*omega_n
        theta_old = theta
        theta = theta_n
        omega_old = omega
        omega = omega_n
        time = time + tau
        timesL.append(time)
        thetasL.append(theta)
        
    plt.plot(timesL, thetasL)

   
def midpoint(theta, tau, time, omega, nstep):
    """Determines thetas using the Midpoint method"""
    timesM = [0]
    thetasM = [theta * 180/math.pi]
    for istep in range(1, nstep):
        accel = -g_over_l * math.sin(theta)
        omega_n = omega + tau*accel
        theta_n = theta + tau*(omega + omega_n)/2
        omega = omega_n
        theta = theta_n
      
        time = time + tau
        timesM.append(time)
        thetasM.append(theta * 180/math.pi)
        
    
    plt.figure(4)
    plt.plot(timesM, thetasM)
    plt.show()  
    print(timesM[-1], thetasM[-1])
    return [timesM, thetasM]
    
    
def main():
    method = input("Which method would you like to use? Enter \n 1 for \
                   Small-Angle theory \n 2 for Euler-Cromer \n 3 for\
                   Leap-Frog \n 4 for Midpoint Method\n 5 for all \n")
    theta = float(input("Enter initial angle (in degrees): ")) * math.pi/180
    tau = float(input("Enter time step: "))
    nstep = int(input("Enter number of steps: "))
    
    if method == '1':
        theoretical(theta, tau, time, nstep)
    if method == '2':
        euler_cromer(theta, tau, time, omega, nstep)
    if method == '3':
        leapfrog(theta, tau, time, omega, nstep)
    if method == '4':
        midpoint(theta, tau, time, omega, nstep)
    if method == '5':
        timesT, thetasT = theoretical(theta, tau, time, nstep)
        timesE, thetasE = euler_cromer(theta, tau, time, omega, nstep)
        timesL, thetasL = leapfrog(theta, tau, time, omega, nstep)
        timesM, thetasM = midpoint(theta, tau, time, omega, nstep)
        plt.figure(5)
        plt.plot(timesT, thetasT, color = 'blue')
        plt.plot(timesE, thetasE, color = 'red')
        plt.plot(timesL, thetasL, color = 'green')
        plt.plot(timesM, thetasM, color = 'purple')
        plt.show()
    

main()