#!/bin/bash
# -*- coding: utf-8 -*-
# 20180428T0206Z
# Tested on Python 3.6.5
# the formulas of the blackhole() function were obtained from the site: http://xaonon.dyndns.org/hawking/

from os import system
import math, sys

# global constants
lum_sol = 3.828e+26 # solar luminosity
rad_sol = 695500000 # solar radius in meters
temp_sol = 5778 # solar temperature in Kelvin
mass_sol = 1.98855e+30 # solar mass
sig = 5.670367e-8 # Stefan-Boltzmann constant
AU = 149597870691 # Astronomical Unit
c = 299792458 # speed of light
hbar = 1.0545718e-34 # reduced Planck constant
G = 6.67408e-11 # gravitational constant
k = 1.38066e-23 # Boltzmann constant

#
# luminosity and habitable zone of a star
def luminosity():
    print()

    # luminosity of star calculus
    rad_star = float(input("Enter the radius of the star in unit of solar radius: "))
    temp_star = float(input("Enter the star temperature in Kelvin: "))
    radstar1 = rad_star*rad_sol
    luminosity_in_watts = ((radstar1/rad_sol)**2*(temp_star/temp_sol)**4)*lum_sol
    luminosity_in_lum_sol = (((radstar1/rad_sol)**2*(temp_star/temp_sol)**4)*lum_sol)/lum_sol

    print("The luminosity of the star is:", luminosity_in_watts, "W")
    print("The luminosity of star in a unit of solar luminosity:", luminosity_in_lum_sol)
    print()
    print("Habitable Zone calculus:")
    print()

    # habitable zone calculus
    boiling_temperature = float(input("Enter the boiling temperature in Kelvin of chemical element: "))
    melting_temperature = float(input("Enter the melting temperature in Kelvin of chemical element: "))
    inner_border_HZ = (((0.5*luminosity_in_watts)/(4*math.pi*sig*boiling_temperature**4))**(0.5))/AU
    outside_border_HZ = (((0.5*luminosity_in_watts)/(4*math.pi*sig*melting_temperature**4))**(0.5))/AU

    print()
    print("Inner border HZ:" , inner_border_HZ , "AU")
    print("Outside border HZ:" , outside_border_HZ , "AU")

# characteristics of a black hole
def blackhole():
    print()
    mass = float(input("Enter the mass of black hole in unit of solar mass: "))
    bhmass = mass*mass_sol
    radius = (2*G*bhmass)/(pow(c,2))
    radius_AU = radius/AU
    diameter_AU = 2*radius_AU
    diameter = 2*radius
    surface_area = ((16*math.pi*pow(G,2)*pow(bhmass,2)))/(pow(c,4))
    surface_gravity = (pow(c,4))/(bhmass*4*G)
    surface_tides = (pow(c,6))/(pow(bhmass,2)*pow(G,2)*4)
    entropia = (4*math.pi*G*pow(bhmass,2))/(hbar*c*math.log(10))
    temp = (hbar*pow(c,3))/(bhmass*8*k*math.pi*G)
    luminosidade = (hbar*pow(c,6))/(pow(bhmass,2)*15360*math.pi*pow(G,2))
    life_time = (pow(bhmass,3)*5120*math.pi*pow(G,2))/(hbar*pow(c,4))
    billion_years = life_time/3.154e+16
    
    print()
    print("Schwarzchild radius of black hole in AU: ", radius_AU, "AU")
    print("Diameter in AU: ", diameter_AU, "AU")
    print("Radius in meters: ", radius, "m")
    print("Diameter in meters: ", diameter, "m")
    print("Surface area: ", surface_area, "m²")
    print("Gravitational acceleration: ", surface_gravity, "m/s²")
    print("Tidal force: ", surface_tides, "m/s²/m")
    print("Entropy: ", entropia)
    print("Temperature: ", temp, "K")
    print("Luminosity: ", luminosidade, "W")
    print("Lifetime in billion of years: ", billion_years)

# light-year to parsec
def ly_to_pc():
    print()
    ly = float(input("Enter the distance n light-year: "))
    ly_to_pc = ly/3.2615637769
    print()
    print(ly, "ly =", ly_to_pc, "pc")

# parsec to light-year
def pc_to_ly():
    print()
    pc = float(input("Enter the distance in parsec: "))
    pc_to_ly = pc*3.2615637769
    print()
    print(pc, "pc =", pc_to_ly, "ly")

# meters to light-second
def m_to_ls():
    print()
    m = float(input("Enter the distance in meters: "))
    m_to_ls = m/299792458
    print()
    print(m, "m =", m_to_ls, "light-second")

# light-second to meters
def ls_to_m():
    print()
    light_second = float(input("Enter the distance in light-second: "))
    ls_to_m = light_second*c
    print()
    print(light_second, "light-second =", ls_to_m, "m")

# light-year to meter
def ly_to_m():
    print()
    light_year = float(input("Enter the distance in light-year: "))
    l_to_m = 9460730472580800*light-year
    print()
    print(light-year, "ly =", l_to_m, "m")

# meters to light-year
def m_to_ly():
    print()
    dist_m = float(input("Enter the distance in meter: "))
    m_to_ly = dist_m/1.0570008340246154e-16
    print()
    print(dist_m, "m =", m_to_ly, "ly")

# time dilation
def time_dilation():
    print()
    time = float(input("Enter travel time in seconds: "))
    velocity = float(input("Enter speed in meters per second: "))
    result = (time)/((1-(velocity**2/c**2)))**(1/2)
    print()
    print("The dilation of time traveling to", velocity, "m/s =", result, "s")

def mass_to_energy():
    print()
    mass = float(input("Enter the mass in kg: "))
    mass_to_energy = mass*pow(c,2)
    print(mass, "kg =", mass_to_energy, "J")
    
def energy_mass():
    print()
    energy = float(input("Digite o valor da energia em Joule: "))
    energy_to_mass = energy/(pow(c,2))
    print(energy, "Joule é igual a", energy_to_mass, "kg")
    
# Scientific calculator functions

# greatest common divisor
def mcd():
    print()
    first_value = int(input("Enter the first value: "))
    second_value = int(input("Enter the second value: "))
    result = math.gcd(first_value,second_value)
    print("The greatest common divisor", first_value, "and", second_value, " = ", result)

# square root
def sqrt():
    print()
    rooting = int(input("Enter the value of rooting: "))
    result = math.sqrt(rooting)
    print()
    print("Sware root of", rooting," = ", result)

# logarithm 
def log():
    print()
    logging = input("Enter the value of logging: ")
    base = input("Enter the value of base: ")
    result = math.log(logging[base])
    print("Logarithm of", logging," base ", base, " = ", result)

# power
def power():
    print()
    base = float(input("Enter the base: "))
    exponent = float(input("Enter the expoent: "))
    result = math.pow(base, exponent)
    print(base, "^", exponent, "=", result)
    
def complex_multiplication():
    print()
    complex_0 = complex(input("Enter the first complex number in the form a+bj: "))
    complex_1 = complex(input("Enter the second complex number in the form a+bj: "))
    result = complex_0 * complex_1
    print(complex_0 ,'*',complex_1 ,'=', result)
    
def quaternion_multiplication():
    print()
    print("Multiplication of two Quaternions in the form q = u+xi+yj+zk and p = a+bi+cj+dk")
    u = float(input("Enter the variable u: "))
    x = float(input("Enter the variable x: "))
    y = float(input("Enter the variable y: "))
    z = float(input("Enter the variable z: "))
    a = float(input("Enter the variable a: "))
    b = float(input("Enter the variable b: "))
    c = float(input("Enter the variable c: "))
    d = float(input("Enter the variable d: "))
    real = (u*a-x*b-y*c-z*d)
    i = (u*b+x*a+y*d-z*c)
    j = (u*c-x*d+y*a+z*b)
    k = (u*d+x*c-y*b+z*a)
    norm = (real**2+i**2+j**2+k**2)**(1/2)
    print()
    print('(',real,',' ,i,'i,' ,j,'j,',k,'k)')
    print()
    print("Norma: ", norm)
    
def quaternion_division():
    print("Division of two quaternions in the form q^-1 * p, q0+iq1+jq2+kq3 and r0+ir1+jr2+kr3")
    print()
    q0 = float(input("Enter the variable u q0: "))
    q1 = float(input("Enter the variable u q1: "))
    q2 = float(input("Enter the variable u q2: "))
    q3 = float(input("Enter the variable u q3: "))
    r0 = float(input("Enter the variable u r0: "))
    r1 = float(input("Enter the variable u r1: "))
    r2 = float(input("Enter the variable u r2: "))
    r3 = float(input("Enter the variable u r3: "))
    print()
    t0 = (r0*q0+r1*q1+r2*q2+r3*q3)/(r0**2+r1**2+r2**2+r3**2)
    t1 = ((r0*q1-r1*q0-r2*q3+r3*q2)/(r0**2+r1**2+r2**2+r3**2))*(-1)
    t2 = ((r0*q2+r1*q3-r2*q0-r3*q1)/(r0**2+r1**2+r2**2+r3**2))*(-1)
    t3 = ((r0*q3-r1*q2+r2*q1-r3*q0)/(r0**2+r1**2+r2**2+r3**2))*(-1)
    print("Result:",t0,',',t1,'i,',t2,'j,',t3,'k')
    
# function to exit the program
def exit():
    print()
    exit = input("Exit? y ou n: ")
    if exit == "y":
        sys.exit(0)
    elif exit == "n":
        main()

# main function
def main():
    print()
    print("Enter 1 to Astronomical Calculator")
    print("Enter 2 to Scientific calculator")
    print()
    var = int(input("Type an option: "))

    # Astronomical Calculator
    if var == 1:
        print()
        print("[1] Calculate luminosity and the habitable zone of a star")
        print("[2] Characteristics of a black hole based on its mass")
        print("[3] Convert light-year to parsec")
        print("[4] Convert parsec to light-year")
        print("[5] Convert meter to light-second")
        print("[6] Convert light-second to meter")
        print("[7] Convert light-year to meter")
        print("[8] Convert meter to light-year")
        print("[9] Calculate time dilation")
        print("[10] Convert mass to energy")
        print("[11] Convert energy in mass")
        print()
        opcao = int(input("Type an option: "))
        if opcao == 1:
            luminosity()
        elif opcao == 2:
            blackhole()
        elif opcao == 3:
            ly_to_pc()
        elif opcao == 4:
            pc_to_ly()
        elif opcao == 5:
            m_to_ls()
        elif opcao == 6:
            ls_to_m()
        elif opcao == 7:
            ly_to_m()
        elif opcao == 8:
            m_to_ly()
        elif opcao == 9:
            time_dilation()
        elif opcao == 10:
            mass_energy()
        elif opcao == 11:
            energy_mass()

        exit()

    # scientific calculator
    elif var == 2:
        print()
        print("[1] Greatest common divisor")
        print("[2] Square root")
        print("[3] Logarithm")
        print("[4] Power")
        print("[5] Complex multiplication")
        print("[6] Quaternion multiplication")
        print("[7] Quaternion division")
        print()
        opcao = int(input("Type an option: "))
        if opcao == 1:
            mdc()
        elif opcao == 2:
            sqrt()
        elif opcao == 3:
            log()
        elif opcao == 4:
            exponenciacao()
        elif opcao == 5:
            complex_multiplication()
        elif opcao == 6:
            quaternion_multiplication()
        elif opcao == 7:
            quaternion_division()

        exit()
# end of field functions

# function of beginning of the program
print()
start = input("Enter the calculator? y ou n: ")
if start == "y":
    main()
elif start == "n":
    sys.exit(0)
