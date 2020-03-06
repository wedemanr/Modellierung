#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Berueberstoemung
================
Darstellung der analytischen Loesung nach Long (1953) und Lilly und Klemp (1979).
Dieses Programm ist im Rahmen der Uebung zur Lehrveranstaltung "Meteorologische
Modellierung I" 2013/2014 entstanden.

Lizenz
------
Copyright 2013 Marek Jacob, marek.jacob@zmaw.de

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

########
## 
## Parameter
## 
L_zb = 300 # m Scheitelhoehe
L_xb = 3000 # m Halbwertsbreite
theta_0 = 290 # K # Bodentemperatur
gamma = 0.005 # K/m Schichtung
g = 9.81 # m/s^2 
U = 10. # m/s Grossskaligert Wind

N = np.sqrt(g/theta_0  * gamma)	


########
##
## Hilfsfunktionen
##
def z_s(x):
	"""Berghoehe nach A-1"""
	return L_zb * np.power(L_xb,2) / (np.power(L_xb,2) + np.power(x,2))

def f_2(x):
	return -(x/L_xb) * z_s(x)


########
##
## Felder der Lösung im Form vom Feldfunktionen
##
def delta(x,z):
	first_term = z_s(x) * np.cos((N/U) * (z-z_s(x)))
	second_term = f_2(x) * np.sin((N/U) * (z-z_s(x)))
	return  first_term + second_term 

def theta(x,z):
	return theta_0 * (1+N*N * (z-delta(x,z)) / g)

def DdeltaDz(x,z):
	"""analytische Ableitung von delta nach z"""
	first_term = -(N/U) * z_s(x) * np.sin((N/U) * (z-z_s(x)))
	second_term = (N/U) * f_2(x) * np.cos((N/U) * (z-z_s(x)))
	return  first_term + second_term

def u(x,z):
	return U * (1 - DdeltaDz(x,z) )

def w(x,z):
	return U * DdeltaDz(x,z)

# Ein PDF fuer alle Plots
pdf = matplotlib.backends.backend_pdf.PdfPages('bergueberstroemung.pdf')


########
##
## Zeichnen des Berges
##
fig = plt.figure(figsize=(4, 4), dpi=110)
x=np.linspace(-5000,5000,100) # m, -5km ... 5km in 100 schritten
ax = fig.add_subplot(1,1,1)
ax.plot(x,z_s(x)) # x gegen Berghoehe
pdf.savefig(fig)
plt.close(fig)

# Variablenfelder X und Z
x = np.linspace(-80e3, 80e3, 100) # m, -80km ... 80km in 100 schritten
z = np.linspace(0,10e3,100) # m, Höhe 0km ... 10km in 100 Schritten
X,Z=np.meshgrid(x,z)


########
##
## Zeichnen der Felder
##
plots = {
	('delta',delta),
	('mean u',u),
	('mean w',w),
	('Theta [K]',theta),
}
for title, function in plots:
	fig = plt.figure(figsize=(10, 6), dpi=110)
	ax = fig.add_subplot(1,1,1)

	# bunte Konturflaechen
	CS = ax.contourf(X,Z,function(X,Z),40,zorder=0)
	cbar = fig.colorbar(CS) # Colorbar
	# schwarze Konturlinien
	ax.contour(X,Z,function(X,Z),10,colors='k',zorder=1)

	#Berg einmalen.
	ax.fill(x,z_s(x),'w',zorder=2)

	#label
	ax.set_xlabel("$x$")
	ax.set_ylabel("$z$")
	ax.set_title(title)
	
	pdf.savefig(fig)
	plt.close(fig)


pdf.close()
