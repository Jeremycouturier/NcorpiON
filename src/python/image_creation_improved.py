########################################################################################
########################################################################################
########################################################################################
######### @file    image_creation_improved.py                                  #########
######### @brief   This script produces images of the simulation               #########
######### @author  Jérémy COUTURIER <jeremycouturier.com>                      #########
#########                                                                      #########
######### @section 	LICENSE                                                #########
######### Copyright (c) 2023 Jérémy COUTURIER                                  #########
#########                                                                      #########
######### This file is part of NcorpiON                                        #########
#########                                                                      #########
######### NcorpiON is free software. You can redistribute it and/or modify     #########
######### it under the terms of the GNU General Public License as published by #########
######### the Free Software Foundation, either version 3 of the License, or    #########
######### (at your option) any later version.                                  #########
#########                                                                      #########
######### NcorpiON is distributed in the hope that it will be useful,          #########
######### but WITHOUT ANY WARRANTY; without even the implied warranty of       #########
######### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #########
######### GNU General Public License for more details.                         #########
#########                                                                      #########
######### You should have received a copy of the GNU General Public License    #########
######### along with rebound.  If not, see <http://www.gnu.org/licenses/>.     #########
########################################################################################
########################################################################################
########################################################################################

###### This python script returns images that can be used to build an animated gif of the trajectories           ######
###### First run "ffmpeg -i %d.png -vf palettegen palette.png" to create a palette that is used to build the gif ######
###### Then run "ffmpeg -framerate 25 -i "%d.png" -i palette.png -lavfi paletteuse output.gif"  to builf the gif ######

import cmath as cm
import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import matplotlib.image as image
import sys

path           = str(sys.argv[4])   #The path indicated in the file src/parameters.h
inner_bl       = int(sys.argv[5])   #Boolean indicating whether or not the simulation featured an inner fluid disk
frag_bl        = int(sys.argv[6])   #Boolean indicating whether or not the simulation used the fragmentation model of NcorpiON
sideral_period = float(sys.argv[7]) #Sideral period of the central body. Needed to draw the correct figure of the central body.
surface_orbital_period_in_days = 0.0585745105636311 #Proportionality constant between the simulation's unit of time and 24 hours : surface orbital period/1 day.
                                                    #This value is for the Earth.

acosi_min = 0.0 #Left border of the image
acosi_max = 16.0 #Right border of the image
asini_min = 0.0 #Bottom border of the image
asini_max = 8.0 #Top border of the image

patha   = path + "a.txt"
pathe   = path + "e.txt"
pathi   = path + "i.txt"
pathrad = path + "radius.txt"
pathstat= path + "stat.txt"

a = open(patha,    "r")
e = open(pathe,    "r")
i = open(pathi,    "r")
R = open(pathrad,  "r")
S = open(pathstat, "r")


def f(e):
      return -1.0/9.0*(16.0*e**2-32.0*e+7.0)

def g(e):
      return 4.0*e*(1.0-e)
      
def h(e):
      return 64.0*e*(0.25-e)


def P2costheta(theta): #The second Legendre polynomial. Used to draw the correct figure of the central body
      return 0.5*(3.0*np.cos(theta)**2-1.0)


def draw_orbit(rad, xx, yy, linewidth, color, npoints): #draws a circle of radius rad at (xx, yy)
      tps=np.linspace(0.0,2.0*m.pi,npoints)
      x=rad*np.cos(tps)+xx
      y=rad*np.sin(tps)+yy
      py.plot(x,y,linestyle="-",color=color, linewidth=linewidth)

      
def draw_oblate_Earth(rad, sideral_period): #The sideral period is given in units of the surface orbital period
      eps_20 = -5.0/6.0/(sideral_period**2)
      tps=np.linspace(0.0,2.0*m.pi,250)
      x=rad*np.sin(tps)*(1.0+eps_20*P2costheta(tps))
      y=rad*np.cos(tps)*(1.0+eps_20*P2costheta(tps))
      py.plot(x,y,linestyle="-",color="grey", linewidth=4.0)


def artisanal_colorbar(bottom,top,left,right):

      ecc = np.linspace(0.0,1.0,500)
      where_to_plot = np.linspace(bottom,top,500)
      for p in range(500):
            if (p < 125):
                  color = (0.0,g(ecc[p]),h(ecc[p]))
            else:
                  color = (f(ecc[p]),g(ecc[p]),0.0)
            py.plot([left,right], [where_to_plot[p],where_to_plot[p]], '-', color=color, linewidth=1)
      py.text(acosi_max*0.915,0.5*(bottom+top)-0.07,r"$e$", color="black", fontsize=25, alpha = 0.5)
      py.plot([left,left],[bottom,top],'-',color='black',linewidth=1)
      py.plot([right,right],[bottom,top],'-',color='black',linewidth=1)
      py.plot([left,right],[bottom,bottom],'-',color='black',linewidth=1)
      py.plot([left,right],[top,top],'-',color='black',linewidth=1)
      py.text(acosi_max*0.963,asini_min+0.038*(asini_max-asini_min),r"$0$", fontsize=15, alpha = 0.5)
      py.text(acosi_max*0.963,asini_max-0.06 *(asini_max-asini_min),r"$1$", fontsize=15, alpha = 0.5)
      py.text(acosi_max*0.963,0.5*(bottom+top)-0.07,r"$1/2$", fontsize=15, alpha = 0.5)
      for d in range(1,10):
            py.plot([left,right], [where_to_plot[50*d],where_to_plot[50*d]], '-', color="black", linewidth=0.5, alpha = 0.4)

      
def draw_moonlet(sma, ecc, inc, rad, p, maxR1, maxR2, maxR3, largest_index_1, largest_index_2, largest_index_3):
      radius = rad[p]
      if (radius > maxR1):
            maxR3 = maxR2
            maxR2 = maxR1
            maxR1 = radius
            largest_index_3 = largest_index_2
            largest_index_2 = largest_index_1
            largest_index_1 = p
      elif (radius > maxR2):
            maxR3 = maxR2
            maxR2 = radius
            largest_index_3 = largest_index_2
            largest_index_2 = p
      elif (radius > maxR3):
            maxR3 = radius
            largest_index_3 = p
      eccentricity = ecc[p]
      if (eccentricity < 0.25):
            color = (0.0, g(eccentricity), h(eccentricity))
      elif (eccentricity <= 1.0):
            color = (f(eccentricity), g(eccentricity), 0.0)
      else:
            color = 'yellow'
      X = sma[p] * m.cos(inc[p])
      Y = sma[p] * m.sin(inc[p])
      if (X <= acosi_max + 0.2 and Y <= asini_max + 0.2):
            if (radius >= 0.12):
                  draw_orbit(radius, X, Y, 3, color, 150)
            elif (radius >= 0.06):
                  draw_orbit(radius, X, Y, 2, color, 100)
            elif (radius >= 0.02):
                  draw_orbit(radius, X, Y, 1, color, 20)
            elif (radius >= 0.001):
                  draw_orbit(radius, X, Y, 1, color, 10)
            else:
                  py.plot(X, Y, marker=",", color = color)
      return (maxR1, maxR2, maxR3, largest_index_1, largest_index_2, largest_index_3)

      
def make_image(sma, ecc, inc, rad, sta, q):
      py.xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
      draw_oblate_Earth(1.0,sideral_period)
      artisanal_colorbar(asini_min+0.05*(asini_max-asini_min),asini_max-0.05*(asini_max-asini_min),acosi_max*0.94,acosi_max*0.96)
      maxR1 = 0.0
      maxR2 = 0.0
      maxR3 = 0.0
      largest_index_1 = 0
      largest_index_2 = 0
      largest_index_3 = 0
      for p in range(len(sma)):
            (maxR1,maxR2,maxR3,largest_index_1,largest_index_2,largest_index_3) = draw_moonlet(sma,ecc,inc,rad,p,maxR1,maxR2,maxR3,largest_index_1,largest_index_2,largest_index_3)
      draw_moonlet(sma, ecc, inc, rad, largest_index_1, 0.0, 0.0, 0.0, 0, 0, 0)
      draw_moonlet(sma, ecc, inc, rad, largest_index_2, 0.0, 0.0, 0.0, 0, 0, 0)
      draw_moonlet(sma, ecc, inc, rad, largest_index_3, 0.0, 0.0, 0.0, 0, 0, 0)
      py.text(0.08,asini_max*0.94,"t = " + str(round(float(sta[0]),3)) + " = " + str(round(float(sta[0])*surface_orbital_period_in_days,3)) + " days", alpha=0.4, fontsize=20)
      py.text(0.08,asini_max*0.88,"N = " + sta[1], alpha=0.4, fontsize=20)
      #if (frag_bl and inner_bl):
      #      py.text(0.08,asini_max*0.82,"Collisions = " + sta[2] + ",  (m, sc, hf, ff) = (" + sta[6]+", "+ sta[7]+", "+ sta[8]+", "+ sta[9]+")", alpha=0.4, fontsize=20)
      #elif (frag_bl):
      #      py.text(0.08,asini_max*0.82,"Collisions = " + sta[2] + ",  (m, sc, hf, ff) = (" + sta[5]+", "+ sta[6]+", "+ sta[7]+", "+ sta[8]+")", alpha=0.4, fontsize=20)
      #else:
      #      py.text(0.08,asini_max*0.82,"Collisions = " + sta[2], alpha=0.4, fontsize=20)
      py.text(0.08,asini_max*0.82,"Collisions = " + sta[2], alpha=0.4, fontsize=20)
      py.text(0.08,asini_max*0.76,"Largest radii = (" + str(round(maxR1,3)) + ", " + str(round(maxR2,3)) + ", " + str(round(maxR3,3)) + ")" + r" $R_{\oplus}$", alpha=0.4, fontsize=20)
      if (inner_bl):
            py.text(0.08,asini_max*0.70,"(Moonlet mass, Inner disk mass) = (" + str(round(float(sta[4]),4))+", "+str(round(float(sta[5]),4))+")"+r" $M_{\oplus}$", alpha=0.4, fontsize=20)
      else:
            py.text(0.08,asini_max*0.70,"Moonlet mass = " + str(round(float(sta[4]),4)) + r" $M_{\oplus}$", alpha=0.4, fontsize=20)
      py.text(0.08,asini_max*0.65,"thin: R < 0.06",  alpha=0.4, fontsize=8)
      py.text(0.08,asini_max*0.62,"thick: R > 0.06",  alpha=0.4, fontsize=8)
      py.text(0.08,asini_max*0.59,"thickest: R > 0.12",  alpha=0.4, fontsize=8)
      py.text(0.08,asini_max*0.56,"single pixel (barely visible): R < 0.001", alpha=0.4, fontsize=8)
      py.text(0.08,asini_max*0.53,"Everything is at scale",  alpha=0.4, fontsize=8)
      py.axis('square')
      py.xlim([acosi_min,acosi_max])
      py.ylim([asini_min,asini_max])
      py.xticks(fontsize = 15)
      py.yticks(fontsize = 15)
      py.xlabel(r"$a\,\cos\,i$ (Earth radii)", fontsize = 25)
      py.ylabel(r"$a\,\sin\,i$ (Earth radii)", fontsize = 25)
      figure = py.gcf() ##get current figure
      figure.set_size_inches(16, 9)
      #py.axis('off')
      py.savefig(path + "gif/" + str(q) + '.png', dpi = 160)
      py.clf()


mod      = int(sys.argv[1])  # Only treats images equal to mod modulo n_thread
n_thread = int(sys.argv[2])  # Number of threads (to be modified in image_creation.sh)

sma = np.float64(np.array(a.readline().strip().split()))
ecc = np.float64(np.array(e.readline().strip().split()))
inc = np.float64(np.array(i.readline().strip().split()))
rad = np.float64(np.array(R.readline().strip().split()))
sta = np.array(S.readline().strip().split())

for j in range(int(sys.argv[3])):
      if (j % n_thread == mod):
            print ("Producing image n° ", j)
            make_image(sma, ecc, inc, rad, sta, j)
      sma = np.float64(np.array(a.readline().strip().split()))
      ecc = np.float64(np.array(e.readline().strip().split()))
      inc = np.float64(np.array(i.readline().strip().split()))
      rad = np.float64(np.array(R.readline().strip().split()))
      sta = np.array(S.readline().strip().split())
           
     
            
