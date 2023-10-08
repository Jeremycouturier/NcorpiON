###### This python script returns images that can be used to build an animated gif of the trajectories           ######
###### First run "ffmpeg -i %d.png -vf palettegen palette.png" to create a palette that is used to build the gif ######
###### Then run "ffmpeg -framerate 25 -i "%d.png" -i palette.png -lavfi paletteuse output.gif"  to builf the gif ######

sideral_period = 3.4076 #To draw the correct figure of the Earth. Unit of time is the surface orbital period

import cmath as cm
import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import matplotlib.image as image
import sys

path = str(sys.argv[3])

a_min = 0.0 #Left border of the gif
a_max = 16.0 #Right border of the gif
asini_min = 0.0 #Bottom border of the gif
asini_max = 8.0 #Top border of the gif

patha=path+"a.txt"
pathe=path+"e.txt"
pathi=path+"i.txt"
pathrad=path+"radius.txt"
pathstat=path+"stat.txt"

a = open(patha,    "r")
e = open(pathe,    "r")
i = open(pathi,    "r")
R = open(pathrad,  "r")
S = open(pathstat, "r")


def g(e):
      return 4.0*e*(1.0-e)


def P2costheta(theta):
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
            if (p <= 249):
                  color = (0.0,0.0,g(ecc[p]))
            else:
                  color = (g(ecc[p]-0.5),0.0,g(ecc[p]))
            py.plot([left,right], [where_to_plot[p],where_to_plot[p]], '-', color=color, linewidth=1)
      py.text(a_max*0.915,0.5*(bottom+top),r"$e$", color="black", fontsize=25)
      py.plot([left,left],[bottom,top],'-',color='black',linewidth=1)
      py.plot([right,right],[bottom,top],'-',color='black',linewidth=1)
      py.plot([left,right],[bottom,bottom],'-',color='black',linewidth=1)
      py.plot([left,right],[top,top],'-',color='black',linewidth=1)
      py.text(a_max*0.963,asini_min+0.04*(asini_max-asini_min),r"$0$", fontsize=15)
      py.text(a_max*0.963,asini_max-0.06*(asini_max-asini_min),r"$1$", fontsize=15)
      py.text(a_max*0.963,0.5*(bottom+top),r"$1/2$", fontsize=15)

      
def draw_moonlet(sma, ecc, inc, rad, p, maxR, largest_index):
      radius = rad[p]
      if (radius > maxR):
            maxR = radius
            largest_index = p
      eccentricity = ecc[p]
      if (eccentricity < 0.5):
            color = (0.0, 0.0, g(eccentricity))
      elif (eccentricity <= 1.0):
            color = (g(eccentricity-0.5), 0.0, g(eccentricity))
      else:
            color = 'green'
      X = sma[p] * m.cos(inc[p])
      Y = sma[p] * m.sin(inc[p])
      if (X <= a_max+0.2 and Y <= asini_max+0.2):
            if (radius >= 0.08):
                  draw_orbit(radius, X, Y, 1.5, color, 100)
            elif (radius >= 0.02):
                  draw_orbit(radius, X, Y, 1, color, 20)
            elif (radius >= 0.001):
                  draw_orbit(radius, X, Y, 1, color, 10)
            else:
                  py.plot(X, Y, marker=",", color=color)
      return (maxR, largest_index)

      
def make_image(sma, ecc, inc, rad, sta, q):
      draw_oblate_Earth(1.0,sideral_period)
      maxR = 0.0
      largest_index = 0
      for p in range(len(sma)):
            (maxR, largest_index) = draw_moonlet(sma, ecc, inc, rad, p, maxR, largest_index)
      draw_moonlet(sma, ecc, inc, rad, largest_index, 0.0, 0)
      py.text(0.08,asini_max*0.94,"t = "+str(round(float(sta[0]),3)), fontsize=20)
      py.text(0.08,asini_max*0.88,"N = "+sta[1], fontsize=20)
      py.text(0.08,asini_max*0.82,"Collisions = "+sta[2], fontsize=20)
      py.text(0.08,asini_max*0.76,"Largest moonlet = "+str(round(float(sta[3]),3))+r" $R_{\oplus}$", fontsize=20)
      py.text(0.08,asini_max*0.70,"Orbiting mass = "+str(round(float(sta[4]),4))+r" $M_{\oplus}$", fontsize=20)
      artisanal_colorbar(asini_min+0.05*(asini_max-asini_min),asini_max-0.05*(asini_max-asini_min),a_max*0.94,a_max*0.96)
      py.axis('square')
      py.xlim([a_min,a_max])
      py.ylim([asini_min,asini_max])
      py.xticks(fontsize=15)
      py.yticks(fontsize=15)
      py.xlabel(r"$a\,\cos i$ (Earth radii)", fontsize=25)
      py.ylabel(r"$a\,\sin i$ (Earth radii)", fontsize=25)
      figure = py.gcf() ##get current figure
      figure.set_size_inches(16, 9)
      #py.axis('off')
      py.savefig(path+"gif/"+str(q)+'.png', dpi=160)
      py.clf()


how_many_to_skip = int(sys.argv[1]) #Number of images to skip
how_many_to_make = int(sys.argv[2]) #Number of images to make

for j in range(how_many_to_skip+1):      
      sma = np.float64(np.array(a.readline().strip().split()))
      ecc = np.float64(np.array(e.readline().strip().split()))
      inc = np.float64(np.array(i.readline().strip().split()))
      rad = np.float64(np.array(R.readline().strip().split()))
      sta = np.array(S.readline().strip().split())
      
counter = how_many_to_skip

for j in range(how_many_to_make):
      print("t = "+str(float(sta[0])))
      make_image(sma, ecc, inc, rad, sta, counter)
      sma = np.float64(np.array(a.readline().strip().split()))
      ecc = np.float64(np.array(e.readline().strip().split()))
      inc = np.float64(np.array(i.readline().strip().split()))
      rad = np.float64(np.array(R.readline().strip().split()))
      sta = np.array(S.readline().strip().split())
      counter += 1

                  
     
            
