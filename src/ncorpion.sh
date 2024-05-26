########################################################################################
########################################################################################
########################################################################################
######### @file    ncorpion.sh                                                 #########
######### @brief   This bash script compiles and runs NcorpiON                 #########
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

make clean && make && make clean

read foo <makefile

visualization="${foo: -1}"

if [ $visualization = 1 ]; then
      cd 3D
      if [ ! -f ./librebound.so ]; then
            make clean
      fi
      make && mpirun -n 1 ../ncorpion : -n 1 ./rebound && cd ../
else
      ./ncorpion
fi
