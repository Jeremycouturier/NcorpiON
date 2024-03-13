########################################################################################
########################################################################################
########################################################################################
######### @file    ncorpion_animation.sh                                       #########
######### @brief   This bash script calls ffmpeg to make an animation          #########
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

cd $1gif &&
ffmpeg -y -i %d.png -vf palettegen palette.png &&
ffmpeg -framerate 30 -y -i "%d.png" -i palette.png -lavfi paletteuse ncorpion.gif &&
ffmpeg -framerate 30 -y -i "%d.png" -i palette.png -lavfi paletteuse ncorpion.mp4
