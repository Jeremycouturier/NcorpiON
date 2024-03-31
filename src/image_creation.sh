########################################################################################
########################################################################################
########################################################################################
######### @file    image_creation.sh                                           #########
######### @brief   This bash script calls python to produce images             #########
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

n_thread=8
mkdir -p $2gif &&
python3 python/image_creation_improved.py 0 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 1 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 2 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 3 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 4 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 5 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 6 $n_thread $1 $2 $3 $4 $5 $6 &
python3 python/image_creation_improved.py 7 $n_thread $1 $2 $3 $4 $5 $6 &
wait
