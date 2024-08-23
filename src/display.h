/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    display.h                                                   ********/
/******** @brief   Header file to display.c                                    ********/
/******** @author  Jérémy COUTURIER <jeremycouturier.com>                      ********/
/********                                                                      ********/
/******** @section 	LICENSE                                                ********/
/******** Copyright (c) 2023 Jérémy COUTURIER                                  ********/
/********                                                                      ********/
/******** This file is part of NcorpiON                                        ********/
/********                                                                      ********/
/******** NcorpiON is free software. You can redistribute it and/or modify     ********/
/******** it under the terms of the GNU General Public License as published by ********/
/******** the Free Software Foundation, either version 3 of the License, or    ********/
/******** (at your option) any later version.                                  ********/
/********                                                                      ********/
/******** NcorpiON is distributed in the hope that it will be useful,          ********/
/******** but WITHOUT ANY WARRANTY; without even the implied warranty of       ********/
/******** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ********/
/******** GNU General Public License for more details.                         ********/
/********                                                                      ********/
/******** You should have received a copy of the GNU General Public License    ********/
/******** along with NcorpiON.  If not, see <http://www.gnu.org/licenses/>.    ********/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/


#ifndef _DISPLAY_H_
#define _DISPLAY_H_

#include "parameters.h"


/******** Defining a structure holding the data of a collision. Used by function display to output collision data to files ********/
struct collisionData{
      typ time;         //Simulation time when the collision occured
      typ m1;           //Mass   of the impactor
      typ m2;           //Mass   of the target
      typ R1;           //Radius of the impactor
      typ R2;           //Radius of the target
      typ DeltaV;       //Relative velocity
      typ impact_angle; //Impact angle
      typ m_tilde;      //Mass of the largest remnant
};
extern struct collisionData * collisionDatas; //Array containing the datas of the collisions
extern int indexCollision;                    //Index in the array


FILE ** file_opening();


void file_closing();


void display(struct moonlet * moonlets);


void readme();


void resumeFile();


void rebound(struct moonlet * moonlets);


void resume(struct moonlet * moonlets);


#endif
