/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    rk4.h                                                       ********/
/******** @brief   Header file to rk4.c                                        ********/
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


#ifndef _RK4_H_
#define _RK4_H_

#include "parameters.h"
#include "structure.h"
#include "ffm.h"


extern struct node * FlatTree;


void kick(struct moonlet * X, struct moonlet * C, void (*F)(struct moonlet *));


void drift(struct moonlet * X, struct moonlet * C);


void end_of_timestep(struct moonlet * moonlets, int progressed);


void integration_tree(typ t);


void integration_mesh(typ t);


void integration_brute_force_SABA1(typ t);


#endif
