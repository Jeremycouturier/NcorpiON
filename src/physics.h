/*
       NNNNNNNN        NNNNNNNN        CCCCCCCCCCCCC     OOOOOOOOO     RRRRRRRRRRRRRRRRR   PPPPPPPPPPPPPPPPP   IIIIIIIIII     OOOOOOOOO     NNNNNNNN        NNNNNNNN
       N:::::::N       N::::::N     CCC::::::::::::C   OO:::::::::OO   R::::::::::::::::R  P::::::::::::::::P  I::::::::I   OO:::::::::OO   N:::::::N       N::::::N
       N::::::::N      N::::::N   CC:::::::::::::::C OO:::::::::::::OO R::::::RRRRRR:::::R P::::::PPPPPP:::::P I::::::::I OO:::::::::::::OO N::::::::N      N::::::N
       N:::::::::N     N::::::N  C:::::CCCCCCCC::::CO:::::::OOO:::::::ORR:::::R     R:::::RPP:::::P     P:::::PII::::::IIO:::::::OOO:::::::ON:::::::::N     N::::::N
       N::::::::::N    N::::::N C:::::C       CCCCCCO::::::O   O::::::O  R::::R     R:::::R  P::::P     P:::::P  I::::I  O::::::O   O::::::ON::::::::::N    N::::::N
       N:::::::::::N   N::::::NC:::::C              O:::::O     O:::::O  R::::R     R:::::R  P::::P     P:::::P  I::::I  O:::::O     O:::::ON:::::::::::N   N::::::N
       N:::::::N::::N  N::::::NC:::::C              O:::::O     O:::::O  R::::RRRRRR:::::R   P::::PPPPPP:::::P   I::::I  O:::::O     O:::::ON:::::::N::::N  N::::::N
       N::::::N N::::N N::::::NC:::::C              O:::::O     O:::::O  R:::::::::::::RR    P:::::::::::::PP    I::::I  O:::::O     O:::::ON::::::N N::::N N::::::N
       N::::::N  N::::N:::::::NC:::::C              O:::::O     O:::::O  R::::RRRRRR:::::R   P::::PPPPPPPPP      I::::I  O:::::O     O:::::ON::::::N  N::::N:::::::N
       N::::::N   N:::::::::::NC:::::C              O:::::O     O:::::O  R::::R     R:::::R  P::::P              I::::I  O:::::O     O:::::ON::::::N   N:::::::::::N
       N::::::N    N::::::::::NC:::::C              O:::::O     O:::::O  R::::R     R:::::R  P::::P              I::::I  O:::::O     O:::::ON::::::N    N::::::::::N
       N::::::N     N:::::::::N C:::::C       CCCCCCO::::::O   O::::::O  R::::R     R:::::R  P::::P              I::::I  O::::::O   O::::::ON::::::N     N:::::::::N
       N::::::N      N::::::::N  C:::::CCCCCCCC::::CO:::::::OOO:::::::ORR:::::R     R:::::RPP::::::PP          II::::::IIO:::::::OOO:::::::ON::::::N      N::::::::N
       N::::::N       N:::::::N   CC:::::::::::::::C OO:::::::::::::OO R::::::R     R:::::RP::::::::P          I::::::::I OO:::::::::::::OO N::::::N       N:::::::N
       N::::::N        N::::::N     CCC::::::::::::C   OO:::::::::OO   R::::::R     R:::::RP::::::::P          I::::::::I   OO:::::::::OO   N::::::N        N::::::N
       NNNNNNNN         NNNNNNN        CCCCCCCCCCCCC     OOOOOOOOO     RRRRRRRR     RRRRRRRPPPPPPPPPP          IIIIIIIIII     OOOOOOOOO     NNNNNNNN         NNNNNNN
*/

/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    physics.h                                                   ********/
/******** @brief   Header file to physics.c                                    ********/
/******** @author  Jérémy COUTURIER <jeremycouturier.com>                      ********/
/********                                                                      ********/
/******** @section 	LICENSE                                                    ********/
/******** Copyright (c) 2023 Jérémy COUTURIER                                  ********/
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


#ifndef _PHYSICS_H_
#define _PHYSICS_H_


#include "parameters.h"
#include "structure.h"


void vector_field(struct moonlet * moonlets);


void KelvinVoigtDamping(struct moonlet * X);


void tides(struct moonlet * bodies);


void collision(struct moonlet * moonlets, int a, int b, typ f);


void merger(struct moonlet * moonlets, int a, int b);


void fragmentation(struct moonlet * moonlets, int a, int b);


void collision_treatment(struct moonlet * moonlets, int a, int b, int type_of_collision);


void get_neighbours_mesh(struct moonlet * moonlets);


int willMergeWithDisk(struct moonlet * bodies, int id);


void innerFluidDiskAngularMomentum(typ m, typ a, typ e, typ cosi, typ g1, typ m1);


int evectionResonance(typ * X, int inner);

#endif
