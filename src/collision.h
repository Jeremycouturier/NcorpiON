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
/******** @file    collision.h                                                 ********/
/******** @brief   Header file to collision.c                                  ********/
/******** @author  Jérémy COUTURIER <jeremycouturier.com>                      ********/
/********                                                                      ********/
/******** @section  LICENSE                                                    ********/
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


#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "parameters.h"
#include "structure.h"
#include "ffm.h"


typ * closest_approach(struct moonlet * moonlets, int a, int b);


void hash_table_cell_index(struct moonlet * moonlets, int a);


void hash_table_cell_index_z_order(struct moonlet * moonlets, int a);


void neighbours(struct moonlet * moonlets, int a);


void mesh(struct moonlet * moonlets);


void brute_force(struct moonlet * moonlets);


void get_center_and_maxR(struct node * FlatTree, struct moonlet * moonlets, int a);


void get_rmax_and_rcrit(struct node * FlatTree, struct moonlet * moonlets, int a);


void get_center_and_maxR_from_children(struct node * FlatTree, int a);


void get_rmax_and_rcrit_from_children(struct node * FlatTree, int a);


void center_and_maxR_flattree(struct node * FlatTree, struct moonlet * moonlets);


void rmax_and_rcrit_flattree(struct node * FlatTree, struct moonlet * moonlets);


void collision_flattree(struct node * FlatTree, struct moonlet * moonlets);


void standard_tree_collision(struct node * FlatTree, struct moonlet * moonlets, int b);


#endif
