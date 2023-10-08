#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "parameters.h"
#include "structure.h"
#include "collision.h"
#include "physics.h"
#include "ffm.h"
#include "rk4.h"
#include <errno.h>
#include <math.h>



void vector_field(struct moonlet * moonlets){

      /******** The vector field passed as argument to the integrator         ********/
      /******** Stores the accelerations into the velocity fields of moonlets ********/
      
      
      int k,p;
      typ rk,r2,r3,r5;
      typ X,Y,Z,aX,aY,aZ,K;
      typ mu,mk;
      typ Xp,Yp,Zp,mp,DX,DY,DZ,D,D3,rp,rp3,softening,Rk,Rp;
      typ X0,Y0,Z0,m0,R0,X1,Y1,Z1,m1,R1,X2,Y2,Z2,m2,R2;
      
      
      /******** Getting the coordinates of the three largest moonlets. Useful only when treating mutual gravitational interactions with the O(N) mesh algorithm ********/
      if (mutual_bool && mesh_bool && !force_naive_bool){
            X0 =  moonlets -> x;
            Y0 =  moonlets -> y;
            Z0 =  moonlets -> z;
            m0 =  moonlets -> mass;
            R0 =  moonlets -> radius;
            X1 =  (moonlets+1) -> x;
            Y1 =  (moonlets+1) -> y;
            Z1 =  (moonlets+1) -> z;
            m1 =  (moonlets+1) -> mass;
            R1 =  (moonlets+1) -> radius;
            X2 =  (moonlets+2) -> x;
            Y2 =  (moonlets+2) -> y;
            Z2 =  (moonlets+2) -> z;
            m2 =  (moonlets+2) -> mass;
            R2 =  (moonlets+2) -> radius;
      }
      
      for (k = 0; k <= largest_id; k++){
            if (*(exists+k)){ //Checking whether or not there is a moonlet in the kth cell of the moonlet array
            
                  /******** Contribution from the Earth's center of mass ********/
                  X  = (moonlets+k)-> x;
                  Y  = (moonlets+k)-> y;
                  Z  = (moonlets+k)-> z;
                  mk = (moonlets+k)->mass;
                  rk = sqrt(X*X+Y*Y+Z*Z);
                  r3 = rk*rk*rk;
                  
                  mu = G*Mearth;
                  aX = -mu*X/r3;
                  aY = -mu*Y/r3;
                  aZ = -mu*Z/r3;
                  
                  
                  /******** Contribution from the Earth symmetrical equatorial bulge ********/
                  if (J2_bool){
                        r2 = rk*rk;
                        r5 = r3*r2;
                        K  = G*Mearth*Rearth*Rearth*J2/r5;
                        /******** Updating acceleration of moonlet k ********/
                        aX += K*(-1.5*X+7.5/r2*Z*Z*X);
                        aY += K*(-1.5*Y+7.5/r2*Z*Z*Y);
                        aZ += K*(-4.5*Z+7.5/r2*Z*Z*Z);
                  }
                  
                  
                  /******** Contribution from the Sun ********/
                  if (Sun_bool){
                  
                  }
                  
                  
                  /******** Contribution from the inner fluid disk ********/
                  if (disk_bool){
                  
                  }
                  
                  /******** Mutual gravitational interactions with the brute-force O(N^2) algorithm ********/
                  if (mutual_bool && (brute_force_bool || force_naive_bool)){
                        for (p=0; p < k; p++){
                              if (*(exists+p)){
                              
                                    /******** Getting the positions, masses and radii ********/
                                    mp = (moonlets+p)-> mass;
                                    Xp = (moonlets+p)-> x;
                                    Yp = (moonlets+p)-> y;
                                    Zp = (moonlets+p)-> z;
                                    Rp = (moonlets+p)-> radius;
                                    Rk = (moonlets+k)-> radius;
                                    
                                    DX = X-Xp;
                                    DY = Y-Yp;
                                    DZ = Z-Zp;
                                    softening = softening_parameter*(Rk + Rp);
                                    D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                                    D3 = D*D*D;
                                    rp = sqrt(Xp*Xp+Yp*Yp+Zp*Zp);
                                    rp3= rp*rp*rp;

                                    /******** Updating acceleration of moonlet k ********/
                                    aX -= G*mp*DX/D3;
                                    aY -= G*mp*DY/D3;
                                    aZ -= G*mp*DZ/D3;
                                    
                                    /******** Updating acceleration of moonlet p ********/
                                    (moonlets+p)->vx  += G*mk*DX/D3; //dV/dt=A
                                    (moonlets+p)->vy  += G*mk*DY/D3;
                                    (moonlets+p)->vz  += G*mk*DZ/D3;
                              }
                        }
                  }
                  
                  /******** Mutual gravitational interactions with the mesh O(N) algorithm ********/
                  if (mutual_bool && mesh_bool && !force_naive_bool){
                  
                        /******** For now, we only consider gravitational interactions between pairs containing one of the three largest moonlets and another moonlet ********/
                        if (k > 2){
                              Rk = (moonlets+k)-> radius;
                              
                              /******** pair (0,k) ********/
                              if (exists[0]){
                                    DX = X-X0;
                                    DY = Y-Y0;
                                    DZ = Z-Z0;
                                    softening = softening_parameter*(Rk + R0);
                                    D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                                    D3 = D*D*D;
                                    /******** Updating acceleration of moonlet k ********/
                                    aX -= G*m0*DX/D3;
                                    aY -= G*m0*DY/D3;
                                    aZ -= G*m0*DZ/D3;
                                    /******** Updating acceleration of moonlet 0 ********/
                                    moonlets -> vx += G*mk*DX/D3;
                                    moonlets -> vy += G*mk*DY/D3;
                                    moonlets -> vz += G*mk*DZ/D3;
                              }
                              
                              /******** pair (1,k) ********/
                              if (exists[1]){
                                    DX = X-X1;
                                    DY = Y-Y1;
                                    DZ = Z-Z1;
                                    softening = softening_parameter*(Rk + R1);
                                    D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                                    D3 = D*D*D;
                                    /******** Updating acceleration of moonlet k ********/
                                    aX -= G*m1*DX/D3;
                                    aY -= G*m1*DY/D3;
                                    aZ -= G*m1*DZ/D3;
                                    /******** Updating acceleration of moonlet 1 ********/
                                    (moonlets+1) -> vx += G*mk*DX/D3;
                                    (moonlets+1) -> vy += G*mk*DY/D3;
                                    (moonlets+1) -> vz += G*mk*DZ/D3;
                              }
                              
                              /******** pair (2,k) ********/
                              if (exists[2]){
                                    DX = X-X2;
                                    DY = Y-Y2;
                                    DZ = Z-Z2;
                                    softening = softening_parameter*(Rk + R2);
                                    D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                                    D3 = D*D*D;
                                    /******** Updating acceleration of moonlet k ********/
                                    aX -= G*m2*DX/D3;
                                    aY -= G*m2*DY/D3;
                                    aZ -= G*m2*DZ/D3;
                                    /******** Updating acceleration of moonlet 2 ********/
                                    (moonlets+2) -> vx += G*mk*DX/D3;
                                    (moonlets+2) -> vy += G*mk*DY/D3;
                                    (moonlets+2) -> vz += G*mk*DZ/D3;
                              }
                        }      
                  }
                  
                  /******** Actualizing the acceleration of moonlet k ********/
                  (moonlets+k)->vx = aX; //dV/dt=A
                  (moonlets+k)->vy = aY;
                  (moonlets+k)->vz = aZ;
                  
            }     
      }
      
      
      /******** Mutual gravitational interactions with a tree algorithm (falcON or the standard tree code) ********/
      if (mutual_bool && (falcON_bool || standard_tree_bool) && !force_naive_bool){
            if (falcON_bool){
                  Cm_flattree(FlatTree, moonlets);
                  Cm_downtree(FlatTree, moonlets);
            }
            else if (standard_tree_bool){
                  for (k = 0; k <= largest_id; k++){
                        if (exists[k]){
                              standard_tree_acceleration(FlatTree, moonlets, k);
                        }
                  }
            }
            for (k = 0; k <= largest_id; k++){
                  if (exists[k]){
                        (moonlets+k) -> vx += C1Moonlets[3*k]  ;
                        (moonlets+k) -> vy += C1Moonlets[3*k+1];
                        (moonlets+k) -> vz += C1Moonlets[3*k+2];
                  }
            }
      }

      
      /******** Mutual gravitational interactions with the mesh O(N) algorithm ********/
      if (mutual_bool && mesh_bool && !brute_force_bool && !force_naive_bool){
      
            /******** Gravitational interactions between the three largest moonlets ********/
            /******** pair (0,1) ********/
            if (exists[0] && exists[1]){
                  DX = X0-X1;
                  DY = Y0-Y1;
                  DZ = Z0-Z1;
                  softening = softening_parameter*(R0 + R1);
                  D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                  D3 = D*D*D;
                  moonlets->vx  -= G*m1*DX/D3;
                  moonlets->vy  -= G*m1*DY/D3;
                  moonlets->vz  -= G*m1*DZ/D3;
                  (moonlets+1)->vx  += G*m0*DX/D3;
                  (moonlets+1)->vy  += G*m0*DY/D3;
                  (moonlets+1)->vz  += G*m0*DZ/D3;
            }
            
            /******** pair (0,2) ********/
            if (exists[0] && exists[2]){
                  DX = X0-X2;
                  DY = Y0-Y2;
                  DZ = Z0-Z2;
                  softening = softening_parameter*(R0 + R2);
                  D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                  D3 = D*D*D;
                  moonlets->vx  -= G*m2*DX/D3;
                  moonlets->vy  -= G*m2*DY/D3;
                  moonlets->vz  -= G*m2*DZ/D3;
                  (moonlets+2)->vx  += G*m0*DX/D3;
                  (moonlets+2)->vy  += G*m0*DY/D3;
                  (moonlets+2)->vz  += G*m0*DZ/D3;
            }
            
            /******** pair (1,2) ********/
            if (exists[1] && exists[2]){
                  DX = X1-X2;
                  DY = Y1-Y2;
                  DZ = Z1-Z2;
                  softening = softening_parameter*(R1 + R2);
                  D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                  D3 = D*D*D;
                  (moonlets+1)->vx  -= G*m2*DX/D3;
                  (moonlets+1)->vy  -= G*m2*DY/D3;
                  (moonlets+1)->vz  -= G*m2*DZ/D3;
                  (moonlets+2)->vx  += G*m1*DX/D3;
                  (moonlets+2)->vy  += G*m1*DY/D3;
                  (moonlets+2)->vz  += G*m1*DZ/D3;
            }
      
            /******** We now consider gravitational interactions between pairs in the same neighbourhood, or between pairs containing a big moonlet. ********/
            int j;
            
            for (j=0; j<how_many_pairs; j++){ //We go over all such pairs. The expected value of how_many_pairs is 0.5*N*how_many_neighbours = O(N)
                  
                  k = (pairs+j)->fst; //The array "pairs" was updated by the function mesh in collision.c
                  p = (pairs+j)->snd;
                  
                  if (*(exists+k) && *(exists+p) && k > 2 && p > 2){
                  
                        /******** Getting the positions, masses and radii ********/
                        X  = (moonlets+k)-> x;
                        Y  = (moonlets+k)-> y;
                        Z  = (moonlets+k)-> z;
                        mk = (moonlets+k)-> mass;
                        Rk = (moonlets+k)-> radius;
                        Xp = (moonlets+p)-> x;
                        Yp = (moonlets+p)-> y;
                        Zp = (moonlets+p)-> z;
                        mp = (moonlets+p)-> mass;
                        Rp = (moonlets+p)-> radius;
                  
                        DX = X-Xp;
                        DY = Y-Yp;
                        DZ = Z-Zp;
                        softening = softening_parameter*(Rk + Rp);
                        D  = sqrt(DX*DX+DY*DY+DZ*DZ+softening*softening);
                        D3 = D*D*D;
                        rk = sqrt(X*X+Y*Y+Z*Z);
                        r3 = rk*rk*rk;
                        rp = sqrt(Xp*Xp+Yp*Yp+Zp*Zp);
                        rp3= rp*rp*rp;
                        /******** Updading the acceleration of moonlet k ********/
                        (moonlets+k)->vx -= G*mp*DX/D3;
                        (moonlets+k)->vy -= G*mp*DY/D3;
                        (moonlets+k)->vz -= G*mp*DZ/D3;
                        /******** Updading the acceleration of moonlet p ********/
                        (moonlets+p)->vx += G*mk*DX/D3;
                        (moonlets+p)->vy += G*mk*DY/D3;
                        (moonlets+p)->vz += G*mk*DZ/D3;           
                  }
            }       
      }
}


void collision(struct moonlet * moonlets, int a, int b, typ f){

      /******** Treats the collision between moonlets a and b. The collision's elasticity is determined     ********/
      /******** by the parameter 1<f<2, as explained is the PDF draft of the project                        ********/
      /******** The array approach given in argument and returned by the function closest_approach contains ********/
      /******** the positions of the moonlets at the collision. See top of function closest_approach        ********/

      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;     //Cartesian speeds of the moonlets.
      typ xa, ya, za, xb, yb, zb;                 //Cartesian positions at the collision
      typ R_a=(moonlets+a)->radius;               //The moonlets' radii
      typ R_b=(moonlets+b)->radius;
      typ R=R_a+R_b;                              //Sum of the radii;
      typ m_a=(moonlets+a)->mass;                 //The moonlets's masses
      typ m_b=(moonlets+b)->mass;              
      
      
      /******** Getting the speeds ********/
      vx_a=(moonlets+a)->vx;
      vy_a=(moonlets+a)->vy;
      vz_a=(moonlets+a)->vz;
      vx_b=(moonlets+b)->vx;
      vy_b=(moonlets+b)->vy;
      vz_b=(moonlets+b)->vz;
      
      
      /******** Getting the positions ********/
      xa=*(approach+0);
      ya=*(approach+1);
      za=*(approach+2);
      xb=*(approach+3);
      yb=*(approach+4);
      zb=*(approach+5);
      
      
      /******** Calculating the new speeds after impact ********/
      typ dr_dot_dv = (xa-xb)*(vx_a-vx_b)+(ya-yb)*(vy_a-vy_b)+(za-zb)*(vz_a-vz_b); //Scalar product dv.dr where dr=r_b-r_a and dv=v_b-v_a are the relative position and speed
      typ alpha = f*m_a*m_b/((m_a+m_b)*R*R);                                       //Factor such that J=alpha*(dv.dr)*dr
      typ Jx = alpha*dr_dot_dv*(xb-xa);                                            //x component of J
      typ Jy = alpha*dr_dot_dv*(yb-ya);                                            //y component of J
      typ Jz = alpha*dr_dot_dv*(zb-za);                                            //z component of J
      
      
      /******** Calculating the speeds after impact ********/
      vx_a+=Jx/m_a;  //v_a=v_a+J/m_a
      vy_a+=Jy/m_a;
      vz_a+=Jz/m_a;
      vx_b-=Jx/m_b;  //v_b=v_b-J/m_b
      vy_b-=Jy/m_b;
      vz_b-=Jz/m_b;
      
      
      /******** Actualizing the speeds and positions in the array moonlets                   ********/
      /******** Moonlets are put back where they were at the beginning of the timestep, but  ********/
      /******** with their post-impact velocity, then they drift with all the other moonlets ********/
      (moonlets+a)->x  = xa - time_until_collision * vx_a;  //Actualizing moonlet a
      (moonlets+a)->y  = ya - time_until_collision * vy_a;
      (moonlets+a)->z  = za - time_until_collision * vz_a;
      (moonlets+a)->vx = vx_a;
      (moonlets+a)->vy = vy_a;
      (moonlets+a)->vz = vz_a;
      (moonlets+b)->x  = xb - time_until_collision * vx_b;  //Actualizing moonlet b
      (moonlets+b)->y  = yb - time_until_collision * vy_b;
      (moonlets+b)->z  = zb - time_until_collision * vz_b;
      (moonlets+b)->vx = vx_b;
      (moonlets+b)->vy = vy_b;
      (moonlets+b)->vz = vz_b;
      

}


void merger(struct moonlet * moonlets, int a, int b){

      /******** Merges moonlet a and moonlet b together ********/


      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;     //Cartesian speeds of the moonlets.
      typ xa, ya, za, xb, yb, zb;                 //Cartesian positions at the collision
      typ m_a = (moonlets+a)->mass;               //The moonlets's masses
      typ m_b = (moonlets+b)->mass;
      typ m = m_a+m_b;                            //Sum of the masses
      typ r_tilde[3]; //Position of the merger
      typ v_tilde[3]; //Speed of the merger
      
      
      /******** Getting the speeds ********/
      vx_a = (moonlets+a)->vx;
      vy_a = (moonlets+a)->vy;
      vz_a = (moonlets+a)->vz;
      vx_b = (moonlets+b)->vx;
      vy_b = (moonlets+b)->vy;
      vz_b = (moonlets+b)->vz;
      
      
      /******** Getting the positions ********/
      xa = *(approach+0);
      ya = *(approach+1);
      za = *(approach+2);
      xb = *(approach+3);
      yb = *(approach+4);
      zb = *(approach+5);
      
      if (!tam_bool){ //The total momentum is conserved
            /******** The speed after merging ********/ 
            v_tilde[0] = (m_a*vx_a+m_b*vx_b)/m;
            v_tilde[1] = (m_a*vy_a+m_b*vy_b)/m;
            v_tilde[2] = (m_a*vz_a+m_b*vz_b)/m; 
                  
            /******** The position after merging ********/
            r_tilde[0] = (m_a*xa+m_b*xb)/m;
            r_tilde[1] = (m_a*ya+m_b*yb)/m;
            r_tilde[2] = (m_a*za+m_b*zb)/m;
      }
      
      else{ //The total angular momentum is conserved
            typ dx=xa-xb, dy=ya-yb, dz=za-zb, dvx=vx_a-vx_b, dvy=vy_a-vy_b, dvz=vz_a-vz_b;
            
            /******** Getting the position of the merger ********/
            typ ra_x_va[3]; //r_a x v_a
            typ rb_x_vb[3]; //r_b x v_b
            cross_product(xa,ya,za,vx_a,vy_a,vz_a,ra_x_va);
            cross_product(xb,yb,zb,vx_b,vy_b,vz_b,rb_x_vb);
            
            typ amG[3]; //Angular momentum of the pair
            amG[0] = m_a * ra_x_va[0] + m_b * rb_x_vb[0];
            amG[1] = m_a * ra_x_va[1] + m_b * rb_x_vb[1];
            amG[2] = m_a * ra_x_va[2] + m_b * rb_x_vb[2];
            
            typ dv_x_dr[3]; //dv x dr
            cross_product(dvx,dvy,dvz,dx,dy,dz,dv_x_dr);
            
            typ G2 = amG[0]*amG[0] + amG[1]*amG[1] + amG[2]*amG[2]; //G**2
            typ delta_r[3];
            typ rb_dv_x_dr = xb*dv_x_dr[0] + yb*dv_x_dr[1] + zb*dv_x_dr[2]; // r_b.(dv x dr)
            typ in_front_of = m_a*m_b/(m*G2)*rb_dv_x_dr;
            delta_r[0] = in_front_of*amG[0];
            delta_r[1] = in_front_of*amG[1];
            delta_r[2] = in_front_of*amG[2];
            
            r_tilde[0] = (m_a*xa + m_b*xb)/m + delta_r[0];
            r_tilde[1] = (m_a*ya + m_b*yb)/m + delta_r[1];
            r_tilde[2] = (m_a*za + m_b*zb)/m + delta_r[2];
            
            /******** Verification step. To be removed later ********/
            /*typ r_tilde_amG = amG[0]*r_tilde[0] + amG[1]*r_tilde[1] + amG[2]*r_tilde[2];
            if (absolute(r_tilde_amG) > 1.0e-13){
                  fprintf(stderr, "Error : r_tilde is not orthogonal to G. |r_tilde . G| = %.15lf\n",absolute(r_tilde_amG));
                  abort();
            }*/
            
            /******** Geting the speed of the merger ********/
            typ r_tilde_vcm = (r_tilde[0]*(m_a*vx_a+m_b*vx_b)+r_tilde[1]*(m_a*vy_a+m_b*vy_b)+r_tilde[2]*(m_a*vz_a+m_b*vz_b))/m; // r_tilde.v_cm
            typ r_tilde2 = r_tilde[0]*r_tilde[0]+r_tilde[1]*r_tilde[1]+r_tilde[2]*r_tilde[2]; //r_tilde**2
            in_front_of = r_tilde_vcm/r_tilde2;
            typ G_x_r_tilde[3]; // G x r_tilde;
            cross_product(amG[0],amG[1],amG[2],r_tilde[0],r_tilde[1],r_tilde[2],G_x_r_tilde);
            typ in_front_of_2 = 1.0/(m*r_tilde2);
            
            v_tilde[0] = in_front_of*r_tilde[0] + in_front_of_2*G_x_r_tilde[0];
            v_tilde[1] = in_front_of*r_tilde[1] + in_front_of_2*G_x_r_tilde[1];
            v_tilde[2] = in_front_of*r_tilde[2] + in_front_of_2*G_x_r_tilde[2];
            
            /******** Verification step. To be removed later ********/
            /*typ G_after[3];
            cross_product(m*r_tilde[0],m*r_tilde[1],m*r_tilde[2],v_tilde[0],v_tilde[1],v_tilde[2],G_after);
            typ variation = sqrt((G_after[0]-amG[0])*(G_after[0]-amG[0])+(G_after[1]-amG[1])*(G_after[1]-amG[1])+(G_after[2]-amG[2])*(G_after[2]-amG[2]));
            if (variation > 1.0e-13){
                  fprintf(stderr, "Error : The total angular momentum changed upon merging. |G_1-G_0| = %.15lf\n",variation);
                  abort();
            }*/
      }
      
      /******** Actualizing the position ********/
      (moonlets+a) -> x  = r_tilde[0] - time_until_collision * v_tilde[0];
      (moonlets+a) -> y  = r_tilde[1] - time_until_collision * v_tilde[1];
      (moonlets+a) -> z  = r_tilde[2] - time_until_collision * v_tilde[2];
      
      /******** Actualizing the position ********/
      (moonlets+a) -> vx = v_tilde[0];
      (moonlets+a) -> vy = v_tilde[1];
      (moonlets+a) -> vz = v_tilde[2];
      
      /******** Actualizing the mass and radius ********/
      (moonlets+a) -> mass = m;
      (moonlets+a) -> radius = pow(3.0*m/(4.0*M_PI*density),1.0/3.0);
      
      /******** Moonlet b does not exist anymore and we disallow moonlet a to collide again for that timestep ********/        
      lose_moonlet(b);
      *(did_collide+a) = 1;

}

void fragmentation(struct moonlet * moonlets, int a, int b){

      /******** Treats the fragmentation due to the collision between moonlets a and b ********/
      
      typ C1_3mu = pow(C1_parameter, 3.0*mu_parameter);
      typ stigma = (3.0*mu_parameter-1.0)/(3.0*mu_parameter);
      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;                  //Cartesian speeds of the moonlets.
      typ xa, ya, za, xb, yb, zb;                              //Cartesian positions at the collision
      typ dx, dy, dz, dvx, dvy, dvz;                           //dx is xa-xb, and so on.
      typ R_a = (moonlets+a)->radius;                          //The moonlets' radii
      typ R_b = (moonlets+b)->radius;
      typ R = R_a+R_b;                                         //Sum of the radii;
      typ m_a = (moonlets+a) -> mass;                          //The moonlets's masses
      typ m_b = (moonlets+b) -> mass;
      typ M = m_a+m_b;                                         //Sum of the masses
      typ m_1;                                                 //Mass of the impactor
      
      /******** Defining the impactor and the target. ********/
      if (R_a > R_b){ // a is the target and b the impactor
            m_1 = m_b;
      }
      else { // a is the impactor and b the target
            m_1 = m_a;
      }            
      
      /******** Getting the speeds at the impact ********/
      vx_a = (moonlets+a) -> vx;
      vy_a = (moonlets+a) -> vy;
      vz_a = (moonlets+a) -> vz;
      vx_b = (moonlets+b) -> vx;
      vy_b = (moonlets+b) -> vy;
      vz_b = (moonlets+b) -> vz;
      
      /******** Getting the positions at the impact ********/
      xa = * approach   ;
      ya = *(approach+1);
      za = *(approach+2);
      xb = *(approach+3);
      yb = *(approach+4);
      zb = *(approach+5);
      
      /******** Getting the relative positions and velocities ********/
      dx  = xa-xb;
      dy  = ya-yb;
      dz  = za-zb;
      dvx = vx_a-vx_b;
      dvy = vy_a-vy_b;
      dvz = vz_a-vz_b;
      
      typ dr_dot_dv = dx*dvx+dy*dvy+dz*dvz; //  (r_1-r_2).(v_1-v_2)
      typ dv_norm = sqrt(dvx*dvx+dvy*dvy+dvz*dvz); //  ||v_1-v_2||
      typ costheta = -dr_dot_dv/(R*dv_norm); //Cosine of impact angle;
      typ R_eq = pow(3.0*M/(4.0*M_PI*density),1.0/3.0); //Hypothetical radius if a and b were to merge
      typ vesc = sqrt(2.0*G*M/R_eq); //Escape velocity at the surface if a and b were to merge
      typ m_check = 3.0*k_parameter*C1_3mu*m_1/(4.0*M_PI)*pow(dv_norm*costheta/vesc,3.0*mu_parameter); //Mass of the tail (not gravitationally bounded to the largest fragment)
      typ m_tilde = M-m_check; //Mass of the largest fragment
      typ m_tilde_2 = m_check / (typ) N_tilde;
      
      /******** Merger case. Moonlets a and b merge together ********/
      
      if (m_check < frag_threshold && m_tilde >= 0.1 * M){
            merger(moonlets, a ,b);
            merger_count++;
            catastrophicity(m_tilde,M);
            return;
      }
      
      
      /******** Some data used for the fragmentation case ********/
      
      typ ra_x_va[3]; //r_a x v_a
      typ rb_x_vb[3]; //r_b x v_b
      typ amG[3]; //Angular momentum to be preserved
      typ v_cm[3], r_cm[3]; //Velocity and position of the center of mass before impact
      typ r_tilde[3]; //Position of the largest fragment
      typ v_tilde[3]; //Speed of the largest fragment
      
      cross_product(xa,ya,za,vx_a,vy_a,vz_a,ra_x_va);
      cross_product(xb,yb,zb,vx_b,vy_b,vz_b,rb_x_vb);
      amG[0] = m_a * ra_x_va[0] + m_b * rb_x_vb[0];
      amG[1] = m_a * ra_x_va[1] + m_b * rb_x_vb[1];
      amG[2] = m_a * ra_x_va[2] + m_b * rb_x_vb[2];
      v_cm[0]=(m_a*vx_a + m_b*vx_b)/M;  v_cm[1]=(m_a*vy_a + m_b*vy_b)/M;  v_cm[2]=(m_a*vz_a + m_b*vz_b)/M;
      r_cm[0]=(m_a*xa   + m_b*  xb)/M;  r_cm[1]=(m_a*ya   + m_b*  yb)/M;  r_cm[2]=(m_a*za   + m_b*  zb)/M;
      
      
      /******** Super-catastrophic fragmentation. The mass of the largest fragment is less than 10 % of the total mass. The ejecta is discarded ********/
      
      if (m_tilde < 0.1*M){ //m_tilde is not proportionnal to Qr/Qr* in this regime
            
            m_tilde = 0.1* M * pow(3.0*k_parameter*C1_3mu/(2.0*1.8*M_PI)*m_1/M,-1.5) * pow(costheta*dv_norm/vesc,-9.0*mu_parameter/2.0); //Eq. (44) of Leinhardt and Stewart (2012)
            
            /******** Only the largest fragment remains in the super-catastrophic regime ********/
            /******** Actualizing its speed ********/
            (moonlets+a) -> vx = (m_a*vx_a+m_b*vx_b)/M;
            (moonlets+a) -> vy = (m_a*vy_a+m_b*vy_b)/M;
            (moonlets+a) -> vz = (m_a*vz_a+m_b*vz_b)/M;
            /******** Actualizing its position ********/
            (moonlets+a) -> x = (m_a*xa+m_b*xb)/M - time_until_collision * (moonlets+a)->vx;
            (moonlets+a) -> y = (m_a*ya+m_b*yb)/M - time_until_collision * (moonlets+a)->vy;
            (moonlets+a) -> z = (m_a*za+m_b*zb)/M - time_until_collision * (moonlets+a)->vz;
            /******** Actualizing its mass and radius ********/
            (moonlets+a) -> mass   = m_tilde;
            (moonlets+a) -> radius = pow(3.0*m_tilde/(4.0*M_PI*density),1.0/3.0);
            /******** Moonlet b does not exist anymore and we disallow moonlet a to collide again for that timestep ********/
            lose_moonlet(b);
            *(did_collide+a) = 1;
            super_catastrophic_count++;
            catastrophicity(m_tilde,M);
            
            /******** To be removed later ********/
            /*typ angular_momentum_loss[3];
            typ X,Y,Z,vX,vY,vZ,m;
            X=(moonlets+a)-> x;
            Y=(moonlets+a)-> y;
            Z=(moonlets+a)-> z;
            vX=(moonlets+a)-> vx;
            vY=(moonlets+a)-> vy;
            vZ=(moonlets+a)-> vz;
            m=(moonlets+a)-> mass;
            cross_product(m*X,m*Y,m*Z,vX,vY,vZ,angular_momentum_loss);
            *tam_loss += amG[0]-angular_momentum_loss[0];
            *(tam_loss+1) += amG[1]-angular_momentum_loss[1];
            *(tam_loss+2) += amG[2]-angular_momentum_loss[2];*/
            
            return;
      }
      
      
      /******** Partial fragmentation. The tail is reunited into a single moonlet. ********/
      
      if (m_check >= frag_threshold && m_tilde_2 < frag_threshold){
            typ r_k[3]; //Position of the tail with respect to the largest fragment
            typ v_k[3]; //Velocity of the tail with respect to the largest fragment
            typ v_k_scalar = vesc/stigma; //Scalar velocity of the tail with respect to the largest fragment
            typ in_front_of = v_k_scalar/R;
            if (R_a > R_b){
                  r_k[0] = -dx;  r_k[1] = -dy;  r_k[2] = -dz;
            }
            else {
                  r_k[0] = dx;  r_k[1] = dy;  r_k[2] = dz;
            }
            v_k[0] = in_front_of * r_k[0];  v_k[1] = in_front_of * r_k[1];  v_k[2] = in_front_of * r_k[2];
            typ s_tilde[3]; //See Sect. about conservation of the total angular momentum of the PDF draft for definition of s_tilde, u_tilde and g_tilde
            typ u_tilde[3];
            typ g_tilde[3];
            s_tilde[0] = m_check*r_k[0];  s_tilde[1] = m_check*r_k[1];  s_tilde[2] = m_check*r_k[2];
            u_tilde[0] = m_check*v_k[0];  u_tilde[1] = m_check*v_k[1];  u_tilde[2] = m_check*v_k[2];
            cross_product(s_tilde[0],s_tilde[1],s_tilde[2],v_k[0],v_k[1],v_k[2],g_tilde);
            
            /******** Determination of r_tilde ********/
            typ mathfrak_a[3];
            typ s_tilde_x_u_tilde[3];
            cross_product(s_tilde[0],s_tilde[1],s_tilde[2],u_tilde[0],u_tilde[1],u_tilde[2],s_tilde_x_u_tilde);
            mathfrak_a[0] = M*(amG[0]-g_tilde[0]) + s_tilde_x_u_tilde[0];
            mathfrak_a[1] = M*(amG[1]-g_tilde[1]) + s_tilde_x_u_tilde[1];
            mathfrak_a[2] = M*(amG[2]-g_tilde[2]) + s_tilde_x_u_tilde[2];
            typ r_cm_dot_mathfrak_a = r_cm[0]*mathfrak_a[0] + r_cm[1]*mathfrak_a[1] + r_cm[2]*mathfrak_a[2];
            in_front_of = -r_cm_dot_mathfrak_a/(mathfrak_a[0]*mathfrak_a[0] + mathfrak_a[1]*mathfrak_a[1] + mathfrak_a[2]*mathfrak_a[2]);
            
            typ delta_r_tilde[3]; // \delta \tilde{r}
            delta_r_tilde[0] = in_front_of * mathfrak_a[0];
            delta_r_tilde[1] = in_front_of * mathfrak_a[1];
            delta_r_tilde[2] = in_front_of * mathfrak_a[2];
            
            typ vector_a[3];
            vector_a[0] = M*(r_cm[0]+delta_r_tilde[0]);  vector_a[1] = M*(r_cm[1]+delta_r_tilde[1]);  vector_a[2] = M*(r_cm[2]+delta_r_tilde[2]);
            r_tilde[0] = (vector_a[0]-s_tilde[0])/M;
            r_tilde[1] = (vector_a[1]-s_tilde[1])/M;
            r_tilde[2] = (vector_a[2]-s_tilde[2])/M;
            
            /******** Determination of v_tilde ********/
            if (tam_bool){ //Total angular momentum is conserved
                  typ vector_b[3];
                  typ r_tilde_x_u_tilde[3];
                  cross_product(r_tilde[0],r_tilde[1],r_tilde[2],u_tilde[0],u_tilde[1],u_tilde[2],r_tilde_x_u_tilde);
                  vector_b[0] = amG[0] - r_tilde_x_u_tilde[0] - g_tilde[0];
                  vector_b[1] = amG[1] - r_tilde_x_u_tilde[1] - g_tilde[1];
                  vector_b[2] = amG[2] - r_tilde_x_u_tilde[2] - g_tilde[2];
                  typ b_x_a[3];
                  cross_product(vector_b[0],vector_b[1],vector_b[2],vector_a[0],vector_a[1],vector_a[2],b_x_a);
                  typ alpha = (v_cm[0]-u_tilde[0]/M)*vector_a[0] + (v_cm[1]-u_tilde[1]/M)*vector_a[1] + (v_cm[2]-u_tilde[2]/M)*vector_a[2];
                  in_front_of = 1.0/(vector_a[0]*vector_a[0] + vector_a[1]*vector_a[1] + vector_a[2]*vector_a[2]);
                  alpha *= in_front_of;
                  v_tilde[0] = in_front_of * b_x_a[0] + alpha * vector_a[0];
                  v_tilde[1] = in_front_of * b_x_a[1] + alpha * vector_a[1];
                  v_tilde[2] = in_front_of * b_x_a[2] + alpha * vector_a[2];
            }
            else{ //Total momentum is conserved
                  v_tilde[0] = v_cm[0] - m_check/M*v_k[0];
                  v_tilde[1] = v_cm[1] - m_check/M*v_k[1];
                  v_tilde[2] = v_cm[2] - m_check/M*v_k[2];
            }
            
            /******** Actualizing the moonlets ********/
            typ R_tilde = pow(3.0*m_tilde/(4.0*M_PI*density),1.0/3.0); //Radius of the largest fragment
            typ R_check = pow(3.0*m_check/(4.0*M_PI*density),1.0/3.0); //Radius of the tail fragment
            (moonlets+a) -> x = r_tilde[0] - time_until_collision * v_tilde[0];
            (moonlets+a) -> y = r_tilde[1] - time_until_collision * v_tilde[1];
            (moonlets+a) -> z = r_tilde[2] - time_until_collision * v_tilde[2];
            (moonlets+a) -> vx = v_tilde[0];
            (moonlets+a) -> vy = v_tilde[1];
            (moonlets+a) -> vz = v_tilde[2];
            (moonlets+a) -> mass = m_tilde;
            (moonlets+a) -> radius = R_tilde;
            (moonlets+b) -> x  = r_tilde[0]+r_k[0] - time_until_collision * (v_tilde[0]+v_k[0]);
            (moonlets+b) -> y  = r_tilde[1]+r_k[1] - time_until_collision * (v_tilde[1]+v_k[1]);
            (moonlets+b) -> z  = r_tilde[2]+r_k[2] - time_until_collision * (v_tilde[2]+v_k[2]);
            (moonlets+b) -> vx = v_tilde[0]+v_k[0];
            (moonlets+b) -> vy = v_tilde[1]+v_k[1];
            (moonlets+b) -> vz = v_tilde[2]+v_k[2];
            (moonlets+b) -> mass = m_check;
            (moonlets+b) -> radius = R_check;
            
            /******** Verification step. To be removed later ********/
            /*typ G_after_tilde[3];
            typ G_after_check[3];
            typ G_after[3];
            cross_product(m_tilde*r_tilde[0],m_tilde*r_tilde[1],m_tilde*r_tilde[2],v_tilde[0],v_tilde[1],v_tilde[2],G_after_tilde);
            cross_product(m_check*(r_tilde[0]+r_k[0]),m_check*(r_tilde[1]+r_k[1]),m_check*(r_tilde[2]+r_k[2]),v_tilde[0]+v_k[0],v_tilde[1]+v_k[1],v_tilde[2]+v_k[2],G_after_check);
            G_after[0] = G_after_tilde[0]+G_after_check[0];
            G_after[1] = G_after_tilde[1]+G_after_check[1];
            G_after[2] = G_after_tilde[2]+G_after_check[2];
            typ variation = sqrt((G_after[0]-amG[0])*(G_after[0]-amG[0])+(G_after[1]-amG[1])*(G_after[1]-amG[1])+(G_after[2]-amG[2])*(G_after[2]-amG[2]));
            if (variation > 1.0e-13){
                  fprintf(stderr, "Error : The total angular momentum changed upon fragmentation. |G_1-G_0| = %.15lf\n",variation);
                  abort();
            }*/
            
            if (mutual_bool && mesh_bool && !force_naive_bool){ //If the mutual gravitational interactions are considered, we register the pair (a,b) to be taken care of
                  (pairs+how_many_pairs)->fst = a;
                  (pairs+how_many_pairs)->snd = b;
                  how_many_pairs++;
            }
            
            /******** We disallow moonlets a and be to collide again during that timestep ********/
            *(did_collide+a) = 1;
            *(did_collide+b) = 1;
            half_fragmentation_count++;
            catastrophicity(m_tilde,M);
            return;
      }
      
      
      /******** Full fragmentation. The tail is made up of N_tilde moonlets ********/

      if (m_tilde_2 >= frag_threshold){
      
            /******** Determination of the r_k' and v_k' ********/ 
            typ r_k[3*N_tilde]; //Position of the fragments of the tail with respect to the largest fragment
            typ v_k[3*N_tilde]; //Velocity of the fragments of the tail with respect to the largest fragment
            typ dr[3]; //r_1-r_2
            typ dv[3]; //v_1-v_2
            typ R_tilde_2 = pow(3.0*m_tilde_2/(4.0*M_PI*density),1.0/3.0); //Radius of the fragments of the tail
            typ u[3]; //Vectors used to locate the fragments of the tail
            typ v[3];
            if (R_a > R_b){
                  dr[0] = -dx;   dr[1] = -dy;   dr[2] = -dz;
                  dv[0] = -dvx;  dv[1] = -dvy;  dv[2] = -dvz;
            }
            else {
                  dr[0] = dx;   dr[1] = dy;   dr[2] = dz;
                  dv[0] = dvx;  dv[1] = dvy;  dv[2] = dvz;
            }
            typ dr_x_dv[3]; // dr x dv
            typ dr_x_dv_norm; // ||dr x dv||
            int pq[4] = pq_min_max; // {p_k_min,p_k_max,q_k_min,q_k_max}
            int p,q;
            cross_product(dr[0],dr[1],dr[2],dv[0],dv[1],dv[2],dr_x_dv);
            dr_x_dv_norm = sqrt(dr_x_dv[0]*dr_x_dv[0] + dr_x_dv[1]*dr_x_dv[1] + dr_x_dv[2]*dr_x_dv[2]);
            if (dr_x_dv_norm > 1.0e-5){ //Oblique collision
                  /******** v = (dr x dv) / ||dr x dv|| ********/
                  v[0] = dr_x_dv[0]/dr_x_dv_norm;  v[1] = dr_x_dv[1]/dr_x_dv_norm;  v[2] = dr_x_dv[2]/dr_x_dv_norm; //Defining the unit vector v
            }
            else { //Nearly Frontal collision
                  /******** v is any unit vector orthogonal to dr ********/
                  int m = maximum(absolute(dr[0]), absolute(dr[1]), absolute(dr[2]));
                  typ numerator;
                  if (m == 0){
                        numerator = dr[1]+dr[2];
                        v[1] = 1.0/sqrt(2.0+numerator*numerator/(dr[0]*dr[0])); //Defining the unit vector v
                        v[2] = v[1];
                        v[0] = -v[1]*numerator/dr[0];
                  }
                  else if (m == 1){
                        numerator = dr[0]+dr[2];
                        v[0] = 1.0/sqrt(2.0+numerator*numerator/(dr[1]*dr[1])); //Defining the unit vector v
                        v[2] = v[0];
                        v[1] = -v[0]*numerator/dr[1];
                  }
                  else {
                        numerator = dr[0]+dr[1];
                        v[0] = 1.0/sqrt(2.0+numerator*numerator/(dr[2]*dr[2])); //Defining the unit vector v
                        v[1] = v[0];
                        v[2] = -v[0]*numerator/dr[2];
                  }
            }
            typ v_x_dr[3]; // v x dr
            cross_product(v[0],v[1],v[2],dr[0],dr[1],dr[2],v_x_dr);
            u[0] = v_x_dr[0]/R;  u[1] = v_x_dr[1]/R;  u[2] = v_x_dr[2]/R; //Defining the unit vector u
            int n = 0;
            typ v_k_scalar, two_p_R_tilde_2, two_q_R_tilde_2, in_front_of;
            for (p = pq[0]; p <= pq[1]; p++){ //We travel along the rectangle of integer coordinate points to define the position and speeds of the fragments of the tail
                  for (q = pq[2]; q <= pq[3]; q++){
                        two_p_R_tilde_2 = ((typ) (2*p))*R_tilde_2;  two_q_R_tilde_2 = ((typ) (2*q))*R_tilde_2;
                        r_k[3*n] = dr[0] + two_p_R_tilde_2*u[0] + two_q_R_tilde_2*v[0];   // x-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        r_k[3*n+1] = dr[1] + two_p_R_tilde_2*u[1] + two_q_R_tilde_2*v[1]; // y-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        r_k[3*n+2] = dr[2] + two_p_R_tilde_2*u[2] + two_q_R_tilde_2*v[2]; // z-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        v_k_scalar = vesc * N_tilde/stigma *(pow(1.0-((typ) n)/N_tilde,stigma)-pow(1.0-((typ) (n+1))/N_tilde,stigma)); //Scalar velocity of the (n+1)^th fragment of the tail
                        in_front_of = v_k_scalar/sqrt(R*R+4.0*R_tilde_2*R_tilde_2*(p*p+q*q));
                        v_k[3*n] = in_front_of * r_k[3*n];  v_k[3*n+1] = in_front_of * r_k[3*n+1];  v_k[3*n+2] = in_front_of * r_k[3*n+2];
                        n++;
                  }
            }

            
            /******** Determination of u_tilde, s_tilde and g_tilde. See Sect. about conservation of the total angular momentum of the PDF draft for their definition ********/
            
            typ s_tilde[3] = {0.0, 0.0, 0.0};
            typ u_tilde[3] = {0.0, 0.0, 0.0};
            typ g_tilde[3] = {0.0, 0.0, 0.0};
            typ rk_x_vk[3];
            for (n = 0; n < N_tilde; n++){
                  s_tilde[0] += m_tilde_2*r_k[3*n];  s_tilde[1] += m_tilde_2*r_k[3*n+1];  s_tilde[2] += m_tilde_2*r_k[3*n+2];
                  u_tilde[0] += m_tilde_2*v_k[3*n];  u_tilde[1] += m_tilde_2*v_k[3*n+1];  u_tilde[2] += m_tilde_2*v_k[3*n+2];
                  cross_product(r_k[3*n],r_k[3*n+1],r_k[3*n+2],v_k[3*n],v_k[3*n+1],v_k[3*n+2],rk_x_vk);
                  g_tilde[0] += m_tilde_2*rk_x_vk[0];  g_tilde[1] += m_tilde_2*rk_x_vk[1];  g_tilde[2] += m_tilde_2*rk_x_vk[2];
            }

            
            /******** Determination of r_tilde ********/
            
            typ mathfrak_a[3];
            typ s_tilde_x_u_tilde[3];
            cross_product(s_tilde[0],s_tilde[1],s_tilde[2],u_tilde[0],u_tilde[1],u_tilde[2],s_tilde_x_u_tilde);
            mathfrak_a[0] = M*(amG[0]-g_tilde[0]) + s_tilde_x_u_tilde[0];
            mathfrak_a[1] = M*(amG[1]-g_tilde[1]) + s_tilde_x_u_tilde[1];
            mathfrak_a[2] = M*(amG[2]-g_tilde[2]) + s_tilde_x_u_tilde[2];
            typ r_cm_dot_mathfrak_a = r_cm[0]*mathfrak_a[0] + r_cm[1]*mathfrak_a[1] + r_cm[2]*mathfrak_a[2];
            in_front_of = -r_cm_dot_mathfrak_a/(mathfrak_a[0]*mathfrak_a[0] + mathfrak_a[1]*mathfrak_a[1] + mathfrak_a[2]*mathfrak_a[2]);
            
            typ delta_r_tilde[3]; // \delta \tilde{r}
            delta_r_tilde[0] = in_front_of * mathfrak_a[0];
            delta_r_tilde[1] = in_front_of * mathfrak_a[1];
            delta_r_tilde[2] = in_front_of * mathfrak_a[2];
            
            typ vector_a[3];
            vector_a[0] = M*(r_cm[0]+delta_r_tilde[0]);  vector_a[1] = M*(r_cm[1]+delta_r_tilde[1]);  vector_a[2] = M*(r_cm[2]+delta_r_tilde[2]);
            r_tilde[0] = (vector_a[0]-s_tilde[0])/M;
            r_tilde[1] = (vector_a[1]-s_tilde[1])/M;
            r_tilde[2] = (vector_a[2]-s_tilde[2])/M;

            
            /******** Determination of v_tilde ********/
            
            if (tam_bool){ //Total angular momentum is conserved
                  typ vector_b[3];
                  typ r_tilde_x_u_tilde[3];
                  cross_product(r_tilde[0],r_tilde[1],r_tilde[2],u_tilde[0],u_tilde[1],u_tilde[2],r_tilde_x_u_tilde);
                  vector_b[0] = amG[0] - r_tilde_x_u_tilde[0] - g_tilde[0];
                  vector_b[1] = amG[1] - r_tilde_x_u_tilde[1] - g_tilde[1];
                  vector_b[2] = amG[2] - r_tilde_x_u_tilde[2] - g_tilde[2];
                  typ b_x_a[3];
                  cross_product(vector_b[0],vector_b[1],vector_b[2],vector_a[0],vector_a[1],vector_a[2],b_x_a);
                  typ alpha = (v_cm[0]-u_tilde[0]/M)*vector_a[0] + (v_cm[1]-u_tilde[1]/M)*vector_a[1] + (v_cm[2]-u_tilde[2]/M)*vector_a[2];
                  in_front_of = 1.0/(vector_a[0]*vector_a[0] + vector_a[1]*vector_a[1] + vector_a[2]*vector_a[2]);
                  alpha *= in_front_of;
                  v_tilde[0] = in_front_of * b_x_a[0] + alpha * vector_a[0];
                  v_tilde[1] = in_front_of * b_x_a[1] + alpha * vector_a[1];
                  v_tilde[2] = in_front_of * b_x_a[2] + alpha * vector_a[2];
            }
            else{ //Total momentum is conserved
                  typ corr[3] = {0.0, 0.0, 0.0};
                  for (n = 0; n < N_tilde; n++){
                        corr[0] += v_k[3*n]  ;
                        corr[1] += v_k[3*n+1];
                        corr[2] += v_k[3*n+2];
                  }
                  v_tilde[0] = v_cm[0] - m_tilde_2/M * corr[0];
                  v_tilde[1] = v_cm[1] - m_tilde_2/M * corr[1];
                  v_tilde[2] = v_cm[2] - m_tilde_2/M * corr[2];
            }

            
            /******** Managing the indexes of all the moonlets ********/
            
            int id[N_tilde+1]; //The N_tilde+1 indexes where the N_tilde+1 moonlets will be stored in the array moonlets
            id[0] = a; //Putting the largest fragment there
            for (n = 1; n <= N_tilde; n++){
                  if (mesh_bool && !force_naive_bool){ //By construction of the mesh algorithm, the fragments must be put at the end 
                        id[n] = get_free_index(1); //Putting the remaining N_tilde fragments of the tail there
                  }
                  else{ //No need to put the fragments at the end
                        id[n] = get_free_index(0); //Putting the remaining N_tilde fragments of the tail there
                  }
            }
            if (mutual_bool && mesh_bool && !force_naive_bool){ //Adding the N_tilde(N_tilde+1)/2 pairs to be taken into account for gravitational interactions
                  for (p = 0; p <= N_tilde; p++){
                        for (q = 0; q < p; q++){
                              (pairs+how_many_pairs) -> fst = id[p];
                              (pairs+how_many_pairs) -> snd = id[q];
                              how_many_pairs ++;
                        }
                  }
            }
            
            lose_moonlet(b); //Moonlet b does not exist anymore
            for (n = 0; n <= N_tilde; n++){
                  *(did_collide+id[n]) = 1; //Moonlets won't be able to collide during that timestep
            }

            /******** Actualizing the largest fragment ********/
            (moonlets+id[0]) -> x      = r_tilde[0] - time_until_collision * v_tilde[0];
            (moonlets+id[0]) -> y      = r_tilde[1] - time_until_collision * v_tilde[1];
            (moonlets+id[0]) -> z      = r_tilde[2] - time_until_collision * v_tilde[2];
            (moonlets+id[0]) -> vx     = v_tilde[0];
            (moonlets+id[0]) -> vy     = v_tilde[1];
            (moonlets+id[0]) -> vz     = v_tilde[2];
            (moonlets+id[0]) -> mass   = m_tilde;
            (moonlets+id[0]) -> radius = pow(3.0*m_tilde/(4.0*M_PI*density),1.0/3.0);

            /******** Actualizing the fragments of the tail ********/
            for (n = 0; n < N_tilde; n++){
                  (moonlets+id[n+1]) -> x      = r_tilde[0]+r_k[3*n]   - time_until_collision * (v_tilde[0]+v_k[3*n])  ;
                  (moonlets+id[n+1]) -> y      = r_tilde[1]+r_k[3*n+1] - time_until_collision * (v_tilde[1]+v_k[3*n+1]);
                  (moonlets+id[n+1]) -> z      = r_tilde[2]+r_k[3*n+2] - time_until_collision * (v_tilde[2]+v_k[3*n+2]);
                  (moonlets+id[n+1]) -> vx     = v_tilde[0]+v_k[3*n];
                  (moonlets+id[n+1]) -> vy     = v_tilde[1]+v_k[3*n+1];
                  (moonlets+id[n+1]) -> vz     = v_tilde[2]+v_k[3*n+2];
                  (moonlets+id[n+1]) -> mass   = m_tilde_2;
                  (moonlets+id[n+1]) -> radius = R_tilde_2;
            }
            full_fragmentation_count++;
            catastrophicity(m_tilde, M);
            return;
      }
}


void collision_treatment(struct moonlet * moonlets, int a, int b, int type_of_collision){

      /******** Treats the collision between moonlet a and moonlet b            ********/
      /******** The integer type_of_collision determines the type of collision. ********/
      /******** 0 --> Elastic collision                                         ********/
      /******** 1 --> Inelastic collision                                       ********/
      /******** 2 --> Instant merger                                            ********/
      /******** 3 --> Fragmentation                                             ********/
      
      if (!(*(exists+a) && *(exists+b))){
            return;
      }
      
      if (type_of_collision == 0){
            collision(moonlets, a, b, 2.0);
            *(did_collide+a) = 1;
            *(did_collide+b) = 1;
            if (mutual_bool && mesh_bool && !force_naive_bool){ //If the mutual gravity is considered, we register the pair (a,b) to be taken care of
                  (pairs+how_many_pairs)->fst = a;
                  (pairs+how_many_pairs)->snd = b;
                  how_many_pairs++;
            }
      }
      
      else if (type_of_collision == 1){
            collision(moonlets, a, b, collision_parameter);
            *(did_collide+a) = 1;
            *(did_collide+b) = 1;
            if (mutual_bool && mesh_bool && !force_naive_bool){ //If the mutual gravity is considered, we register the pair (a,b) to be taken care of
                  (pairs+how_many_pairs)->fst = a;
                  (pairs+how_many_pairs)->snd = b;
                  how_many_pairs++;
            }
      }
      
      else if (type_of_collision == 2){
            merger(moonlets, a, b);
      }
      
      else if (type_of_collision == 3){
            fragmentation(moonlets, a, b);
      }
      collision_count++;
}


void get_neighbours_mesh(struct moonlet * moonlets){


      /******** When mutual interactions are treated with the mesh algorithm and neighbouring pairs ********/
      /******** were not previously found by the function mesh (e.g. if collision_bool is 0)        ********/


      int k,p;
      typ * the_approach;
      int current_largest_id = largest_id;
      
      for (k = 0; k <= current_largest_id; k++){
            if(*(exists+k)){ //Checking if there is a moonlet in the k^th cell of the array moonlets
                  neighbours(moonlets, k); //Adding the moonlet to the hash table, and retrieving its neighbours in the chain nghb
                        
                  while(nghb -> how_many > 0){ //If the moonlet is inside the collision cube and has at least one neighbour
                        p = (nghb -> ids)[nghb -> how_many - 1]; // p is the id of a moonlet neighbour to k
                        (pairs + how_many_pairs) -> fst = k;
                        (pairs + how_many_pairs) -> snd = p;
                        how_many_pairs ++;
                        nghb = partial_delete(nghb);
                        total_neighbours ++;
                  }
            }
      }
}





