//Definition of member function in LBM.h
#include<iostream>
#include<stdio.h>
#include<cmath>
#include<omp.h>
#include<vector>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <bits/stdc++.h>
#include<string>
#include "eigen-3.4.0/Eigen/Dense"
#include "LBM.h"
using namespace std;
using namespace Eigen;

LBM::LBM(int Nx, int Ny, int Nz, double tau) : Nx(Nx), Ny(Ny), Nz(Nz), tau(tau) {
	double omega = dt/tau;
	//Memory allocation for LATTICE **fluid
	fluid = new LATTICE **[Nx+2];
	
	for (int i = 0; i<Nx+2; i++) {
		fluid[i] = new LATTICE *[Ny+2];
		for(int j = 0; j<Ny+2; j++) {
			fluid[i][j] = new LATTICE [Nz+2];
		}
	}
}


void LBM::Init() { //Initialization using equilibrium 
	
	#pragma omp parallel for 
	for (int i = 0; i<=Nx+1; i++) { //x direction
		for (int j = 0; j<=Ny+1; j++) {  //y direction
			for(int k = 0; k<=Nz+1; k++) { // z direction
				// if(fluid[i][j][k].type == TYPE_F || fluid[i][j][k].type == TYPE_OUT || fluid[i][j][k].type == TYPE_IN || fluid[i][j][k].type == TYPE_L)  {
					
					for (int l = 0; l<npop; l++) {
						double f_eq = Calc_feq(i,j,k,l);
						fluid[i][j][k].f[l] = f_eq;
						fluid[i][j][k].f_star[l] = f_eq;

						#ifdef ENERGY
							vector<double> uwall(3,0);
							double g_eq = Calc_geq(i,j,k,l, fluid[i][j][k].type, uwall);
							fluid[i][j][k].g[l] = g_eq;
							fluid[i][j][k].g_star[l] = g_eq;
						#endif
					}
				// }

			}
		}
	}
	
}

//Terminating in Collision function
void LBM::Collision() {
	#pragma omp parallel for
	for(int i = 0; i<=Nx+1; i++) {
		for (int j = 0; j<=Ny+1; j++) {
			for(int k = 0; k<=Nz+1; k++) {
				if(fluid[i][j][k].type == TYPE_F) {
					#if defined Conventional_LBM
						ConvLBM(i,j,k);
					#endif
				}
			}
		}
	}
}

void LBM::ConvLBM(int i, int j, int k) {
	double omega = dt/tau;
	//Macroscopic property of the fluid 
	double rho = fluid[i][j][k].density;
	double ux = fluid[i][j][k].ux;
	double uy = fluid[i][j][k].uy;
	double uz = fluid[i][j][k].uz;
	double Fx = fluid[i][j][k].F[0];
	double Fy = fluid[i][j][k].F[1];
	double Fz = fluid[i][j][k].F[2];

	#ifdef ENERGY
		double T = fluid[i][j][k].T;
	#endif
	

	//Velocity Equilibrium distribution 
	double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
	double F_dot_u = Fx*ux + Fy*uy + Fz*uz;
	double c_dot_u[npop];
	for(int l = 0; l<npop; l++) {
		//Preparation for collision term
		c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
		double f_eq = Calc_feq(i,j,k,l);	

		//Preparation for source term
		double F_dot_c = Fx*cx[l] + Fy*cy[l] + Fz*cz[l];
		double Source = (1-omega*0.5)*w[l]*(F_dot_c/cs2 + (F_dot_c)*(c_dot_u[l])/(cs2*cs2) - F_dot_u/cs2 );

		//Collision process
		fluid[i][j][k].f_star[l] = fluid[i][j][k].f[l]*(1-omega) + omega*f_eq + Source*dt;

		#ifdef Cylinder_IBM_Moving_Arbie
			double radius = D/2;
			if(i == floor(center_x - radius - 1) or i == ceil(center_x + radius + 1)) {
				//Finding f with condition : 
				if(fluid[i][j][k].f[l] < 0 or fluid[i][j][k].f[l] >= 1.5) {cout<<"f condition happens at "<<i<<endl;}
			}
		#endif
		if(f_eq <0) {
			cout<<"Terminating in Collision function"<<endl;
			terminate();
		}
	}

	//Energy part
	#ifdef ENERGY
		//Viscous dissipation part : 
		double S[3][3];

		for(int row = 0; row<3; row++) {
			for(int col = 0; col<3; col++) {
				S[row][col] = 0;

				for(int l = 0; l<npop; l++) {
					//Preparation for collision term
					double f_eq = Calc_feq(i,j,k,l);

					//Calculating Q
					if(row == 0 and col == 0) S[row][col] += cx[l]*cx[l]*(fluid[i][j][k].f[l] - f_eq);
					else if(row == 0 and col == 1) S[row][col] += cx[l]*cy[l]*(fluid[i][j][k].f[l] - f_eq);
					else if(row == 0 and col == 2) S[row][col] += cx[l]*cz[l]*(fluid[i][j][k].f[l] - f_eq);
					else if(row == 1 and col == 1) S[row][col] += cy[l]*cy[l]*(fluid[i][j][k].f[l] - f_eq);
					else if(row == 1 and col == 2) S[row][col] += cy[l]*cz[l]*(fluid[i][j][k].f[l] - f_eq);
					else if(row == 2 and col == 2) S[row][col] += cz[l]*cz[l]*(fluid[i][j][k].f[l] - f_eq);
					
				}
			}
		}

		S[1][0] = S[0][1]; S[2][0] = S[0][2]; S[2][1] = S[1][2];

		for(int row = 0; row<3; row++) {
			for(int col = 0; col<3; col++) {
				S[row][col] = -S[row][col]/(2*rho*tau*cs2);
			}
		}

		//Calculate Phi
		double Phi = 0;

		for(int row = 0; row<3; row++) {
			for(int col = 0; col<3; col++) {
				Phi += S[row][col]*S[row][col];
			}
		}

		Phi = 2*nu*rho*Phi;

		//Collision part : 
		for(int l = 0; l<npop; l++) {
			double F_before = fluid[i][j][k].F_vis[l];
			double F = w[l]*Phi/cp*(1 + c_dot_u[l]*(tau_g - 0.5)/(cs2*tau_g));
			fluid[i][j][k].F_vis[l] = F;

			//Calculate dF/dt : 
			double dF_dt = (F - F_before)/dt;

			//Calculate equilibrium term : 
			vector<double> uwall(3,0);
			double g_eq = Calc_geq(i,j,k, l, fluid[i][j][k].type, uwall);
			fluid[i][j][k].g_star[l] = fluid[i][j][k].g[l] - (fluid[i][j][k].g[l] - g_eq)/tau_g + dt*F + 0.5*pow(dt,2)*dF_dt;
		}
	#endif
}


void LBM::Streaming() {	//Streaming and Boundary condition
	#pragma omp parallel for 
	//Streaming for fluid part : 
	for(int i = 0; i<=Nx+1; i++) {
		for(int j = 0; j<=Ny+1; j++) {
			for(int k = 0; k<=Nz+1; k++) {
				if(fluid[i][j][k].type == TYPE_F) {
					int i_nb, j_nb, k_nb;

					double ux = fluid[i][j][k].ux,  uy = fluid[i][j][k].uy, uz = fluid[i][j][k].uz;
					double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);

					//Check the surrounding of fluid whether there is solid or not 
					bool solid = 0;
					for(int l = 0; l<npop; l++) {
						i_nb = (i + cx[l] + (Nx+2))%(Nx+2);
						j_nb = (j + cy[l] + (Ny+2))%(Ny+2);
						k_nb = (k + cz[l] + (Nz+2))%(Nz+2);

						if(fluid[i_nb][j_nb][k_nb].type == TYPE_S) {
							solid = 1; l = npop+1;
						}
					}

					//If there is solid, find the direction with c.n <=0 and c.n > 0
					vector<int> cdotn_pos;
					vector<int> cdotn_notpos;

					//Initialize cdotn_notpos : 
					for(int l = 0; l<npop; l++) cdotn_notpos.push_back(l);

					if(solid) {
						//For cdotn_pos
						for(int l = 0; l<npop; l++) {
							i_nb = (i + cx[l] + (Nx+2))%(Nx+2);
							j_nb = (j + cy[l] + (Ny+2))%(Ny+2);
							k_nb = (k + cz[l] + (Nz+2))%(Nz+2);

							if(fluid[i_nb][j_nb][k_nb].type == TYPE_S) {
								cdotn_pos.push_back(opposite[l]);
							}
							
						}	

						//For cdotn_notpos : 
						for(auto element:cdotn_pos) {
							int valtobeDel = element;
							auto it = find(cdotn_notpos.begin(), cdotn_notpos.end(), valtobeDel);

							if(it != cdotn_notpos.end()) {
								cdotn_notpos.erase(it);
							}	
						}
						

						// cout<<"cdotn_pos (i,j,k) = ("<<i<<", "<<j<<", "<<k<<") : ";
						// for(int l = 0; l<cdotn_pos.size(); l++) {
						// 	cout<<cdotn_pos[l]<<" ";
						// }
						// cout<<"\n";

						// cout<<"cdotn_notpos (i,j,k) = ("<<i<<", "<<j<<", "<<k<<") : ";
						// for(int l = 0; l<cdotn_notpos.size(); l++) {
						// 	cout<<cdotn_notpos[l]<<" ";
						// }
						// cout<<"\n\n";
						
					}


					for (int l = 0; l<npop; l++) {
						i_nb = (i + cx[l] + (Nx+2))%(Nx+2);
						j_nb = (j + cy[l] + (Ny+2))%(Ny+2);
						k_nb = (k + cz[l] + (Nz+2))%(Nz+2);
						
						//Boundary condition : 
						if(fluid[i_nb][j_nb][k_nb].type == TYPE_F) {
							//Velocity BC : 
							fluid[i_nb][j_nb][k_nb].f[l] = fluid[i][j][k].f_star[l];

							//Thermal BC : 
							#ifdef ENERGY
								fluid[i_nb][j_nb][k_nb].g[l] = fluid[i][j][k].g_star[l];
							#endif
						}
						else if(fluid[i_nb][j_nb][k_nb].type == TYPE_S) {
							//Determining uwall 
							vector<double> uwall(3,0.0);
							Determine_Wall_Velocity(i, j, k, uwall);

							double uwall_dot_c = uwall[0]*cx[l] + uwall[1]*cy[l] + uwall[2]*cz[l];
							double uwall_dot_uwall = uwall[0]*uwall[0] + uwall[1]*uwall[1] + uwall[2]*uwall[2];
							fluid[i][j][k].f[opposite[l]] = fluid[i][j][k].f_star[l] - 2*w[l]*uwall_dot_c/cs2;

							#ifdef ENERGY
								if(fluid[i_nb][j_nb][k_nb].type_Energy == TYPE_T) {
									double g_eq = Calc_geq(i,j,k,l, fluid[i_nb][j_nb][k_nb].type, uwall);
									fluid[i][j][k].g[opposite[l]] = -fluid[i][j][k].g_star[l] + 2*g_eq;
								}
								else if(fluid[i_nb][j_nb][k_nb].type_Energy == TYPE_Q) {
									//Determine wall temperature 
									double Twall = Determine_Twall_Neumann(i, j, k, uwall, cdotn_pos, cdotn_notpos, j_dot_n);
									fluid[i_nb][j_nb][k_nb].T = Twall;

									//Finding the corresponding g at unknown direction (approximated as equilibrium)
									fluid[i][j][k].g[opposite[l]] = w[l]*fluid[i][j][k].density*Twall*(1 + uwall_dot_c/cs2 + pow(uwall_dot_c,2)/(2*cs2*cs2) - uwall_dot_uwall/(2*cs2));
								}	
							#endif
						}
					}
				}
			}
		}
	}



	//Boundary condition : 
	//Left (i = 0): 
	for(int j = 0; j<=Ny+1; j++) {
		for (int k = 0; k<=Nz+1; k++) {
			//Velocity boundary condition : 
			if(fluid[0][j][k].type == TYPE_OUT) {//for outflow
				fluid[0][j][k].f[1] = fluid[1][j][k].f[1];
				fluid[0][j][k].f[5] = fluid[1][j][k].f[5];
				fluid[0][j][k].f[8] = fluid[1][j][k].f[8];

				#ifdef ENERGY
					fluid[0][j][k].g[1] = fluid[1][j][k].g[1];
					fluid[0][j][k].g[5] = fluid[1][j][k].g[5];
					fluid[0][j][k].g[8] = fluid[1][j][k].g[8];
				#endif
			}
			// else if(fluid[0][j][k].type == TYPE_L) { //slip BC 
			// 	fluid[0][j][k].f[1] = fluid[0][j][k].f_star[3];
			// 	fluid[0][j][k].f[5] = fluid[0][j][k].f_star[6];
			// 	fluid[0][j][k].f[8] = fluid[0][j][k].f_star[7];

			// 	#ifdef ENERGY
			// 		fluid[0][j][k].g[1] = fluid[0][j][k].g_star[3];
			// 		fluid[0][j][k].g[5] = fluid[0][j][k].g_star[6];
			// 		fluid[0][j][k].g[8] = fluid[0][j][k].g_star[7];
			// 	#endif
			// }
			else if(fluid[0][j][k].type == TYPE_IN) { //inflow BC 
				double rho = fluid[0][j][k].density;
				double ux = fluid[0][j][k].ux;
				double uy = fluid[0][j][k].uy;
				double uz = fluid[0][j][k].uz;

				double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
				double c_dot_u[npop];

				for(int l = 0; l<npop; l++) {
					c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
					double f_eq = rho*w[l]*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
					fluid[0][j][k].f[l] = f_eq;
				}


				#ifdef ENERGY
					double T = fluid[0][j][k].T;

					for(int l = 0; l<npop; l++) {
						double g_eq = rho*w[l]*T*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[0][j][k].g[l] = g_eq;
					}
				#endif
			}
		}
	}

	//Right (i = Nx+1): 
	for(int j = 0; j<=Ny+1; j++) {
		for (int k = 0; k<=Nz+1; k++) {
			//Velocity boundary condition : 
			if(fluid[Nx+1][j][k].type == TYPE_OUT) { //for outflow
				fluid[Nx+1][j][k].f[3] = fluid[Nx][j][k].f[3];
				fluid[Nx+1][j][k].f[6] = fluid[Nx][j][k].f[6];
				fluid[Nx+1][j][k].f[7] = fluid[Nx][j][k].f[7];

				#ifdef ENERGY
					fluid[Nx+1][j][k].g[3] = fluid[Nx][j][k].g[3];
					fluid[Nx+1][j][k].g[6] = fluid[Nx][j][k].g[6];
					fluid[Nx+1][j][k].g[7] = fluid[Nx][j][k].g[7];
				#endif
			}
			// else if(fluid[Nx+1][j][k].type == TYPE_L) { //slip BC 
			// 	fluid[Nx+1][j][k].f[3] = fluid[Nx+1][j][k].f_star[1];
			// 	fluid[Nx+1][j][k].f[6] = fluid[Nx+1][j][k].f_star[5];
			// 	fluid[Nx+1][j][k].f[7] = fluid[Nx+1][j][k].f_star[8];

			// 	#ifdef ENERGY
			// 		fluid[Nx+1][j][k].g[3] = fluid[Nx+1][j][k].g_star[1];
			// 		fluid[Nx+1][j][k].g[6] = fluid[Nx+1][j][k].g_star[5];
			// 		fluid[Nx+1][j][k].g[7] = fluid[Nx+1][j][k].g_star[8];
			// 	#endif
			// }	
			else if(fluid[Nx+1][j][k].type == TYPE_IN) { //inflow BC 
				double rho = fluid[Nx+1][j][k].density;
				double ux = fluid[Nx+1][j][k].ux;
				double uy = fluid[Nx+1][j][k].uy;
				double uz = fluid[Nx+1][j][k].uz;

				double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
				double c_dot_u[npop];

				for(int l = 0; l<npop; l++) {
					c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
					double f_eq = rho*w[l]*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
					fluid[Nx+1][j][k].f[l] = f_eq;
				}


				#ifdef ENERGY
					double T = fluid[Nx+1][j][k].T;

					for(int l = 0; l<npop; l++) {
						double g_eq = rho*w[l]*T*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[Nx+1][j][k].g[l] = g_eq;
					}
				#endif
			}
		}
	}


	//Top (j = Ny+1):
	for(int i = 0; i<=Nx+1; i++) {
		for (int k = 0; k<=Nz+1; k++) {
			//Velocity boundary condition : 
			if(fluid[i][Ny+1][k].type == TYPE_OUT) {//for outflow
				fluid[i][Ny+1][k].f[4] = fluid[i][Ny][k].f[4];
				fluid[i][Ny+1][k].f[7] = fluid[i][Ny][k].f[7];
				fluid[i][Ny+1][k].f[8] = fluid[i][Ny][k].f[8];

				#ifdef ENERGY
					fluid[i][Ny+1][k].g[4] = fluid[i][Ny][k].g[4];
					fluid[i][Ny+1][k].g[7] = fluid[i][Ny][k].g[7];
					fluid[i][Ny+1][k].g[8] = fluid[i][Ny][k].g[8];
				#endif
			}
			// else if(fluid[i][Ny+1][k].type == TYPE_L) { //slip BC 
			// 	fluid[i][Ny+1][k].f[4] = fluid[i][Ny+1][k].f_star[2];
			// 	fluid[i][Ny+1][k].f[7] = fluid[i][Ny+1][k].f_star[6];
			// 	fluid[i][Ny+1][k].f[8] = fluid[i][Ny+1][k].f_star[5];

			// 	#ifdef ENERGY
			// 		fluid[i][Ny+1][k].g[4] = fluid[i][Ny+1][k].g_star[2];
			// 		fluid[i][Ny+1][k].g[7] = fluid[i][Ny+1][k].g_star[6];
			// 		fluid[i][Ny+1][k].g[8] = fluid[i][Ny+1][k].g_star[5];
			// 	#endif
			// }
			else if(fluid[i][Ny+1][k].type == TYPE_IN) { //inflow BC 
				double rho = fluid[i][Ny+1][k].density;
				double ux = fluid[i][Ny+1][k].ux;
				double uy = fluid[i][Ny+1][k].uy;
				double uz = fluid[i][Ny+1][k].uz;

				double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
				double c_dot_u[npop];

				for(int l = 0; l<npop; l++) {
					c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
					double f_eq = rho*w[l]*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
					fluid[i][Ny+1][k].f[l] = f_eq;
				}


				#ifdef ENERGY
					double T = fluid[i][Ny+1][k].T;

					for(int l = 0; l<npop; l++) {
						double g_eq = rho*w[l]*T*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[i][Ny+1][k].g[l] = g_eq;
					}
				#endif
			}
		}
	}


	//Bottom (j = 0): 
	for(int i = 0; i<=Nx+1; i++) {
		for (int k = 0; k<=Nz+1; k++) {
			//Velocity boundary condition : 
			if(fluid[i][0][k].type == TYPE_OUT) {//for outflow
				fluid[i][0][k].f[2] = fluid[i][1][k].f[2];
				fluid[i][0][k].f[5] = fluid[i][1][k].f[5];
				fluid[i][0][k].f[6] = fluid[i][1][k].f[6];
				
				#ifdef ENERGY
					fluid[i][0][k].g[2] = fluid[i][1][k].g[2];
					fluid[i][0][k].g[5] = fluid[i][1][k].g[5];
					fluid[i][0][k].g[6] = fluid[i][1][k].g[6];
				#endif
			}
			// else if(fluid[i][0][k].type == TYPE_L) { //slip BC 
			// 	fluid[i][0][k].f[2] = fluid[i][0][k].f_star[4];
			// 	fluid[i][0][k].f[5] = fluid[i][0][k].f_star[8];
			// 	fluid[i][0][k].f[6] = fluid[i][0][k].f_star[7];

			// 	#ifdef ENERGY
			// 		fluid[i][0][k].g[2] = fluid[i][0][k].g_star[4];
			// 		fluid[i][0][k].g[5] = fluid[i][0][k].g_star[8];
			// 		fluid[i][0][k].g[6] = fluid[i][0][k].g_star[7];
			// 	#endif
			// }
			else if(fluid[i][0][k].type == TYPE_IN) { //inflow BC 
				double rho = fluid[i][0][k].density;
				double ux = fluid[i][0][k].ux;
				double uy = fluid[i][0][k].uy;
				double uz = fluid[i][0][k].uz;

				double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
				double c_dot_u[npop];

				for(int l = 0; l<npop; l++) {
					c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
					double f_eq = rho*w[l]*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
					fluid[i][0][k].f[l] = f_eq;
				}


				#ifdef ENERGY
					double T = fluid[i][0][k].T;

					for(int l = 0; l<npop; l++) {
						double g_eq = rho*w[l]*T*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[i][0][k].g[l] = g_eq;
					}
				#endif
			}
		}
	}

	//Front (k = Nz+1) :
	for(int i = 0; i<=Nx+1; i++) {
		for (int j = 0; j<=Ny+1; j++) {
			//Velocity boundary condition : 
			if(fluid[i][j][Nz+1].type == TYPE_OUT) {//for outflow
				fluid[i][j][Nz+1].f[4] = fluid[i][j][Nz].f[4];
				fluid[i][j][Nz+1].f[7] = fluid[i][j][Nz].f[7];
				fluid[i][j][Nz+1].f[8] = fluid[i][j][Nz].f[8];

				#ifdef ENERGY
					fluid[i][j][Nz+1].g[4] = fluid[i][j][Nz].g[4];
					fluid[i][j][Nz+1].g[7] = fluid[i][j][Nz].g[7];
					fluid[i][j][Nz+1].g[8] = fluid[i][j][Nz].g[8];
				#endif
			}
			// else if(fluid[i][j][Nz+1].type == TYPE_L) { //slip BC 
			// 	fluid[i][j][Nz+1].f[4] = fluid[i][j][Nz+1].f_star[2];
			// 	fluid[i][j][Nz+1].f[7] = fluid[i][j][Nz+1].f_star[6];
			// 	fluid[i][j][Nz+1].f[8] = fluid[i][j][Nz+1].f_star[5];

			// 	#ifdef ENERGY
			// 		fluid[i][j][Nz+1].g[4] = fluid[i][j][Nz+1].g_star[2];
			// 		fluid[i][j][Nz+1].g[7] = fluid[i][j][Nz+1].g_star[6];
			// 		fluid[i][j][Nz+1].g[8] = fluid[i][j][Nz+1].g_star[5];
			// 	#endif
			// }
			else if(fluid[i][j][Nz+1].type == TYPE_IN) { //inflow BC 
				double rho = fluid[i][j][Nz+1].density;
				double ux = fluid[i][j][Nz+1].ux;
				double uy = fluid[i][j][Nz+1].uy;
				double uz = fluid[i][j][Nz+1].uz;

				double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
				double c_dot_u[npop];

				for(int l = 0; l<npop; l++) {
					c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
					double f_eq = rho*w[l]*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
					fluid[i][j][Nz+1].f[l] = f_eq;
				}


				#ifdef ENERGY
					double T = fluid[i][j][Nz+1].T;

					for(int l = 0; l<npop; l++) {
						double g_eq = rho*w[l]*T*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[i][j][Nz+1].g[l] = g_eq;
					}
				#endif
			}
		}
	}

	//Back (k = 0) :
	for(int i = 0; i<=Nx+1; i++) {
		for (int j = 0; j<=Ny+1; j++) {
			//Velocity boundary condition : 
			if(fluid[i][j][0].type == TYPE_OUT) {//for outflow
				fluid[i][j][0].f[4] = fluid[i][j][1].f[4];
				fluid[i][j][0].f[7] = fluid[i][j][1].f[7];
				fluid[i][j][0].f[8] = fluid[i][j][1].f[8];

				#ifdef ENERGY
					fluid[i][j][0].g[4] = fluid[i][j][1].g[4];
					fluid[i][j][0].g[7] = fluid[i][j][1].g[7];
					fluid[i][j][0].g[8] = fluid[i][j][1].g[8];
				#endif
			}
			// else if(fluid[i][j][0].type == TYPE_L) { //slip BC 
			// 	fluid[i][j][0].f[4] = fluid[i][j][1].f_star[2];
			// 	fluid[i][j][0].f[7] = fluid[i][j][1].f_star[6];
			// 	fluid[i][j][0].f[8] = fluid[i][j][1].f_star[5];

			// 	#ifdef ENERGY
			// 		fluid[i][j][0].g[4] = fluid[i][j][1].g_star[2];
			// 		fluid[i][j][0].g[7] = fluid[i][j][1].g_star[6];
			// 		fluid[i][j][0].g[8] = fluid[i][j][1].g_star[5];
			// 	#endif
			// }
			else if(fluid[i][j][0].type == TYPE_IN) { //inflow BC 
				double rho = fluid[i][j][0].density;
				double ux = fluid[i][j][0].ux;
				double uy = fluid[i][j][0].uy;
				double uz = fluid[i][j][0].uz;

				double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
				double c_dot_u[npop];

				for(int l = 0; l<npop; l++) {
					c_dot_u[l] = ux*cx[l] + uy*cy[l] + uz*cz[l];
					double f_eq = rho*w[l]*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
					fluid[i][j][Nz+1].f[l] = f_eq;
				}


				#ifdef ENERGY
					double T = fluid[i][j][0].T;

					for(int l = 0; l<npop; l++) {
						double g_eq = rho*w[l]*T*(1 + c_dot_u[l]/cs2 + pow(c_dot_u[l],2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[i][j][Nz+1].g[l] = g_eq;
					}
				#endif
			}
		}
	}
	
}

double LBM::Calc_feq(int i, int j, int k, int l) {
	double u_dot_u = pow(fluid[i][j][k].ux,2) + pow(fluid[i][j][k].uy,2) + pow(fluid[i][j][k].uz,2);
	double c_dot_u = fluid[i][j][k].ux*cx[l] + fluid[i][j][k].uy*cy[l] + fluid[i][j][k].uz*cz[l];
	double f_eq = fluid[i][j][k].density*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2) + pow(c_dot_u,3)/(6*pow(cs2,3)) - c_dot_u*u_dot_u/(2*cs2*cs2));	
	return f_eq;
}

#ifdef ENERGY
	double LBM::Calc_geq(int i, int j, int k, int l, short type, vector<double> uwall) {
		double g_eq;
		if(type == TYPE_S) {
			//Determine Twall :
			double Twall = Determine_Wall_Temperature(i, j, k);
			double u_dot_u = pow(uwall[0],2) + pow(uwall[1],2) + pow(uwall[2],2);
			double c_dot_u = uwall[0]*cx[l] + uwall[1]*cy[l] + uwall[2]*cz[l];
			g_eq = fluid[i][j][k].density*w[l]*Twall*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2) + pow(c_dot_u,3)/(6*pow(cs2,3)) - c_dot_u*u_dot_u/(2*cs2*cs2));	
		}

		else {
			double u_dot_u = pow(fluid[i][j][k].ux,2) + pow(fluid[i][j][k].uy,2) + pow(fluid[i][j][k].uz,2);
			double c_dot_u = fluid[i][j][k].ux*cx[l] + fluid[i][j][k].uy*cy[l] + fluid[i][j][k].uz*cz[l];
			g_eq = fluid[i][j][k].density*w[l]*fluid[i][j][k].T*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2) + pow(c_dot_u,3)/(6*pow(cs2,3)) - c_dot_u*u_dot_u/(2*cs2*cs2));	
		}
		
		return g_eq;
	}

#endif

void LBM::MacroProp() {
	#ifdef ENERGY
		double k_cond = k;
	#endif
	
	#pragma omp parallel for 
	for (int i = 0; i<=Nx+1; i++) {
		for(int j = 0; j<=Ny+1; j++) {
			for(int k = 0; k<=Nz+1; k++) {
				
				// if(fluid[i][j][k].type == TYPE_F or fluid[i][j][k].type == TYPE_IN or fluid[i][j][k].type == TYPE_OUT) {
					double rho = 0, rho_ux = 0, rho_uy = 0, rho_uz = 0;
					for(int l = 0; l<npop; l++) {
						rho += fluid[i][j][k].f[l];
						rho_ux += fluid[i][j][k].f[l]*cx[l];
						rho_uy += fluid[i][j][k].f[l]*cy[l];
						rho_uz += fluid[i][j][k].f[l]*cz[l];
					}
					
					fluid[i][j][k].density = rho;
					fluid[i][j][k].ux = (rho_ux + 0.5*dt*fluid[i][j][k].F[0])/rho;
					fluid[i][j][k].uy = (rho_uy + 0.5*dt*fluid[i][j][k].F[1])/rho;
					fluid[i][j][k].uz = (rho_uz + 0.5*dt*fluid[i][j][k].F[2])/rho;
					
					#if defined MORE_DETAILS
						//Pressure : 
						fluid[i][j][k].pressure = fluid[i][j][k].density/3.0;


						//Stress tensor: 
						for(int alphaa = 0; alphaa<3; alphaa++) {
                            for (int beta = alphaa; beta<3; beta++) {
                                double sum_sigma = 0;

                                //Calculating summation
                                for(int l = 0; l<npop; l++) {
                                    if(alphaa == 0 and beta == 0)      sum_sigma += fluid[i][j][k].f[l]*(cx[l] - fluid[i][j][k].ux)*(cx[l] - fluid[i][j][k].ux);
                                    else if(alphaa == 0 and beta == 1) sum_sigma += fluid[i][j][k].f[l]*(cx[l] - fluid[i][j][k].ux)*(cy[l] - fluid[i][j][k].uy);
                                    else if(alphaa == 0 and beta == 2) sum_sigma += fluid[i][j][k].f[l]*(cx[l] - fluid[i][j][k].ux)*(cz[l] - fluid[i][j][k].uz);
                                    else if(alphaa == 1 and beta == 1) sum_sigma += fluid[i][j][k].f[l]*(cy[l] - fluid[i][j][k].uy)*(cy[l] - fluid[i][j][k].uy);
                                    else if(alphaa == 1 and beta == 2) sum_sigma += fluid[i][j][k].f[l]*(cy[l] - fluid[i][j][k].uy)*(cz[l] - fluid[i][j][k].uz);
                                    else if(alphaa == 2 and beta == 2) sum_sigma += fluid[i][j][k].f[l]*(cz[l] - fluid[i][j][k].uz)*(cz[l] - fluid[i][j][k].uz);
                                }

                                //Calculating sigma (viscous stress)
                                fluid[i][j][k].sigma[alphaa][beta] = - (tau - 0.5)*sum_sigma/tau;
                                fluid[i][j][k].sigma[beta][alphaa] = fluid[i][j][k].sigma[alphaa][beta];								
                            }
                        }
						
						// //Vorticity :
						// double dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
						
						// //Derivative wrt x : 
						// if(i == 0) {
						// 	dudx = (-fluid[i+2][j][k].ux + 4*fluid[i+1][j][k].ux - 3*fluid[i][j][k].ux)/(2*dx);
						// 	dvdx = (-fluid[i+2][j][k].uy + 4*fluid[i+1][j][k].uy - 3*fluid[i][j][k].uy)/(2*dx);
						// 	dwdx = (-fluid[i+2][j][k].uz + 4*fluid[i+1][j][k].uz - 3*fluid[i][j][k].uz)/(2*dx);
						// }
						// else if (i == Nx+1) {
						// 	dudx = (fluid[i-2][j][k].ux - 4*fluid[i-1][j][k].ux + 3*fluid[i][j][k].ux)/(2*dx);
						// 	dvdx = (fluid[i-2][j][k].uy - 4*fluid[i-1][j][k].uy + 3*fluid[i][j][k].uy)/(2*dx);
						// 	dwdx = (fluid[i-2][j][k].uz - 4*fluid[i-1][j][k].uz + 3*fluid[i][j][k].uz)/(2*dx);
						// }
						// else {
						// 	dudx = (fluid[i+1][j][k].ux - fluid[i-1][j][k].ux)/(2*dx);
						// 	dvdx = (fluid[i+1][j][k].uy - fluid[i-1][j][k].uy)/(2*dx);
						// 	dwdx = (fluid[i+1][j][k].uz - fluid[i-1][j][k].uz)/(2*dx);
						// }

						// //Derivative wrt y : 
						// if(j == 0) {
						// 	dudy = (-fluid[i][j+2][k].ux + 4*fluid[i][j+1][k].ux - 3*fluid[i][j][k].ux)/(2*dy);
						// 	dvdy = (-fluid[i][j+2][k].uy + 4*fluid[i][j+1][k].uy - 3*fluid[i][j][k].uy)/(2*dy);
						// 	dwdy = (-fluid[i][j+2][k].uz + 4*fluid[i][j+1][k].uz - 3*fluid[i][j][k].uz)/(2*dy);
						// }
						// else if (j == Ny+1) {
						// 	dudy = (fluid[i][j-2][k].ux - 4*fluid[i][j-1][k].ux + 3*fluid[i][j][k].ux)/(2*dy);
						// 	dvdy = (fluid[i][j-2][k].uy - 4*fluid[i][j-1][k].uy + 3*fluid[i][j][k].uy)/(2*dy);
						// 	dwdy = (fluid[i][j-2][k].uz - 4*fluid[i][j-1][k].uz + 3*fluid[i][j][k].uz)/(2*dy);
						// }
						// else {
						// 	dudy = (fluid[i][j+1][k].ux - fluid[i][j-1][k].ux)/(2*dy);
						// 	dvdy = (fluid[i][j+1][k].uy - fluid[i][j-1][k].uy)/(2*dy);
						// 	dwdy = (fluid[i][j+1][k].uz - fluid[i][j-1][k].uz)/(2*dy);
						// }
						

						// //Derivative wrt z : 
						// if(ndim == 2) {
						// 	dudz = 0; dvdz = 0; dwdz = 0;
						// }
						// else if (ndim == 3) {
						// 	if(k == 0) {
						// 	dudz = (-fluid[i][j][k+2].ux + 4*fluid[i][j][k+1].ux - 3*fluid[i][j][k].ux)/(2*dz);
						// 	dvdz = (-fluid[i][j][k+2].uy + 4*fluid[i][j][k+1].uy - 3*fluid[i][j][k].uy)/(2*dz);
						// 	dwdz = (-fluid[i][j][k+2].uz + 4*fluid[i][j][k+1].uz - 3*fluid[i][j][k].uz)/(2*dz);
						// 	}
						// 	else if (k == Nz+1) {
						// 		dudz = (fluid[i][j][k-2].ux - 4*fluid[i][j][k-1].ux + 3*fluid[i][j][k].ux)/(2*dz);
						// 		dvdz = (fluid[i][j][k-2].uy - 4*fluid[i][j][k-1].uy + 3*fluid[i][j][k].uy)/(2*dz);
						// 		dwdz = (fluid[i][j][k-2].uz - 4*fluid[i][j][k-1].uz + 3*fluid[i][j][k].uz)/(2*dz);
						// 	}
						// 	else {
						// 		dudz = (fluid[i][j][k+1].ux - fluid[i][j][k-1].ux)/(2*dz);
						// 		dvdz = (fluid[i][j][k+1].uy - fluid[i][j][k-1].uy)/(2*dz);
						// 		dwdz = (fluid[i][j][k+1].uz - fluid[i][j][k-1].uz)/(2*dz);
						// 	}
						// }

						// fluid[i][j][k].vorticity[0] = dwdy - dvdz;
						// fluid[i][j][k].vorticity[1] = dudz - dwdx;
						// fluid[i][j][k].vorticity[2] = dvdx - dudy;
					#endif

					if(fluid[i][j][k].density <0 or fluid[i][j][k].density >100 ) {
						cout<<"Terminating in MacroProp function for density"<<endl;
						terminate();
					}

					#if defined ENERGY
						double T = 0;
						double sum_gx = 0, sum_gy = 0, sum_gz = 0;

						for (int l = 0; l<npop; l++) {
							T += fluid[i][j][k].g[l];
							sum_gx += fluid[i][j][k].g[l]*cx[l];
							sum_gy += fluid[i][j][k].g[l]*cy[l];
							sum_gz += fluid[i][j][k].g[l]*cz[l];
						}

						//For temperature : 
						fluid[i][j][k].T = T;
						fluid[i][j][k].q[0] = sum_gx;// - T*fluid[i][j][k].ux;
						fluid[i][j][k].q[1] = sum_gy;// - T*fluid[i][j][k].uy;
						fluid[i][j][k].q[2] = sum_gz;// - T*fluid[i][j][k].uz;
						
						if(fluid[i][j][k].T <0 or fluid[i][j][k].T >100 ) {
							cout<<"Terminating in MacroProp function for temperature"<<endl;
							terminate();
						}
					#endif

					//Updating force density term : 
					#if defined Poiseuille_Flow
						fluid[i][j][k].F[0] = -dp_dx;
					#elif defined Rayleigh_Benard
						double T_star = (Tb+Tt)/2;
						fluid[i][j][k].F[1] = g_beta*(fluid[i][j][k].T - T_star);
					#endif

				// }				
			}
		}
	}

	//For heat flux 
	// #ifdef ENERGY
	// 	for(int i = 0; i<=Nx+1; i++) {
	// 		for(int j = 0; j<=Ny+1; j++) {
	// 			for(int k = 0; k<=Nz+1; k++) {
	// 				if(fluid[i][j][k].type == TYPE_F or fluid[i][j][k].type == TYPE_IN or fluid[i][j][k].type == TYPE_OUT) {
	// 					double dT_dx, dT_dy, dT_dz;
	// 					//Calculate dT_dx
	// 					if(i == 0) {
	// 						dT_dx = (-3*fluid[i][j][k].T + 4*fluid[i+1][j][k].T - fluid[i+2][j][k].T)/(dx*2);
	// 					}
	// 					else if (i == Nx+1) {
	// 						dT_dx = (3*fluid[i][j][k].T - 4*fluid[i-1][j][k].T + fluid[i-2][j][k].T)/(dx*2);
	// 					}
	// 					else {
	// 						dT_dx = (fluid[i+1][j][k].T - fluid[i-1][j][k].T)/(dx*2);
	// 					}

	// 					//Calculate dT_dy
	// 					if(j == 0) {
	// 						dT_dy = (-3*fluid[i][j][k].T + 4*fluid[i][j+1][k].T - fluid[i][j+2][k].T)/(dy*2);
	// 					}
	// 					else if (j == Ny+1) {
	// 						dT_dy = (3*fluid[i][j][k].T - 4*fluid[i][j-1][k].T + fluid[i][j-2][k].T)/(dy*2);
	// 					}
	// 					else {
	// 						dT_dy = (fluid[i][j+1][k].T - fluid[i][j-1][k].T)/(dy*2);
	// 					}

	// 					//Calculate dT_dz 
	// 					if(k == 0) {
	// 						dT_dz = (-3*fluid[i][j][k].T + 4*fluid[i][j][k+1].T - fluid[i][j][k+2].T)/(dz*2);
	// 					}
	// 					else if (k == Nz+1) {
	// 						dT_dz = (3*fluid[i][j][k].T - 4*fluid[i][j][k-1].T + fluid[i][j][k-2].T)/(dz*2);
	// 					}
	// 					else {
	// 						dT_dz = (fluid[i][j][k+1].T - fluid[i][j][k-1].T)/(dz*2);
	// 					}

	// 					//Heat flux calculation
	// 					fluid[i][j][k].q[0] = -k_cond*dT_dx;
	// 					fluid[i][j][k].q[1] = -k_cond*dT_dy;
	// 					fluid[i][j][k].q[2] = -k_cond*dT_dz;
	// 				}
	// 			}
	// 		}
	// 	}
	// #endif
	
}

#if defined MORE_DETAILS
	void LBM::Sup_MacroProp(int n) {
		double told = (n-Nt)*dt;

		#pragma omp parallel for 
		for(int i = 0; i<=Nx+1; i++) {
			for(int j = 0; j<=Ny+1; j++) {
				for(int k = 0; k<=Nz+1; k++) {
					double ux = fluid[i][j][k].ux, uy = fluid[i][j][k].uy, uz = fluid[i][j][k].uz;
					double wx = fluid[i][j][k].vorticity[0], wy = fluid[i][j][k].vorticity[1], wz = fluid[i][j][k].vorticity[2];
					double p = fluid[i][j][k].pressure;
					//Average field properties : 
						//Velocity
						fluid[i][j][k].avg_velocity[0] = (fluid[i][j][k].avg_velocity[0]*told + ux*dt)/(told+dt);
						fluid[i][j][k].avg_velocity[1] = (fluid[i][j][k].avg_velocity[1]*told + uy*dt)/(told+dt);
						fluid[i][j][k].avg_velocity[2] = (fluid[i][j][k].avg_velocity[2]*told + uz*dt)/(told+dt);

						//Vorticity
						fluid[i][j][k].avg_vorticity[0] = (fluid[i][j][k].avg_vorticity[0]*told + wx*1)/(told+1);
						fluid[i][j][k].avg_vorticity[1] = (fluid[i][j][k].avg_vorticity[1]*told + wy*1)/(told+1);
						fluid[i][j][k].avg_vorticity[2] = (fluid[i][j][k].avg_vorticity[2]*told + wz*1)/(told+1);

						//pressure : 
						fluid[i][j][k].avg_pressure = (fluid[i][j][k].avg_pressure*told + p*1)/(told+1);

						//Reynolds stress : 
						double ux_avg = fluid[i][j][k].avg_velocity[0], uy_avg = fluid[i][j][k].avg_velocity[1], uz_avg = fluid[i][j][k].avg_velocity[2];
						fluid[i][j][k].avg_ReynoldStress[0][0] = (fluid[i][j][k].avg_ReynoldStress[0][0] * told + (ux - ux_avg)*(ux - ux_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[0][1] = (fluid[i][j][k].avg_ReynoldStress[0][1] * told + (ux - ux_avg)*(uy - uy_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[0][2] = (fluid[i][j][k].avg_ReynoldStress[0][2] * told + (ux - ux_avg)*(uz - uz_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[1][0] = (fluid[i][j][k].avg_ReynoldStress[1][0] * told + (uy - uy_avg)*(ux - ux_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[1][1] = (fluid[i][j][k].avg_ReynoldStress[1][1] * told + (uy - uy_avg)*(uy - uy_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[1][2] = (fluid[i][j][k].avg_ReynoldStress[1][2] * told + (uy - uy_avg)*(uz - uz_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[2][0] = (fluid[i][j][k].avg_ReynoldStress[2][0] * told + (uz - uz_avg)*(ux - ux_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[2][1] = (fluid[i][j][k].avg_ReynoldStress[2][1] * told + (uz - uz_avg)*(uy - uy_avg)*1)/(told+1);
						fluid[i][j][k].avg_ReynoldStress[2][2] = (fluid[i][j][k].avg_ReynoldStress[2][2] * told + (uz - uz_avg)*(uz - uz_avg)*1)/(told+1);

					//Rms field properties : 
						//velocity : 
						fluid[i][j][k].rms_velocity[0] = sqrt( (pow(fluid[i][j][k].rms_velocity[0],2)*told + pow(ux - ux_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_velocity[1] = sqrt( (pow(fluid[i][j][k].rms_velocity[1],2)*told + pow(uy - uy_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_velocity[2] = sqrt( (pow(fluid[i][j][k].rms_velocity[2],2)*told + pow(uz - uz_avg,2)*1 )/(told+1));

						//Vorticity : 
						fluid[i][j][k].rms_vorticity[0] = sqrt( (pow(fluid[i][j][k].rms_vorticity[0],2)*told + pow(wx - fluid[i][j][k].avg_vorticity[0],2)*1 )/(told+1));
						fluid[i][j][k].rms_vorticity[1] = sqrt( (pow(fluid[i][j][k].rms_vorticity[1],2)*told + pow(wy - fluid[i][j][k].avg_vorticity[1],2)*1 )/(told+1));
						fluid[i][j][k].rms_vorticity[2] = sqrt( (pow(fluid[i][j][k].rms_vorticity[2],2)*told + pow(wz - fluid[i][j][k].avg_vorticity[2],2)*1 )/(told+1));

						//pressure : 
						fluid[i][j][k].rms_pressure = sqrt( (pow(fluid[i][j][k].pressure,2)*told + pow(p - fluid[i][j][k].avg_pressure,2)*1 )/(told+1));

						//Reynolds stress : 
						double uu = fluid[i][j][k].ux*fluid[i][j][k].ux, uv = fluid[i][j][k].ux*fluid[i][j][k].uy, uw = fluid[i][j][k].ux*fluid[i][j][k].uz;
						double vv = fluid[i][j][k].uy*fluid[i][j][k].uy, vw = fluid[i][j][k].uy*fluid[i][j][k].uz;
						double ww = fluid[i][j][k].uz*fluid[i][j][k].uz;

						double uu_avg = fluid[i][j][k].avg_ReynoldStress[0][0], uv_avg = fluid[i][j][k].avg_ReynoldStress[0][1], uw_avg = fluid[i][j][k].avg_ReynoldStress[0][2];
						double vu_avg = fluid[i][j][k].avg_ReynoldStress[1][0], vv_avg = fluid[i][j][k].avg_ReynoldStress[1][1], vw_avg = fluid[i][j][k].avg_ReynoldStress[1][2];
						double wu_avg = fluid[i][j][k].avg_ReynoldStress[2][0], wv_avg = fluid[i][j][k].avg_ReynoldStress[2][1], ww_avg = fluid[i][j][k].avg_ReynoldStress[2][2];

						fluid[i][j][k].rms_ReynoldStress[0][0] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[0][0],2)*told + pow(uu - uu_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[0][1] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[0][1],2)*told + pow(uv - uv_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[0][2] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[0][2],2)*told + pow(uw - uw_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[1][0] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[1][0],2)*told + pow(uv - vu_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[1][1] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[1][1],2)*told + pow(vv - vv_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[1][2] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[1][2],2)*told + pow(vw - vw_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[2][0] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[2][0],2)*told + pow(uw - wu_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[2][1] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[2][1],2)*told + pow(vw - wv_avg,2)*1 )/(told+1));
						fluid[i][j][k].rms_ReynoldStress[2][2] = sqrt( (pow(fluid[i][j][k].rms_ReynoldStress[2][2],2)*told + pow(ww - ww_avg,2)*1 )/(told+1));
					
				}
			}
		}
	}
#endif

void LBM::Determine_Wall_Velocity(int i, int j, int k, vector<double>& uwall) {
	int num = 0; 	//number of solid type around fluid
	for(int l = 0; l<npop; l++) {
		int i_nb = (i + cx[l] + (Nx+2))%(Nx+2);
		int j_nb = (j + cy[l] + (Ny+2))%(Ny+2);
		int k_nb = (k + cz[l] + (Nz+2))%(Nz+2);

		if(fluid[i_nb][j_nb][k_nb].type == TYPE_S) {
			if(fluid[i_nb][j_nb][k_nb].ux != 0 or fluid[i_nb][j_nb][k_nb].uy != 0 or fluid[i_nb][j_nb][k_nb].uz != 0) {
				uwall[0] += fluid[i_nb][j_nb][k_nb].ux;
				uwall[1] += fluid[i_nb][j_nb][k_nb].uy;
				uwall[2] += fluid[i_nb][j_nb][k_nb].uz;
				num++;
			}
		}
	}

	if(num != 0) {
		uwall[0] /= num; uwall[1] /= num; uwall[2] /= num;
	}
	
}

#ifdef ENERGY
	double LBM::Determine_Wall_Temperature(int i, int j, int k) {
		int num = 0; 	//number of solid type around fluid
		double Twall = 0;
		for(int l = 0; l<npop; l++) {
			int i_nb = (i + cx[l] + (Nx+2))%(Nx+2);
			int j_nb = (j + cy[l] + (Ny+2))%(Ny+2);
			int k_nb = (k + cz[l] + (Nz+2))%(Nz+2);

			if(fluid[i_nb][j_nb][k_nb].type == TYPE_S) {
				Twall += fluid[i_nb][j_nb][k_nb].T;
				num++;
			}
		}

		Twall /= num;
		return Twall;
	}

	double LBM::Determine_Twall_Neumann(int i, int j, int k, vector<double> uwall, vector<int> cdotn_pos, vector<int> cdotn_notpos, double j_dot_n) {
		//Determining unit normal vector using cdotn_pos
		vector<double> normal(3,0);
		for(int l = 0; l<cdotn_pos.size(); l++) {
			normal[0] += cx[cdotn_pos[l]]; normal[1] += cy[cdotn_pos[l]]; normal[2] += cz[cdotn_pos[l]];
			
		}
		normal[0] /= cdotn_pos.size(); normal[1] /= cdotn_pos.size(); normal[2] /= cdotn_pos.size();

		//Normalizing: 
		// double norm = sqrt( pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2) );
		// normal[0] /= norm; normal[1] /= norm; normal[2] /= norm;

		//For cdotn_notpos part : 
		double sum_cdotn_notpos = 0;
		for(int l = 0; l<cdotn_notpos.size(); l++) {
			sum_cdotn_notpos += fluid[i][j][k].g[cdotn_notpos[l]]*( cx[cdotn_notpos[l]]*normal[0] + cy[cdotn_notpos[l]]*normal[1] + cz[cdotn_notpos[l]]*normal[2] );
		}

		//For cdotn_pos part : 
		double uwall_dot_uwall = uwall[0]*uwall[0] + uwall[1]*uwall[1] + uwall[2]*uwall[2];
		double sum_cdotn_pos = 0;
		for(int l = 0; l<cdotn_pos.size(); l++) {
			double uwall_dot_c = uwall[0]*cx[cdotn_pos[l]] + uwall[1]*cy[cdotn_pos[l]] + uwall[2]*cz[cdotn_pos[l]];
			//Equilibrium value : 
			double eq = w[cdotn_pos[l]]*fluid[i][j][k].density*(1 + uwall_dot_c/cs2 + pow(uwall_dot_c,2)/(2*cs2*cs2) - uwall_dot_uwall/(2*cs2) + pow(uwall_dot_c,3)/(6*pow(cs2,3)) - uwall_dot_c*uwall_dot_uwall/(2*cs2*cs2));
			sum_cdotn_pos += eq*( cx[cdotn_pos[l]]*normal[0] + cy[cdotn_pos[l]]*normal[1] + cz[cdotn_pos[l]]*normal[2] );
		}

		//Wall temperature : 
		double numerator = j_dot_n - sum_cdotn_notpos;
		double denum = sum_cdotn_pos;
		double Twall = numerator/denum;
		
		//Applying wall temperature associated with the normal direction (the opposite of c vector will be used)
		for(int l = 0; l<npop; l++) {
			if(cx[opposite[l]] == normal[0] and cy[opposite[l]] == normal[1] and cz[opposite[l]] == normal[2]) {
				int i_nb = (i + cx[l] + (Nx+2))%(Nx+2);
				int j_nb = (j + cy[l] + (Ny+2))%(Ny+2);
				int k_nb = (k + cz[l] + (Nz+2))%(Nz+2);
				fluid[i_nb][j_nb][k_nb].T = Twall;
				break;
			}
		}
		return Twall;
	}
#endif

int LBM::getNx() const {
	return Nx;
}

int LBM::getNy() const {
	return Ny;
}

int LBM::getNz() const {
	return Nz;
}

double LBM::getTAU() const{
	return tau;
}


void MARKER::Normal_Marker(int n, double x1, double x2, double y1, double y2, double z1, double z2) {
	//Finding tangential t1 : 
	double deltax = x2-x1, deltay = y2-y1, deltaz = z2-z1;
	double norm = sqrt( pow(deltax,2) + pow(deltay,2) + pow(deltaz,2) );
	Vector3d t1(deltax/norm, deltay/norm, deltaz/norm);
	

	//Finding tangential t2 : 
	Vector3d t2;
	if(ndim == 2) {
		t2[0] = 0; t2[1] = 0; t2[2] = 1;
	}
	else {

	}

	//Finding normal 
	Vector3d normal2 = t1.cross(t2);
	normal[0] = normal2[0];
	normal[1] = normal2[1];
	normal[2] = normal2[2];
}

//------------------------------------------STEPS FOR IMMERSED BOUNDARY METHOD---------------------------------
#if defined IBM
	//Find nearest coordinate : 
	void NearestLatticeCoordinate(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int& i, int& j, int& k, vector<MARKER>& marker, int n) {
		double dot = 0;
		for(int x = xmin; x<=xmax; x++) {
			for(int y = ymin; y<=ymax; y++) {
				for(int z = zmin; z<=zmax; z++) {
					Vector3d vec (x - marker[n].x, y - marker[n].y, z - marker[n].z);
					Vector3d normal_vec (marker[n].normal[0], marker[n].normal[1], marker[n].normal[2]);
					double dot_trial = vec.dot(normal_vec);

					if(dot_trial>=dot and dot_trial>0) {
						i = x, j = y, k = z;
						dot = dot_trial;
						
					}
				}	
			}	
		}
		// cout<<"xmin = "<<xmin<<", xmax = "<<xmax<<" marker n = "<<n<<endl;
	}

	//Kernel : 
	double kernel(double r) {
		double d;
		if(abs(r)>= 0 and abs(r)<= 1) d = (3-2*abs(r) + sqrt(1 + 4*abs(r) - 4*pow(r,2)) )/8 ;
		else if (abs(r)>=1 and abs(r)<=2) d = (5-2*abs(r) - sqrt(-7 + 12*abs(r) - 4*pow(r,2)) )/8 ;
		else d = 0;

		return d;
	}

	//Calculate average density for added mass
	double Calc_AverageDensity(LBM& lb, vector<MARKER>& marker, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) {
		double rho = 0;
		int number = 0;

		for(int i = xmin; i<=xmax; i++) {
			for(int j = ymin; j<=ymax; j++) {
				for(int k = zmin; k<=zmax; k++) {
					//Condition for adding : 
					bool is_dot_negative = 1;
					for(int n = 0; n<number_nodes; n++) {
						Vector3d vec (i*dx - marker[n].x, j*dy - marker[n].y, k*dz - marker[n].z);
						Vector3d normal_vec (marker[n].normal[0], marker[n].normal[1], marker[n].normal[2]);
						double dot_product = vec.dot(normal_vec);
						if(dot_product > 0) {
							is_dot_negative = 0; n = number_nodes;
						}
					}
					
					if(is_dot_negative == 1) {
						rho += lb.fluid[i][j][k].density; number++;
					}
				}
			}
		}
		rho /= number;
		
		return rho;
	}

	//Steps for IBM (reference from Kang et al Disseration for IB-LBM): 
	void IBLBM(LBM& lb, vector<MARKER>& marker, vector<double>& ds, vector<double>& Cl, vector<double>& Cd, double V_CM[], double t) {
		double Drag = 0, Lift = 0;
		double V_CM_before[3] = {V_CM[0], V_CM[1], V_CM[2]};

		//Update velocity : 
		#if defined Cylinder_IBM
			V_CM[0] = 0; V_CM[1] = 0; V_CM[2] = 0;
		#elif defined Cylinder_Inline_Oscillating
			V_CM[0] = -u_max*cos(2*PI*frequency*t);
			V_CM[1] = 0; V_CM[2] = 0;
		#elif defined Cylinder_Transverse_Oscillation
			V_CM[0] = 0; V_CM[2] = 0;
			V_CM[1] = -2*PI*frequency_e*amplitude*sin(2*PI*frequency_e*t);
		#elif defined Cylinder_IBM_Moving_Arbie
			V_CM[0] = -u_max;
			V_CM[1] = 0; V_CM[2] = 0;
		#endif

		//Update position : 
		for(int n = 0; n<number_nodes; n++) {
			marker[n].x += (V_CM[0] + V_CM_before[0])*dt/2;
			marker[n].y += (V_CM[1] + V_CM_before[1])*dt/2;
			marker[n].z += (V_CM[2] + V_CM_before[2])*dt/2;
			marker[n].u_boundary_des[0] = (V_CM[0]);
			marker[n].u_boundary_des[1] = (V_CM[1]);
			marker[n].u_boundary_des[2] = (V_CM[2]);
		}
		center_x += (V_CM[0] + V_CM_before[0])*dt/2;
		center_y += (V_CM[1] + V_CM_before[1])*dt/2;

		int xmin = floor(center_x - D), xmax = ceil(center_x + D);
		int ymin = floor(center_y - D), ymax = ceil(center_y + D);
		int zmin = 1, zmax = 1;
		

		//IBM step : 
		#if defined Explicit_Diffuse_Interface_Scheme
			//Unforced lattice velocity : 
			#pragma omp parallel for 
			for(int i = 0; i<=Nx+1; i++) {
				for(int j = 0; j<=Ny+1; j++) {
					for(int k = 0; k<=Nz+1; k++) {
						if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
							double rho = 0, rho_u[3] = {0,0,0};

							for(int l = 0; l<npop; l++) {
								rho += lb.fluid[i][j][k].f[l];
								rho_u[0] += lb.fluid[i][j][k].f[l]*cx[l];
								rho_u[1] += lb.fluid[i][j][k].f[l]*cy[l];
								rho_u[2] += lb.fluid[i][j][k].f[l]*cz[l];
							}

							lb.fluid[i][j][k].density = rho;
							lb.fluid[i][j][k].ux = rho_u[0]/rho;
							lb.fluid[i][j][k].uy = rho_u[1]/rho;
							lb.fluid[i][j][k].uz = rho_u[2]/rho;
						}
					}
				}
			}

			//Unforced velocity interpolation on boundary : 
			#pragma omp parallel for 
			for(int n = 0; n<number_nodes; n++) {
				marker[n].u_boundary[0] = 0; marker[n].u_boundary[1] = 0; marker[n].u_boundary[2] = 0;

				for(int i = xmin; i<=xmax; i++) {
					for(int j = ymin; j<=ymax; j++) {
						for(int k = zmin; k<=zmax; k++) {
							if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
								double rx = (i*dx - marker[n].x)/dx, ry = (j*dy - marker[n].y)/dy, rz = (k*dz - marker[n].z)/dz;
								marker[n].u_boundary[0] += lb.fluid[i][j][k].ux*kernel(rx)*kernel(ry)*kernel(rz);
								marker[n].u_boundary[1] += lb.fluid[i][j][k].uy*kernel(rx)*kernel(ry)*kernel(rz);
								marker[n].u_boundary[2] += lb.fluid[i][j][k].uz*kernel(rx)*kernel(ry)*kernel(rz);
							}
						}
					}
				}
			}

			//Boundary force evaluation : 
			#pragma omp parallel for 
			for(int n = 0; n<number_nodes; n++) {
				//Finding nearest coordinate 
				int i = 0,j = 0,k = 1;
				#if defined Cylinder_IBM
					NearestLatticeCoordinate( floor(marker[n].x), ceil(marker[n].x), floor(marker[n].y), ceil(marker[n].y), 1, 1, i, j, k, marker, n);
				#endif

				//Finding rho
				double rho = lb.fluid[i][j][k].density;

				//Finding the boundary force
				marker[n].F_boundary[0] = 2*rho*(marker[n].u_boundary_des[0] - marker[n].u_boundary[0])/dt;
				marker[n].F_boundary[1] = 2*rho*(marker[n].u_boundary_des[1] - marker[n].u_boundary[1])/dt;
				marker[n].F_boundary[2] = 2*rho*(marker[n].u_boundary_des[2] - marker[n].u_boundary[2])/dt;
			}

			//Force distribution and updating forcing velocity in lattice : 
			#pragma omp parallel for 
			for(int i = xmin; i<=xmax; i++) {
				for(int j = ymin; j<=ymax; j++) {
					for(int k = zmin; k<=zmax; k++) {
						if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
							lb.fluid[i][j][k].F[0] = 0; lb.fluid[i][j][k].F[1] = 0; lb.fluid[i][j][k].F[2] = 0; 

							for(int n = 0; n<number_nodes; n++) {
								double rx = (i*dx - marker[n].x)/dx, ry = (j*dy - marker[n].y)/dy, rz = (k*dz - marker[n].z)/dz;
								//Force distribution on lattice : 
								lb.fluid[i][j][k].F[0] += marker[n].F_boundary[0]*ds[n]*kernel(rx)*kernel(ry)*kernel(rz);
								lb.fluid[i][j][k].F[1] += marker[n].F_boundary[1]*ds[n]*kernel(rx)*kernel(ry)*kernel(rz);
								lb.fluid[i][j][k].F[2] += marker[n].F_boundary[2]*ds[n]*kernel(rx)*kernel(ry)*kernel(rz);
							}
							//Update forced velocity
							lb.fluid[i][j][k].ux += dt*lb.fluid[i][j][k].F[0]/(2*lb.fluid[i][j][k].density);
							lb.fluid[i][j][k].uy += dt*lb.fluid[i][j][k].F[1]/(2*lb.fluid[i][j][k].density);
							lb.fluid[i][j][k].uz += dt*lb.fluid[i][j][k].F[2]/(2*lb.fluid[i][j][k].density);
							Drag += -lb.fluid[i][j][k].F[0]*pow(dx,2);
							Lift += -lb.fluid[i][j][k].F[1]*pow(dx,2);
						}
					}
				}
			}

		#elif defined Implicit_Diffuse_Interface_Scheme
			//At m = 0 : 
			//Unforced velocity
			#pragma omp parallel for 
			for(int i = 0; i<=Nx+1; i++) {
				for(int j = 0; j<=Ny+1; j++) {
					for(int k = 0; k<=Nz+1; k++) {
						if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
							double rho = 0, rho_u[3] = {0,0,0};

							for(int l = 0; l<npop; l++) {
								rho += lb.fluid[i][j][k].f[l];
								rho_u[0] += lb.fluid[i][j][k].f[l]*cx[l];
								rho_u[1] += lb.fluid[i][j][k].f[l]*cy[l];
								rho_u[2] += lb.fluid[i][j][k].f[l]*cz[l];
							}

							lb.fluid[i][j][k].density = rho;
							lb.fluid[i][j][k].ux = rho_u[0]/rho;
							lb.fluid[i][j][k].uy = rho_u[1]/rho;
							lb.fluid[i][j][k].uz = rho_u[2]/rho;

							lb.fluid[i][j][k].F[0] = 0; lb.fluid[i][j][k].F[1] = 0; lb.fluid[i][j][k].F[2] = 0; 
						}
					}
				}
			}
			
			//Unforced velocity interpolation on boundary : 
			#pragma omp parallel for 
			for(int n = 0; n<number_nodes; n++) {
				marker[n].u_boundary[0] = 0; marker[n].u_boundary[1] = 0; marker[n].u_boundary[2] = 0;

				for(int i = xmin; i<=xmax; i++) {
					for(int j = ymin; j<=ymax; j++) {
						for(int k = zmin; k<=zmax; k++) {
							if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
								double rx = (i*dx - marker[n].x)/dx, ry = (j*dy - marker[n].y)/dy, rz = (k*dz - marker[n].z)/dz;
								marker[n].u_boundary[0] += lb.fluid[i][j][k].ux*kernel(rx)*kernel(ry)*kernel(rz);
								marker[n].u_boundary[1] += lb.fluid[i][j][k].uy*kernel(rx)*kernel(ry)*kernel(rz);
								marker[n].u_boundary[2] += lb.fluid[i][j][k].uz*kernel(rx)*kernel(ry)*kernel(rz);

								// if(i == ceil(center_x + D/2) and j == int(center_y) and n == 0){
								// 	cout<<"Lattice velocity before = "<<lb.fluid[i][j][k].ux/u_max<<", "<<lb.fluid[i][j][k].uy/u_max<<", "<<lb.fluid[i][j][k].uz/u_max<<endl;
								// 	cout<<"Boundary velocity before = "<<marker[n].u_boundary[0]/u_max<<", "<<marker[n].u_boundary[1]/u_max<<", "<<marker[n].u_boundary[2]/u_max<<endl;
								// }		
							}
						}
					}
				}
			}

			for(int m = 1; m<=m_max; m++) {
				
				//Boundary force evaluation : 
				for(int n = 0; n<number_nodes;n++) {
					int i,j,k = 1;
					NearestLatticeCoordinate(xmin, xmax, ymin, ymax, zmin, zmax, i, j, k, marker, n);
					double rho = lb.fluid[i][j][k].density;
					//Finding the boundary force
					marker[n].F_boundary[0] = 2*rho*(marker[n].u_boundary_des[0] - marker[n].u_boundary[0])/dt;
					marker[n].F_boundary[1] = 2*rho*(marker[n].u_boundary_des[1] - marker[n].u_boundary[1])/dt;
					marker[n].F_boundary[2] = 2*rho*(marker[n].u_boundary_des[2] - marker[n].u_boundary[2])/dt;
					Drag += -marker[n].F_boundary[0]*dx*ds[n];
					Lift += -marker[n].F_boundary[1]*dy*ds[n];

					// if(n == 0) {
					// 	cout<<"marker "<<n<<" at trial "<<m<<", u_boundary = "<< marker[n].u_boundary[0] <<", "<<marker[n].u_boundary[1]<<", "<<marker[n].u_boundary[2]<<endl;
					// }
					
				}

				//Force distribution on lattice and updating velocity : 
				#pragma omp parallel for 
				for(int i = xmin; i<=xmax; i++) {
					for(int j = ymin; j<=ymax; j++) {
						for(int k = zmin; k<=zmax; k++) {
							if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
								lb.fluid[i][j][k].Fm[0] = 0; lb.fluid[i][j][k].Fm[1] = 0; lb.fluid[i][j][k].Fm[2] = 0;
								for(int n = 0; n<number_nodes; n++) {
									double rx = (i*dx - marker[n].x)/dx, ry = (j*dy - marker[n].y)/dy, rz = (k*dz - marker[n].z)/dz;
									//Force distribution on lattice : 
									lb.fluid[i][j][k].Fm[0] += marker[n].F_boundary[0]*ds[n]*kernel(rx)*kernel(ry)*kernel(rz);
									lb.fluid[i][j][k].Fm[1] += marker[n].F_boundary[1]*ds[n]*kernel(rx)*kernel(ry)*kernel(rz);
									lb.fluid[i][j][k].Fm[2] += marker[n].F_boundary[2]*ds[n]*kernel(rx)*kernel(ry)*kernel(rz);
								}

								//Update forced velocity
								lb.fluid[i][j][k].F[0] += -lb.fluid[i][j][k].Fm[0];
								lb.fluid[i][j][k].F[1] += -lb.fluid[i][j][k].Fm[1];
								lb.fluid[i][j][k].F[2] += -lb.fluid[i][j][k].Fm[2];
								lb.fluid[i][j][k].ux += dt*lb.fluid[i][j][k].Fm[0]/(2*lb.fluid[i][j][k].density);
								lb.fluid[i][j][k].uy += dt*lb.fluid[i][j][k].Fm[1]/(2*lb.fluid[i][j][k].density);
								lb.fluid[i][j][k].uz += dt*lb.fluid[i][j][k].Fm[2]/(2*lb.fluid[i][j][k].density);
								// Drag += -lb.fluid[i][j][k].Fm[0]*pow(dx,2);
								// Lift += -lb.fluid[i][j][k].Fm[1]*pow(dx,2);

								
							}
						}
					}
				}

				//Velocity interpolation :
				#pragma omp parallel for 
				for(int n = 0; n<number_nodes; n++) {
					marker[n].u_boundary[0] = 0; marker[n].u_boundary[1] = 0; marker[n].u_boundary[2] = 0;
					for(int i = xmin; i<=xmax; i++) {
						for(int j = ymin; j<=ymax; j++) {
							for(int k = zmin; k<=zmax; k++) {
								if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_IN or lb.fluid[i][j][k].type == TYPE_OUT) {
									double rx = (i*dx - marker[n].x)/dx, ry = (j*dy - marker[n].y)/dy, rz = (k*dz - marker[n].z)/dz;
									marker[n].u_boundary[0] += lb.fluid[i][j][k].ux*kernel(rx)*kernel(ry)*kernel(rz);
									marker[n].u_boundary[1] += lb.fluid[i][j][k].uy*kernel(rx)*kernel(ry)*kernel(rz);
									marker[n].u_boundary[2] += lb.fluid[i][j][k].uz*kernel(rx)*kernel(ry)*kernel(rz);

									// if(i == ceil(center_x + D/2) and j == int(center_y) and n == 0 and m == m_max){
									// 	cout<<"Lattice velocity after = "<<lb.fluid[i][j][k].ux/u_max<<", "<<lb.fluid[i][j][k].uy/u_max<<", "<<lb.fluid[i][j][k].uz/u_max<<endl;
									// 	cout<<"Boundary velocity after = "<<marker[n].u_boundary[0]/u_max<<", "<<marker[n].u_boundary[1]/u_max<<", "<<marker[n].u_boundary[2]/u_max<<endl<<endl;
									// }		
								}
							}
							
						}
					}
					// if(n == 0) {
					// 	cout<<"marker "<<n<<" at trial "<<m<<", u_boundary = "<< marker[n].u_boundary[0]/u_max<<", "<<marker[n].u_boundary[1]/u_max<<", "<<marker[n].u_boundary[2]/u_max<<endl;
					// }
				} 
				
			}
		#endif
		
		//Added mass contribution
		//Calculate average density 
		double rho = Calc_AverageDensity(lb, marker, xmin, xmax, ymin, ymax, zmin, zmax);

		//Calculate for added mass :
		Drag += rho*Vsolid*(V_CM[0] - V_CM_before[0])/dt;
		Lift += rho*Vsolid*(V_CM[1] - V_CM_before[1])/dt;

		//Lift and drag coefficient : 
		#if defined Cylinder_IBM or defined Cylinder_IBM_Moving_Arbie
			double CL = Lift/(0.5*rho0*pow(u_max,2)*D);
			double CD = Drag/(0.5*rho0*pow(u_max,2)*D);
			Cl.push_back(CL);
			Cd.push_back(CD);
		#elif defined Cylinder_Inline_Oscillating
			double CL = Lift/(0.5*rho0*pow(u_max,2)*(D+dx));
			double CD = Drag/(0.5*rho0*pow(u_max,2)*(D+dx));
			Cl.push_back(CL);
			Cd.push_back(CD);
		#elif defined Cylinder_Transverse_Oscillation
			double CL = Lift/(0.5*rho0*pow(u_max,2)*D);
			double CD = Drag/(0.5*rho0*pow(u_max,2)*D);
			Cl.push_back(CL);
			Cd.push_back(CD);
		#endif
	}
	
	
#endif
//-------------------------------------------------------------------------------------------------------------




//------------------------------------------------SETUP FOR SIMULATION ----------------------------------------


#if defined Cylinder_Basic
	void Generate_Cylinder(LBM& lb, int& count) {
		double radius = D/2.0;

		#pragma omp parallel for 
		for(int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 0; k<=Nz; k++) {
					
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					
					if(sqrt(pow(center_x-rx,2.0) + pow(center_y-ry,2.0)) <= radius) {
						lb.fluid[i][j][k].type = TYPE_S;
						count++;
					}

				}
			}	
		}
		cout<<"count = "<<count<<"and Nx*Ny = "<<Nx*Ny<<endl;
	}

	
	LBM main_setup_Cylinder() { //Cylinder 2D flow Bounce Back 
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;

		//generating cylinder and its solid content 
		Generate_Cylinder(lb, count);
		cout<<"Count = "<<count<<endl;
		cout<<"Cylinder generated\n";
		
		
		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_IN; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;

						#if defined MORE_DETAILS
							lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
						#endif
					}
					
					
				}
			}
		}
		
		return lb;
	}
#elif defined Poiseuille_Flow
	LBM main_setup_Poiseuille() {
		LBM lb(Nx, Ny, Nz, tau);

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(i == 0) lb.fluid[i][j][k].type = TYPE_F; //no boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_F; //no boundary					
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //No slip boundary

					#ifdef ENERGY
						if(j == 0 or j == Ny+1) {
							lb.fluid[i][j][k].type_Energy = TYPE_T;
							lb.fluid[i][0][k].T = Tb;
							lb.fluid[i][Ny+1][k].T = Tt;
						}
					#endif
					
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;	
						lb.fluid[i][j][k].F[0] = -dp_dx;			

						#if defined MORE_DETAILS
							lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
						#endif

						#ifdef ENERGY
							lb.fluid[i][j][k].T = Tb;
						#endif
					}
					
					
				}
			}
		}
		
		return lb;		
	}
#elif defined Couette_Flow
	LBM main_setup_Couette() {
		LBM lb(Nx, Ny, Nz, tau);

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(i == 0) lb.fluid[i][j][k].type = TYPE_F; //no boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_F; //no boundary									
					if( (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //no boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //No slip boundary	

					if(j == Ny+1) lb.fluid[i][j][k].ux = uwall_topx;
					if(j == 0) lb.fluid[i][j][k].ux = uwall_botx;

					#ifdef ENERGY
						if(j == 0 or j == Ny+1) {
							lb.fluid[i][j][k].type_Energy = TYPE_T;
							lb.fluid[i][0][k].T = Tb;
							lb.fluid[i][Ny+1][k].T = Tt;
						}
					#endif
					
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;

						

						#if defined MORE_DETAILS
							lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
						#endif

						#ifdef ENERGY
							lb.fluid[i][j][k].T = Tb;
						#endif
					}
					
					
				}
			}
		}
		
		return lb;		
	}
#elif defined Rayleigh_Benard
	LBM main_setup_RayleighBenard1() {
		LBM lb(Nx, Ny, Nz, tau);

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					#ifdef ENERGY
						//From incorpera
						if( k == 0 or  k == Nz+1 ) {
							lb.fluid[i][j][k].type = TYPE_S; //wall BC
							lb.fluid[i][j][k].type_Energy = TYPE_Q; //Neumann BC	
							lb.fluid[i][j][k].T = Tt;
						}
						if(i == 0 or i == Nx+1) {
							lb.fluid[i][j][k].type = TYPE_S; //wall BC
							lb.fluid[i][j][k].type_Energy = TYPE_Q; //Neumann BC	
							lb.fluid[i][j][k].T = Tt;
						}
						if(j == 0 or j == Ny+1) {
							lb.fluid[i][j][k].type = TYPE_S; //wall BC
							lb.fluid[i][j][k].type_Energy = TYPE_T; //Dirichlet BC	
							lb.fluid[i][j][k].T = Tt;
						}
						
						
						lb.fluid[i][Ny+1][k].T = Tb;
						lb.fluid[i][0][k].T = Tt;
						
						if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
							lb.fluid[i][j][k].density = 1.0;
							lb.fluid[i][j][k].ux = 0.0;
							lb.fluid[i][j][k].uy = 0.0;
							lb.fluid[i][j][k].uz = 0.0;

							

							#if defined MORE_DETAILS
								lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
							#endif

							lb.fluid[i][j][k].T = Tt;

							double T_star = (Tb+Tt)/2;
							lb.fluid[i][j][k].F[1] = g_beta*(lb.fluid[i][j][k].T - T_star);
						}
					#endif

					
					
					
				}
			}
		}
		
		return lb;	
	}

#elif defined Cylinder_Tandem_Circle
	void Generate_Circle(LBM& lb, int& count, int center_x, int center_y) {
		double radius = D/2;
		for(int i = center_x-D; i<=center_x+D; i++) {
			for (int j = center_y - D; j<=center_y + D; j++) {
				for (int k = 0; k<=Nz; k++) {
					
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					
					if(sqrt(pow(center_x-rx,2.0) + pow(center_y-ry,2.0)) <= radius) {
						lb.fluid[i][j][k].type = TYPE_S;
						count++;
					}

				}
			}	
		}

	}

	LBM main_setup_Tandem_Circle() {
		LBM lb(Nx, Ny, Nz, tau);
		int count = 0;

		//Generating Circle 
		//#pragma omp parallel for 
		for(int n = 0; n<N; n++) {
			Generate_Circle(lb, count, center_x, center_y - n*L);
		}

		cout<<"Count = "<<count<<endl;
		cout<<"Circles generated\n";

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_IN; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;

						#if defined MORE_DETAILS
							lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
						#endif
					}
					
					
				}
			}
		}

		return lb;
	}
#elif defined Cylinder_Tandem_Triangle
	void Generate_Triangle(LBM& lb, int& count, int tip_x, int tip_y) {

		for(int i = 0; i<=D; i++) {
			for(int j = floor(-i/sqrt(3)); j<=ceil(i/sqrt(3)); j++) {
				for(int k = 0; k<=Nz+1; k++) {
					lb.fluid[tip_x + i][tip_y + j][k].type = TYPE_S; count++;
				}
			}
		}

	}

	LBM main_setup_Tandem_Triangle() {
		LBM lb(Nx, Ny, Nz, tau);
		int count = 0;

		//Generating Triangle 
		//#pragma omp parallel for 
		for(int n = 0; n<N; n++) {
			Generate_Triangle(lb, count, tip_x + n*L, tip_y);
		}

		cout<<"Count = "<<count<<endl;
		cout<<"Triangles generated\n";

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_IN; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(k == 0 or k == Nz+1) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					
					lb.fluid[i][j][k].density = 1.0;
					if(lb.fluid[i][j][k].type == TYPE_S) { lb.fluid[i][j][k].ux = 0;}
					else {lb.fluid[i][j][k].ux = u_max;}
					
					lb.fluid[i][j][k].uy = 0.0;
					lb.fluid[i][j][k].uz = 0.0;

					#if defined MORE_DETAILS
						lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
					#endif
					
				}
			}
		}

		return lb;
	}
#elif defined Sphere_BB
	void Generate_Sphere(LBM& lb, int& count) {
		double radius = D/2.0;

		#pragma omp parallel for 
		for(int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 0; k<=Nz; k++) {
					
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					double rz = dz*k-dz*0.5;
					
					if(sqrt(pow(center_x-rx,2.0) + pow(center_y-ry,2.0) + pow(center_z - rz, 2.0)) <= radius) {
						lb.fluid[i][j][k].type = TYPE_S;
						count++;
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}

				}
			}	
		}
		cout<<"count = "<<count<<"and Nx*Ny*Nz = "<<Nx*Ny*Nz<<endl;
	}

	
	LBM main_setup_Sphere() { //Cylinder 2D flow Bounce Back 
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;

		//generating cylinder and its solid content 
		Generate_Sphere(lb, count);
		cout<<"Count = "<<count<<endl;
		cout<<"Sphere generated\n";
		
		
		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_IN; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary
					if(k == 0 or  k == Nz+1) lb.fluid[i][j][k].type = TYPE_OUT; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;

						#if defined MORE_DETAILS
							lb.fluid[i][j][k].pressure = lb.fluid[i][j][k].density/3.0;
						#endif
					}
					
					
					
				}
			}
		}
		
		return lb;
	}
#elif defined Cylinder_IBM
	LBM main_setup_Cylinder_IBM() { //Cylinder 2D IBM
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_IN; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
					
				}
			}
		}
		
		return lb;
	}

	void Generate_Cylinder_Marker(vector<MARKER>& marker) {
		double radius = D/2.0;
		const double PI = 3.14159265;

		//Setting the markers :
		for (int i = 0; i<number_nodes; i++) {
			//Consider 2D flow
			marker[i].x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; 
			marker[i].y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
			marker[i].z = 1;
		}

		//Calculate normal :
		for (int n = 0; n<number_nodes; n++) {
			int n1 = n, n2 = (n+1)%number_nodes;	//counterclockwise definition of markers
			double x1 = marker[n1].x, x2 = marker[n2].x;
			double y1 = marker[n1].y, y2 = marker[n2].y;
			double z1 = marker[n1].z, z2 = marker[n2].z;
			marker[n].Normal_Marker(n, x1, x2, y1, y2, z1, z2);
		}
	}
#elif defined Cylinder_Inline_Oscillating
	LBM main_setup_Cylinder_InlineOscillating() {
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_OUT; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
					
				}
			}
		}
		
		return lb;
	}

	void Generate_Cylinder_Marker(vector<MARKER>& marker) {
		double radius = D/2.0;
		const double PI = 3.14159265;

		//Setting the markers :
		for (int i = 0; i<number_nodes; i++) {
			//Consider 2D flow
			marker[i].x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; 
			marker[i].y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
			marker[i].z = 1;
		}

		//Calculate normal :
		for (int n = 0; n<number_nodes; n++) {
			int n1 = n, n2 = (n+1)%number_nodes;	//counterclockwise definition of markers
			double x1 = marker[n1].x, x2 = marker[n2].x;
			double y1 = marker[n1].y, y2 = marker[n2].y;
			double z1 = marker[n1].z, z2 = marker[n2].z;
			marker[n].Normal_Marker(n, x1, x2, y1, y2, z1, z2);
		}
	}
#elif defined Cylinder_Transverse_Oscillation
	LBM main_setup_Cylinder_TransverseOscillating() {
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_IN; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
					
				}
			}
		}
		
		return lb;
	}

	void Generate_Cylinder_Marker(vector<MARKER>& marker) {
		double radius = D/2.0;
		const double PI = 3.14159265;

		//Setting the markers :
		for (int i = 0; i<number_nodes; i++) {
			//Consider 2D flow
			marker[i].x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; 
			marker[i].y = radius*sin(2*PI/number_nodes*i) + center_y + amplitude + 0.5*dy;
			marker[i].z = 1;
		}

		//Calculate normal :
		for (int n = 0; n<number_nodes; n++) {
			int n1 = n, n2 = (n+1)%number_nodes;	//counterclockwise definition of markers
			double x1 = marker[n1].x, x2 = marker[n2].x;
			double y1 = marker[n1].y, y2 = marker[n2].y;
			double z1 = marker[n1].z, z2 = marker[n2].z;
			marker[n].Normal_Marker(n, x1, x2, y1, y2, z1, z2);
		}
	}
#elif defined Cylinder_IBM_Moving_Arbie
	LBM main_setup_Cylinder_IBM_Moving() {
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {					
					if(i == 0) lb.fluid[i][j][k].type = TYPE_OUT; //Inlet boundary
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_OUT; //Outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_OUT; //Slip boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_IN || lb.fluid[i][j][k].type == TYPE_OUT || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
					
				}
			}
		}
		
		return lb;
	}

	void Generate_Cylinder_Marker(vector<MARKER>& marker) {
		double radius = D/2.0;
		const double PI = 3.14159265;

		//Setting the markers :
		for (int i = 0; i<number_nodes; i++) {
			//Consider 2D flow
			marker[i].x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; 
			marker[i].y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
			marker[i].z = 1;
		}

		//Calculate normal :
		for (int n = 0; n<number_nodes; n++) {
			int n1 = n, n2 = (n+1)%number_nodes;	//counterclockwise definition of markers
			double x1 = marker[n1].x, x2 = marker[n2].x;
			double y1 = marker[n1].y, y2 = marker[n2].y;
			double z1 = marker[n1].z, z2 = marker[n2].z;
			marker[n].Normal_Marker(n, x1, x2, y1, y2, z1, z2);
		}
	}
#endif

//---------------------------------------------------------------------------------------------------------------
//usual function definition
void print_Info() {
	#if defined Cylinder_Basic
		cout<<"Case                                         : Flow over cylinder (Bounce back)"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Freestream Velocity                          : "<<u_max<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Cylinder geometry                            : Circle"<<endl;
		cout<<"Diameter                                     : "<<D<<endl;
		cout<<"center x                                     : "<<center_x<<endl;
		cout<<"center y                                     : "<<center_y<<endl;
    #elif defined Poiseuille_Flow
		cout<<"Case                                         : Poiseuille Flow 2D"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Pressure_gradient                            : "<<dp_dx<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;        
        cout<<endl;

		#ifdef ENERGY
			cout<<"Bottom wall temperature                      : "<<Tb<<endl;
			cout<<"Top wall temperature                         : "<<Tt<<endl;

			cout<<"Prandtl number                               : "<<Pr<<endl;
			
			cout<<"Thermal diffusivity                          : "<<alpha<<endl;
			cout<<"Thermal relaxation time                      : "<<tau_g<<endl;
			cout<<"Thermal conductivity                         : "<<k<<endl;
			cout<<"Specific heat constant                       : "<<cp<<endl;
		#endif
	#elif defined Couette_Flow
		cout<<"Case                                         : Couette Flow 2D"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Wall Velocity                                : "<<u_wall<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;        
        cout<<endl;

		#ifdef ENERGY
			cout<<"Bottom wall temperature                      : "<<Tb<<endl;
			cout<<"Wall temperature difference                  : "<<dT<<endl;
			cout<<"Top wall temperature                         : "<<Tt<<endl;

			cout<<"Prandtl number                               : "<<Pr<<endl;
			cout<<"Brinkmann number                             : "<<Br<<endl;
			cout<<"Eckart number                                : "<<Ec<<endl;
			
			cout<<"Thermal diffusivity                          : "<<alpha<<endl;
			cout<<"Thermal relaxation time                      : "<<tau_g<<endl;
			cout<<"Thermal conductivity                         : "<<k<<endl;
			cout<<"Specific heat constant                       : "<<cp<<endl;
		#endif
	#elif defined Rayleigh_Benard
		#ifdef ENERGY
			cout<<"Case                                         : Rayleigh-Benard Flow"<<endl<<endl;

			cout<<"Simulation Parameters : "<<endl;
			cout<<"tau                                          : "<<tau<<endl;
			cout<<"nu                                           : "<<nu<<endl;
			cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;        
			cout<<endl;

			cout<<"Bottom wall temperature                      : "<<Tb<<endl;
			cout<<"Wall temperature difference                  : "<<dT<<endl;
			cout<<"Top wall temperature                         : "<<Tt<<endl;

			cout<<"Prandtl number                               : "<<Pr<<endl;
			cout<<"Rayleigh number                              : "<<Ra<<endl;
			cout<<"Thermal diffusivity                          : "<<alpha<<endl;
			cout<<"Thermal relaxation time                      : "<<tau_g<<endl;
			cout<<"Thermal conductivity                         : "<<k<<endl;
			cout<<"Specific heat constant                       : "<<cp<<endl;
		#endif
	#elif defined Cylinder_Tandem_Circle
		cout<<"Case                                         : Tandem Circular Cylinder"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Freestream Velocity                          : "<<u_max<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Cylinder geometry                            : Circle"<<endl;
		cout<<"Horizontal distance of triangle              : "<<D<<endl;
		cout<<"Tip coordinate (x)                           : "<<center_x<<endl;
		cout<<"Tip coordinate (y)                           : "<<center_y<<endl;
		cout<<"Spacing to horizontal ratio                  : "<<LpD<<endl;


	#elif defined Cylinder_Tandem_Triangle
		cout<<"Case                                         : Equilateral Triangle Tandem Cylinder"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Freestream Velocity                          : "<<u_max<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Cylinder geometry                            : Equilateral Triangle"<<endl;
		cout<<"Horizontal distance of triangle              : "<<D<<endl;
		cout<<"Tip coordinate (x)                           : "<<tip_x<<endl;
		cout<<"Tip coordinate (y)                           : "<<tip_y<<endl;
		cout<<"Spacing to horizontal ratio                  : "<<LpD<<endl;
	#elif defined Sphere_BB
		cout<<"Case                                         : Flow over sphere (Bounce back)"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Freestream Velocity                          : "<<u_max<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Geometry                                     : Sphere"<<endl;
		cout<<"Diameter                                     : "<<D<<endl;
		cout<<"center x                                     : "<<center_x<<endl;
		cout<<"center y                                     : "<<center_y<<endl;
		cout<<"center z                                     : "<<center_z<<endl;
	#elif defined Cylinder_IBM
		cout<<"Case                                         : Cylinder IBM"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Freestream Velocity                          : "<<u_max<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Cylinder geometry                            : Circle"<<endl;
		cout<<"Diameter                                     : "<<D<<endl;
		cout<<"center x                                     : "<<center_x<<endl;
		cout<<"center y                                     : "<<center_y<<endl<<endl;

		cout<<"IBM properties                               : \n";
		cout<<"Stencil                                      : 4 \n";
		
		
	#elif defined Cylinder_Inline_Oscillating
		cout<<"Case                                         : Cylinder Inline Oscillating"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
		cout<<"Keulegan-Carpenter number                    : "<<KC<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Cylinder geometry                            : Circle"<<endl;
		cout<<"Diameter                                     : "<<D<<endl;
		cout<<"center x                                     : "<<center_x<<endl;
		cout<<"center y                                     : "<<center_y<<endl<<endl;

		cout<<"External motion                              : \n";
		cout<<"Maximum object velocity                      : "<<u_max<<endl;
		cout<<"Period of oscillation                        : "<<1/frequency<<endl;
		cout<<"Frequency of movement (f)                    : "<<frequency<<endl;
		cout<<"Amplitude of movement (A)                    : "<<amplitude<<endl<<endl;
		

		cout<<"IBM properties                               : \n";
		cout<<"Stencil                                      : 4 \n";
	#elif defined Cylinder_Transverse_Oscillation
		cout<<"Case                                         : Cylinder Inline Oscillating"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      

		cout<<"Cylinder geometry                            : Circle"<<endl;
		cout<<"Diameter                                     : "<<D<<endl;
		cout<<"center x                                     : "<<center_x<<endl;
		cout<<"center y                                     : "<<center_y<<endl<<endl;

		cout<<"External motion                              : \n";
		cout<<"Maximum object velocity                      : "<<u_max<<endl;
		cout<<"Frequency ratio (fe/f0)                      : "<<frequency_ratio<<endl;
		cout<<"Strouhal number for f0                       : "<<Strouhal_number<<endl;
		cout<<"Period of oscillation                        : "<<1/frequency_e<<endl;
		cout<<"Frequency of movement (fe)                   : "<<frequency_e<<endl;
		cout<<"Amplitude ratio (A/d)                        : "<<ratio_Amd<<endl;
		cout<<"Amplitude of movement (A)                    : "<<amplitude<<endl<<endl;
		

		cout<<"IBM properties                               : \n";
		cout<<"Stencil                                      : 4 \n";
	#elif defined Cylinder_IBM_Moving_Arbie
		cout<<"Case                                         : Cylinder IBM (Benchmark)"<<endl<<endl;

        cout<<"Simulation Parameters : "<<endl;
        cout<<"Reynolds number                              : "<<Re<<endl;
        cout<<"tau                                          : "<<tau<<endl;
        cout<<"nu                                           : "<<nu<<endl;
        cout<<"Freestream Velocity                          : "<<u_max<<endl;
		cout<<"Target nondimensional time                   : "<<target_tstar<<endl;  
        cout<<"Total simulation time                        : "<<T_OUT<<endl<<endl;      
		

		cout<<"Cylinder geometry                            : Circle"<<endl;
		cout<<"Diameter                                     : "<<D<<endl;
		cout<<"center x                                     : "<<center_x<<endl;
		cout<<"center y                                     : "<<center_y<<endl<<endl;

		cout<<"IBM properties                               : \n";
		cout<<"Stencil                                      : 4 \n";
	#endif

	#ifdef IBM
		#if defined Explicit_Diffuse_Interface_Scheme
			cout<<"Method                                       : Explicit Diffuse Interface Scheme \n";
		#elif defined Implicit_Diffuse_Interface_Scheme
			cout<<"Method                                       : Implicit Diffuse Interface Scheme \n";
		#endif
	#endif
	
	cout<<"\nAdditional Mode                              : "<<endl;

	#if defined MORE_DETAILS
		cout<<"More details (pressure and viscous)          : ON\n";
	#else 
		cout<<"More details (pressure and viscous)          : OFF\n";
	#endif

	#ifdef ENERGY
		cout<<"Energy                                       : ON\n";
	#else
		cout<<"Energy                                       : OFF\n";
	#endif

	#if defined Dimensional
		cout<<"Dimensional                                  : ON\n\n";
	#else 
		cout<<"Dimensional                                  : OFF\n\n";
	#endif

	

	cout<<"Output result every "<<Nt<<" time step \n\n";
    cout<<"----------------------------------------------------------------------------------------------------"<<endl;
}

//Output fluid properties
void OutputVTK(int &nout, LBM &lb) {
	int		i,j,k;
	char	filename[128];
	FILE	*fp;
	unsigned int array_size;
	unsigned long int offset=0;
	// short  num16; // Int16 2byte
	float  val32; // Float32 4byte

	#if defined Cylinder_Tandem_Triangle
		sprintf(filename,"./output Tandem Triangle/LpD = %.1f Re = %.0f/Flowfield/field%06d.vtr",LpD, Re, nout);
	#else
		sprintf(filename,"./field%06d.vtr",nout);
	#endif
	
	fp=fopen(filename,"wb");

	#if defined D2Q9
		int zcomp = Nz;	
	#else 
		int zcomp = Nz;
	#endif
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
	fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",Nx+2,Ny+2,zcomp);
	fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Nx+2,Ny+2,zcomp);
	fprintf(fp,"    <PointData>\n");
	fprintf(fp,"    </PointData>\n");
	fprintf(fp,"    <CellData>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);

	#if defined Couette_Flow or defined Poiseuille_Flow
		#if defined Velocity_Ratio
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity Ratio\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
		#elif defined Dimensional
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
		#else
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
		#endif
	#else
		#ifndef Dimensional
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
		#else 
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
		#endif
	#endif
	
	#ifdef IBM
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Force density\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
	#endif
	

	#if defined MORE_DETAILS
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Pressure\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Viscous Stress\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
	#endif

	#ifdef ENERGY
		#if defined Couette_Flow or defined Poiseuille_Flow
			#if defined Temperature_Ratio
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Temperature Ratio\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);	
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Heat flux\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			#elif defined Dimensional 
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Temperature (Dimensional)\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);	
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Heat flux (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			#else
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Temperature\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);	
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Heat flux\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			#endif

		#else
			#ifdef Dimensional 
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Temperature (Dimensional)\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);	
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Heat flux (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			#else
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Temperature\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);	
				fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Heat flux\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			#endif
			
		#endif
	#endif

	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
	#ifdef ENERGY
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"EnergyCellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
	#endif
	fprintf(fp,"    </CellData>\n");
	fprintf(fp,"    <Coordinates>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+3);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Ny+3);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(zcomp+1);
	fprintf(fp,"    </Coordinates>\n");
	fprintf(fp,"  </Piece>\n");
	fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"  <AppendedData encoding=\"raw\">");
	fprintf(fp,"_");

   
	// Density (cell)
	array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<=zcomp;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx+2;i++){
				#if defined Dimensional
					val32=(float)lb.fluid[i][j][k].density*rho_dim/rho0; fwrite(&val32,sizeof(float),1,fp);
				#else
					val32=(float)lb.fluid[i][j][k].density; fwrite(&val32,sizeof(float),1,fp);
				#endif
				
			}
		}
	}

	

	// Velocity (cell)
	array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<=zcomp;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx+2;i++){
				#if defined Couette_Flow
					#if defined Velocity_Ratio
						val32=(float)lb.fluid[i][j][k].ux/u_wall; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy/u_wall; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz/u_wall; fwrite(&val32,sizeof(float),1,fp);
					#elif defined Dimensional 
						val32=(float)lb.fluid[i][j][k].ux*V_dim/u_wall; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy*V_dim/u_wall; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz*V_dim/u_wall; fwrite(&val32,sizeof(float),1,fp);
					#else 
						val32=(float)lb.fluid[i][j][k].ux; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz; fwrite(&val32,sizeof(float),1,fp);
					#endif
				#elif defined Poiseuille_Flow
					#if defined Velocity_Ratio
						val32=(float)lb.fluid[i][j][k].ux/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz/u_max; fwrite(&val32,sizeof(float),1,fp);
					#elif defined Dimensional 
						val32=(float)lb.fluid[i][j][k].ux*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
					#else 
						val32=(float)lb.fluid[i][j][k].ux; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz; fwrite(&val32,sizeof(float),1,fp);
					#endif
				#elif defined Dimensional
					val32=(float)lb.fluid[i][j][k].ux*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].uy*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].uz*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
				#else
					val32=(float)lb.fluid[i][j][k].ux/u_max; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].uy/u_max; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].uz/u_max; fwrite(&val32,sizeof(float),1,fp);
				#endif
				
			}
		}
	}

	#ifdef IBM
		// Force density (cell)
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					val32=(float)lb.fluid[i][j][k].F[0]; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].F[1]; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].F[2]; fwrite(&val32,sizeof(float),1,fp);
				}
			}
		}
	#endif
	

	#if defined MORE_DETAILS
		//Pressure (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].pressure*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].pressure; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Viscous stress
		array_size=9*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){

					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].sigma[0][0]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][1]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][2]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][0]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][1]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][2]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][0]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][1]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][2]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].sigma[0][0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][2]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][2]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][2]; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Vorticity
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					val32=(float)lb.fluid[i][j][k].vorticity[0]; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].vorticity[1]; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.fluid[i][j][k].vorticity[2]; fwrite(&val32,sizeof(float),1,fp);
				}
			}
		}

	#endif

	#ifdef ENERGY
		// Temperature (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Couette_Flow
						#if defined Temperature_Ratio
							val32=(float) (lb.fluid[i][j][k].T - Tb)/dT; fwrite(&val32,sizeof(float),1,fp);
						#elif defined Dimensional
							val32=(float)lb.fluid[i][j][k].T/Tb*Tb_dim; fwrite(&val32,sizeof(float),1,fp);
						#else
							val32=(float)lb.fluid[i][j][k].T; fwrite(&val32,sizeof(float),1,fp);
						#endif
					#elif defined Poiseuille_Flow
						#if defined Temperature_Ratio
							val32=(float) (lb.fluid[i][j][k].T)/Tb - 1; fwrite(&val32,sizeof(float),1,fp);
						#elif defined Dimensional
							val32=(float)lb.fluid[i][j][k].T/Tb*Tb_dim; fwrite(&val32,sizeof(float),1,fp);
						#else 
							val32=(float)lb.fluid[i][j][k].T; fwrite(&val32,sizeof(float),1,fp);
						#endif
					#else
						#if defined Dimensional
							val32=(float)lb.fluid[i][j][k].T/Tb*Tb_dim; fwrite(&val32,sizeof(float),1,fp);
						#else
							val32=(float)lb.fluid[i][j][k].T; fwrite(&val32,sizeof(float),1,fp);
						#endif
					#endif
					
					
				}
			}
		}

		// Heat Flux (cell)
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].q[0]*rho_dim*cv_dim*c_dim*Tb_dim; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].q[1]*rho_dim*cv_dim*c_dim*Tb_dim; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].q[2]*rho_dim*cv_dim*c_dim*Tb_dim; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].q[0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].q[1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].q[2]; fwrite(&val32,sizeof(float),1,fp);
					#endif
				}
			}
		}
	#endif
	
    // CellType (cell)
	array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<=zcomp;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx+2;i++){
				val32=(int)lb.fluid[i][j][k].type; fwrite(&val32,sizeof(int),1,fp);
			}
		}
	}

	#ifdef ENERGY
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					val32=(int)lb.fluid[i][j][k].type_Energy; fwrite(&val32,sizeof(int),1,fp);
				}
			}
		}
	#endif
	// Coordinates (vertices)
	array_size=1*4*(Nx+3);
	fwrite(&array_size,sizeof(int),1,fp);
	for(i=0;i<Nx+3;i++){ val32=(float)(i*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Ny+3);
	fwrite(&array_size,sizeof(int),1,fp);
	#if defined Poiseuille_Flow
		double H = Ny+2;
		for(j=0;j<Ny+3;j++){ val32=(float)(j*dy)/H*2-1; fwrite(&val32,sizeof(float),1,fp); }
	#else
		for(j=0;j<Ny+3;j++){ val32=(float)(j*dy); fwrite(&val32,sizeof(float),1,fp); }
	#endif
	

	array_size=1*4*(zcomp+1);
	fwrite(&array_size,sizeof(int),1,fp);
	#if defined Poiseuille_Flow
		for(k=0;k<zcomp+1;k++){ val32=(float)(k*dz)/H*2-1; fwrite(&val32,sizeof(float),1,fp); }
	#else
		for(k=0;k<zcomp+1;k++){ val32=(float)(k*dz); fwrite(&val32,sizeof(float),1,fp); }
	#endif
	

	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

#ifdef MORE_DETAILS
	void OutputAVGVTK(int &nout, LBM &lb) {
		int		i,j,k;
		char	filename[128];
		FILE	*fp;
		unsigned int array_size;
		unsigned long int offset=0;
		// short  num16; // Int16 2byte
		float  val32; // Float32 4byte

		#ifdef Cylinder_Tandem_Triangle
			sprintf(filename,"./output Tandem Triangle/LpD = %.1f Re = %.0f/avgfield.vtr",LpD, Re);
		#else

		#endif
		

		fp=fopen(filename,"wb");

		#if defined D2Q9
			int zcomp = Nz;	
		#else 
			int zcomp = Nz;
		#endif
		fprintf(fp,"<?xml version=\"1.0\"?>\n");
		fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
		fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",Nx+2,Ny+2,zcomp);
		fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Nx+2,Ny+2,zcomp);
		fprintf(fp,"    <PointData>\n");
		fprintf(fp,"    </PointData>\n");
		fprintf(fp,"    <CellData>\n");
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);

		#ifndef Dimensional
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Pressure\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg pressure\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms pressure\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Viscous Stress\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg vorticity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms vorticity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg reynolds Stress\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms reynolds Stress\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
		#else 
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg velocity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms velocity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Pressure (Dimensional)\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg pressure (Dimensional)\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms pressure (Dimensional)\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Viscous Stress (Dimensional)\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Vorticity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg vorticity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms vorticity (Dimensional)\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(zcomp);	
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Avg reynolds Stress (Dimensional)\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
			fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Rms reynolds Stress (Dimensional)\" NumberOfComponents=\"9\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+9*4*(Nx+2)*(Ny+2)*(zcomp);
		#endif
		
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(zcomp);
		fprintf(fp,"    </CellData>\n");
		fprintf(fp,"    <Coordinates>\n");
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+3);
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Ny+3);
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(zcomp+1);
		fprintf(fp,"    </Coordinates>\n");
		fprintf(fp,"  </Piece>\n");
		fprintf(fp,"  </RectilinearGrid>\n");
		fprintf(fp,"  <AppendedData encoding=\"raw\">");
		fprintf(fp,"_");

		// Density (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].density*rho_dim/rho0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].density; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		// Velocity (cell)
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].ux*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].ux/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uy/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].uz/u_max; fwrite(&val32,sizeof(float),1,fp);
					#endif	
				}
			}
		}

		// Average velocity (cell)
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].avg_velocity[0]*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_velocity[1]*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_velocity[2]*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].avg_velocity[0]/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_velocity[1]/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_velocity[2]/u_max; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		// Rms velocity (cell)
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].rms_velocity[0]*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_velocity[1]*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_velocity[2]*V_dim/u_max; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].rms_velocity[0]/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_velocity[1]/u_max; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_velocity[2]/u_max; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Pressure (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].pressure*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].pressure; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Average Pressure (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].avg_pressure*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].avg_pressure; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Rms Pressure (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].rms_pressure*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].rms_pressure; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Viscous stress
		array_size=9*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){

					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].sigma[0][0]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][1]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][2]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][0]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][1]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][2]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][0]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][1]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][2]*pressure_dim/p0; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].sigma[0][0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[0][2]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[1][2]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][0]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][1]; fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].sigma[2][2]; fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Vorticity
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#ifdef Dimensional
						val32=(float)lb.fluid[i][j][k].vorticity[0]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].vorticity[1]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].vorticity[2]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].vorticity[0]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].vorticity[1]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].vorticity[2]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Average Vorticity
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#ifdef Dimensional
						val32=(float)lb.fluid[i][j][k].avg_vorticity[0]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_vorticity[1]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_vorticity[2]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].avg_vorticity[0]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_vorticity[1]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_vorticity[2]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
					#endif
				}
			}
		}

		//Rms Vorticity
		array_size=3*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					#ifdef Dimensional
						val32=(float)lb.fluid[i][j][k].rms_vorticity[0]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_vorticity[1]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_vorticity[2]/(u_max/D)*(V_dim/D_dim); fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].rms_vorticity[0]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_vorticity[1]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_vorticity[2]/(u_max/D); fwrite(&val32,sizeof(float),1,fp);
					#endif
				}
			}
		}

		//Average Reynolds stress
		array_size=9*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){

					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[0][0]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[0][1]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[0][2]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[1][0]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[1][1]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[1][2]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[2][0]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[2][1]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[2][2]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[0][0]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[0][1]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[0][2]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[1][0]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[1][1]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[1][2]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[2][0]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[2][1]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].avg_ReynoldStress[2][2]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		//Rms Reynolds stress
		array_size=9*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){

					#if defined Dimensional
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[0][0]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[0][1]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[0][2]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[1][0]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[1][1]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[1][2]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[2][0]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[2][1]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[2][2]/(pow(u_max,2))*pow(V_dim,2); fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[0][0]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[0][1]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[0][2]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[1][0]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[1][1]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[1][2]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[2][0]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[2][1]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
						val32=(float)lb.fluid[i][j][k].rms_ReynoldStress[2][2]/(pow(u_max,2)); fwrite(&val32,sizeof(float),1,fp);
					#endif
					
				}
			}
		}

		// CellType (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(zcomp);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=zcomp;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					val32=(int)lb.fluid[i][j][k].type; fwrite(&val32,sizeof(int),1,fp);
				}
			}
		}

		// Coordinates (vertices)
		array_size=1*4*(Nx+3);
		fwrite(&array_size,sizeof(int),1,fp);
		for(i=0;i<Nx+3;i++){ val32=(float)(i*dx); fwrite(&val32,sizeof(float),1,fp); }

		array_size=1*4*(Ny+3);
		fwrite(&array_size,sizeof(int),1,fp);
		for(j=0;j<Ny+3;j++){ val32=(float)(j*dy); fwrite(&val32,sizeof(float),1,fp); }
		
		

		array_size=1*4*(zcomp+1);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=0;k<zcomp+1;k++){ val32=(float)(k*dz); fwrite(&val32,sizeof(float),1,fp); }
		

		fprintf(fp,"  </AppendedData>\n");
		fprintf(fp,"</VTKFile>\n");

		fclose(fp);
	}

#endif

#if defined IBM
	void OutputVTK_Marker(int& time, vector<MARKER>& marker) {
		//Create filename 
		stringstream output_filename;
		output_filename << "wall_t "<< int(time/Nt) << ".vtk";
		ofstream output_file;
		output_file.open(output_filename.str().c_str());
		
		// Write VTK header
	    output_file << "# vtk DataFile Version 3.0\n";
	    output_file << "wall_state\n";
	    output_file << "ASCII\n";
	    output_file << "DATASET POLYDATA\n";
	    
	    //Write node positions
	    output_file<< "POINTS " << number_nodes << " float\n";
	    

	    for (int n = 0; n<number_nodes; n++) {
	    	output_file << marker[n].x << " " << marker[n].y<< " " <<marker[n].z << " \n";
		}

		// Write lines between neighboring nodes
		/*
		output_file << "LINES " << number_nodes - 2 << " " << 3 * (number_nodes - 2) << "\n";
		for (int i = 0; i<number_nodes; i++) {
			output_file << "2 " << i%number_nodes << " " << (i + 1) << "\n";
		}
		*/

		output_file << "LINES " << number_nodes << " " << 3 * (number_nodes) << "\n";
      for(int n = 0; n < number_nodes; ++n) {
        output_file << "2 " << n%(number_nodes) << " " << (n + 1)%(number_nodes) << "\n";
      }
      
	
		// Write vertices
		output_file << "VERTICES 1 " << number_nodes + 1 << "\n";
		output_file << number_nodes << " ";

		for(int n = 0; n < number_nodes; n++) {
		output_file << n << " ";
		}

		output_file.close();
	}
	
#endif

//---------------------------------------------------FORCE CALCULATION--------------------------------------------

#if defined MOMENTUM_EXCHANGE_ALGORITHM
	//Momentum exchange Algorithm : to calculate CL and CD
	double MEA_CL(LBM &lb, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) {
		double Py = 0; //Total momentum exchange in y direction
		//#pragma omp parallel for
		for(int i=xmin; i<=xmax; ++i) {
			for(int j=ymin; j<=ymax; ++j) {
				for(int k = zmin; k<=zmax; ++k) {
					
					if(lb.fluid[i][j][k].type==TYPE_F) {
						int i_nb, j_nb, k_nb;
						for (int l=0; l < npop; ++l)
						{
							i_nb = i - cx[l];
							j_nb = j - cy[l];
							k_nb = k - cz[l]; 
							//---- Solid Boundary Condition ----------------------
							
							if(lb.fluid[i_nb][j_nb][k_nb].type==TYPE_S)
							{
								//Py += lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] + lb.fluid[i][j][k].f[l]*cy[l];
								//Py += lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cy[l]; //usually used
								//Py += (lb.fluid[i][j][k].f_star[opposite[l]] + lb.fluid[i][j][k].f[l])*cy[opposite[l]];
								Py += (lb.fluid[i][j][k].f_star[l] + lb.fluid[i][j][k].f[opposite[l]])*cy[l];

							}
						}
					}
					
					
				}
			}
		} 
		
		Py *= pow(dx,3)/dt;
		#if defined Cylinder_Basic or defined Cylinder_Tandem_Circle
			double CL = Py/(0.5*rho0*pow(u_max,2)*D); 
		#elif defined Cylinder_Tandem_Triangle
			double CL = Py/(0.5*rho0*pow(u_max,2)*2*D/sqrt(3)); 
		#elif defined Sphere_BB
			double CL = Py/(0.5*rho0*pow(u_max,2)*D*D*PI/4); 
		#endif
		return CL;
	}

	double MEA_CD(LBM &lb, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) {
		double Px = 0; //Total momentum exchange in x direction
		//#pragma omp parallel for 
		for(int i=xmin; i<=xmax; ++i) {
			for(int j=ymin; j<=ymax; ++j) {
				for(int k = zmin; k<=zmax; ++k) {
					if(lb.fluid[i][j][k].type==TYPE_F) {
						int i_nb, j_nb, k_nb;
						for (int l=0; l < npop; ++l)
						{
							i_nb = i - cx[l];
							j_nb = j - cy[l];
							k_nb = k - cz[l]; 
							//---- Solid Boundary Condition ----------------------
							
							if(lb.fluid[i_nb][j_nb][k_nb].type==TYPE_S)
							{
								//Px += lb.fluid[i][j][k].f[opposite[l]]*cx[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cx[l]; //usually used
								//Px += (lb.fluid[i][j][k].f_star[opposite[l]] + lb.fluid[i][j][k].f[l])*cx[opposite[l]];
								Px += -(lb.fluid[i][j][k].f_star[l] + lb.fluid[i][j][k].f[opposite[l]])*cx[l];
							}
						}
					}
				}
			}
		} 
		Px *= pow(dx,3)/dt;
		#if defined Cylinder_Basic or defined Cylinder_Tandem_Circle
			double CD = Px/(0.5*rho0*pow(u_max,2)*D); 
		#elif defined Cylinder_Tandem_Triangle
			double CD = Px/(0.5*rho0*pow(u_max,2)*2*D/sqrt(3)); 
		#elif defined Sphere_BB
			double CD = Px/(0.5*rho0*pow(u_max,2)*D*D*PI/4); 
		#endif
		return CD;
	}
#endif

//----------------------------------------------------------------------------------------------------------------

//--------------------------------------------------OUTPUTING IN CSV----------------------------------------------
void OutputCSV(vector<double>& error, string name) { 
	ofstream ofs;
	ofs.open(name);
	ofs << "time,error\n";
	for (int i = 0; i<error.size(); i++) {
		#if defined Couette_Flow or defined Poiseuille_Flow or defined Rayleigh_Benard
			ofs << (i+1)*Nt * dt << "," << error[i]<<"\n"; //error here
		#elif defined Cylinder_IBM or defined Cylinder_IBM_Moving_Arbie
			ofs << (i+1) * dt *u_max/(D) << "," << error[i]<<"\n";
		#elif defined Cylinder_Inline_Oscillating
			ofs << (i+1) * dt *frequency << "," << error[i]<<"\n";
		#else
			ofs << (i+1)*Nt * dt *u_max/(D) << "," << error[i]<<"\n";
		#endif
		
	}
	ofs.close();
}

//----------------------------------------------------------------------------------------------------------------



//-----------------------------------------ANALYTICAL SOLUTION ---------------------------------
#if defined Couette_Flow
	//Velocity analytical solution 
    double Couette_Analytic_Velocity(double t, double y) {
        double h = (Ny+1)*dy;
        double sum = 0;
        //Transient term 
        for (int n = 1; n<=100; n++) {
            sum += pow(-1,n)/n*exp(-nu*t*pow(n*PI/h,2))*sin(n*PI*y/h);
        }
        //Steady term 
        double u = u_wall*(y/h + 2*sum/PI);

        return u;
    }

    double Velocity_Error_Calculation(LBM& lb, int step) {
        //Analytical result
        double u_analytic[Ny+1];
        double t = step*dt,y;

        for(int j = 0; j<Ny+1; j++) {
            y = j*dy;
            u_analytic[j] = Couette_Analytic_Velocity(t,y);
        }

        //error calculation
        double sum_error_num = 0;
        double sum_error_denom = 0;
        double error;
        double u;

        for(int j = 0; j<Ny+1; j++) {
            u = lb.fluid[Nx][j][Nz].ux;
			sum_error_num += pow( (u - u_analytic[j]) ,2);
			sum_error_denom += pow( (u_analytic[j]) ,2);
        }

        if(sum_error_denom == 0) error = sqrt(sum_error_num);
		else error = sqrt(sum_error_num/sum_error_denom);

        return error;
    }

    #if defined MORE_DETAILS
        //Viscous stress analytical solution
        double Couette_Analytic_Viscous(double t, double y) {
            double h = (Ny+1)*dy;
            double sum = 0;
            //Transient term 
            for (int n = 1; n<=100; n++) {
                sum += pow(-1,n)*exp(-nu*t*pow(n*PI/h,2))*cos(n*PI*y/h);
            }
            //Steady term 
            double rho = 1.0;
            double mu = nu*rho;
            double s = mu*u_wall/h*(1+2*sum);

            return s;
        }

        double Viscous_Error_Calculation(LBM& lb, int step) {
            //Analytical result
            double s_analytic[Ny+1];
            double t,y;

            for(int j = 1; j<=Ny; j++) {
                t = step*dt;
                y = j*dy;
                s_analytic[j] = Couette_Analytic_Viscous(t,y);
            }

            //error calculation
            double sum_error_num = 0;
            double sum_error_denom = 0;
            double error;
            double s;

            for(int j = 1; j<Ny+1; j++) {
                s = lb.fluid[Nx][j][Nz].sigma[0][1];
                sum_error_num += pow( (s - s_analytic[j]) ,2);
                sum_error_denom += pow( (s_analytic[j]) ,2);
            }

            if(sum_error_denom == 0) error = sqrt(sum_error_num);
            else error = sqrt(sum_error_num/sum_error_denom);

            return error;
        }
    #endif

	#ifdef ENERGY
		double Couette_Analytic_Temperature(double t, double y) {
			double h = (Ny+1)*dy;
			double T = Tb + dT*y/h*(1+0.5*Br*(1-y/h));
			return T;
		}

		double Temperature_Error_Calculation(LBM& lb, int step) {
			//Analytical value calculation 
			double T_analytic[Ny+1];
			double ratio_analytic[Ny+1];
			for (int j = 0; j<Ny+1; j++) {
				double y = j*dy;
				double t = step*dt;
				T_analytic[j] = Couette_Analytic_Temperature(t, y);
				ratio_analytic[j] = (T_analytic[j] - Tb)/dT;
			}

			//Error calculation 
			double sum_error_denom = 0;
			double sum_error_num = 0;
			double error;
			double T; 	//numerical value of temperature
			double ratio;

			for(int j = 0; j<Ny+1; j++) {
				T = lb.fluid[Nx][j][Nz].T;
				ratio = (T-Tb)/dT;
				sum_error_num += pow( (ratio - ratio_analytic[j]) ,2);
				sum_error_denom += pow( (ratio_analytic[j]) ,2);
			}
			if(sum_error_denom == 0) error = sqrt(sum_error_num);
			else error = sqrt(sum_error_num/sum_error_denom);

			return error;
		}

		double Couette_Analytic_HeatFlux(double t, double y) {
			double h = (Ny+1)*dy;
			double qy = -k*(dT*(1+0.5*Br)/h - (Br*dT*y)/pow(h,2));
			return qy;
		}

		double HeatFlux_Error_Calculation(LBM& lb, int step) {
			//Analytical value calculation 
			double q_analytic[Ny+1];
			
			for (int j = 0; j<Ny+1; j++) {
				double y = j*dy;
				double t = step*dt;
				q_analytic[j] = Couette_Analytic_HeatFlux(t, y);
			}

			//Error calculation 
			double sum_error_denom = 0;
			double sum_error_num = 0;
			double error;
			double q; 	//numerical value of heat flux

			for(int j = 2; j<Ny; j++) {
				q = lb.fluid[Nx][j][Nz].q[1];
				sum_error_num += pow( (q - q_analytic[j]) ,2);
				sum_error_denom += pow( (q_analytic[j]) ,2);
			}
			if(sum_error_denom == 0) error = sqrt(sum_error_num);
			else error = sqrt(sum_error_num/sum_error_denom);

			return error;
		}
	#endif

#elif defined Poiseuille_Flow
	double Poiseuille_Analytic_Velocity(double t, double y) {
        double h = (Ny+1)*dy;
        double sum = 0;

        //Transient term 
        for (int k = 1; k<=100; k++) {
            double n = 2*k-1;
            sum += (1-pow(-1,n))*exp(-nu*t*pow(n*PI/h,2))*sin(n*PI*y/h)/pow(n,3);
        }

        //Steady and Transient term 
        double u = 1/nu*dp_dx*(0.5*y*(y-h) + 2*pow(h,2)/pow(PI,3)*sum);
		// double u = 4*u_max*y*(h-y)/pow(h,2);
        return u;
    }

    double Velocity_Error_Calculation(LBM& lb, int step) {
        //Analytical result
        double u_analytic[Ny+2];
        double t,y;

        for(int j = 0; j<=Ny+1; j++) {
            t = step*dt;
            y = j*dy;
            u_analytic[j] = Poiseuille_Analytic_Velocity(t,y);
        }

        //error calculation
        double sum_error_num = 0;
        double sum_error_denom = 0;
        double error;
        double u;

        for(int j = 0; j<Ny+2; j++) {
            u = lb.fluid[Nx][j][Nz].ux;
			sum_error_num += pow( (u - u_analytic[j]) ,2);
			sum_error_denom += pow( (u_analytic[j]) ,2);

            //cout<<"Analytic = "<<u_analytic[j]<<" "<<"Numeric = "<<u<<endl;
        }

        if(sum_error_denom == 0) error = sqrt(sum_error_num);
		else error = sqrt(sum_error_num/sum_error_denom);

        return error;
    }
    
    #if defined MORE_DETAILS
        double Poiseuille_Analytic_Viscous(double t, double y) {
            double h = (Ny+1)*dy;
            double sum = 0;

            //Transient term 
            for (int n = 1; n<=100; n++) {
                sum += (1-pow(-1,n))*exp(-nu*t*pow(n*PI/h,2))*cos(n*PI*y/h)/pow(n,2);
            }

            double s = dp_dx*(y - h/2 + 2*h*sum/pow(PI,2));
            return s;
        }

       double Viscous_Error_Calculation(LBM& lb, int step) {
            //Analytical result
            double s_analytic[Ny+1];
            double t,y;

            for(int j = 1; j<Ny+1; j++) {
                t = step*dt;
                y = j*dy;
                s_analytic[j] = Poiseuille_Analytic_Viscous(t,y);
            }

            //error calculation
            double sum_error_num = 0;
            double sum_error_denom = 0;
            double error;
            double s;

            for(int j = 1; j<Ny+1; j++) {
                s = lb.fluid[Nx][j][Nz].sigma[0][1];
                sum_error_num += pow( (s - s_analytic[j]) ,2);
                sum_error_denom += pow( (s_analytic[j]) ,2);
            }

            if(sum_error_denom == 0) error = sqrt(sum_error_num);
            else error = sqrt(sum_error_num/sum_error_denom);

            return error;
        }

    #endif

    #if defined ENERGY
		//Temperature
        double Poiseuille_Analytic_Temperature(double t, double y) {
            double H = (Ny+1)*dy;
            double T = Tb + Pr*pow(u_max,2)*(1-pow((2*y/H-1),4))/(3*cp) + (Tt - Tb)*y/H;    //steady state solution
            return T;
        }

        double Temperature_Error_Calculation(LBM& lb, int step) {
            //Analytical result
            double T_analytic[Ny+2];
            double t,y;

            for(int j = 0; j<=Ny+1; j++) {
                t = step*dt;
                y = j*dy;
                T_analytic[j] = (Poiseuille_Analytic_Temperature(t,y)-1);
            }

            //error calculation
            double sum_error_num = 0;
            double sum_error_denom = 0;
            double error;
            double T;

            for(int j = 1; j<Ny+1; j++) {
                T = (lb.fluid[Nx][j][Nz].T-1);
                sum_error_num += pow( (T - T_analytic[j]) ,2);
                sum_error_denom += pow( (T_analytic[j]) ,2);
            }

            if(sum_error_denom == 0) error = sqrt(sum_error_num);
            else error = sqrt(sum_error_num/sum_error_denom);

            return error;
        }

		//Heat Flux
		double Poiseuille_Analytic_HeatFlux(double t, double y) {
            double H = (Ny+1)*dy;
            double q = -k*(Pr*pow(u_max,2)/(3*cp))*(-4*pow(2*y/H-1,3)*2/H);    //steady state solution
            return q;
        }

        double HeatFlux_Error_Calculation(LBM& lb, int step) {
            //Analytical result
            double qy_analytic[Ny+2];
            double t,y;

            for(int j = 0; j<=Ny+1; j++) {
                t = step*dt;
                y = j*dy;
                qy_analytic[j] = Poiseuille_Analytic_HeatFlux(t,y);
            }

            //error calculation
            double sum_error_num = 0;
            double sum_error_denom = 0;
            double error;
            double qy;

            for(int j = 2; j<Ny; j++) {
                qy = lb.fluid[Nx][j][Nz].q[1];
                sum_error_num += pow( (qy - qy_analytic[j]) ,2);
                sum_error_denom += pow( (qy_analytic[j]) ,2);
            }

            if(sum_error_denom == 0) error = sqrt(sum_error_num);
            else error = sqrt(sum_error_num/sum_error_denom);

            return error;
        }

    #endif
#elif defined Rayleigh_Benard
	double Find_Maximum_uy(LBM& lb) {
		double maximum = -99;

		#pragma omp parallel for 
		for(int i = 0; i<=Nx+1; i++) {
			for(int j = 0; j<=Ny+1; j++) {
				for(int k = 0; k<=Nz+1; k++) {
					if(maximum < abs(lb.fluid[i][j][k].uy)) {
						maximum = abs(lb.fluid[i][j][k].uy);
					}
				}
			}
		}
		return maximum;
	}

	double Find_Nusselt_number(LBM& lb) {
		#ifdef ENERGY
			//Finding average of uy*T : 
			double sum_uyT = 0;
			int number = 0; 

			#pragma omp parallel for 	
			for(int i = 0; i<=Nx+1; i++) {
				for(int j = 0; j<=Ny+1; j++) {
					for(int k = 0; k<=Nz+1; k++) {
						if(lb.fluid[i][j][k].type == TYPE_F or lb.fluid[i][j][k].type == TYPE_OUT or lb.fluid[i][j][k].type == TYPE_IN ) {
							sum_uyT += lb.fluid[i][j][k].uy*lb.fluid[i][j][k].T;
							number++;
						}
					}
				}
			}

			double avg_uyT = sum_uyT/number;

			double Nusselt = 1 + avg_uyT/(alpha*dT/(Ny+1));
			return Nusselt;
		#endif
	}
#endif
//----------------------------------------------------------------------------------------------------------------


//-----------------------------------------POST PROCESSING ---------------------------------
#if defined Couette_Flow or defined Poiseuille_Flow
	void Postprocessing(LBM& lb, int step) {
		//Defining variables :
		//Numerical and Analytical Value : 
		vector <double> Velocity_numeric; 
		vector <double> Velocity_analytic; 

		#if defined MORE_DETAILS
			vector <double> Viscous_numeric; 
			vector <double> Viscous_analytic; 
		#endif

		#if defined ENERGY
			vector <double> Temperature_numeric; 
			vector <double> Temperature_analytic; 
			vector <double> HeatFlux_numeric; 
			vector <double> HeatFlux_analytic; 
		#endif

		for(int j = 1; j<Ny+1; j++) {
			//Velocity
			Velocity_numeric.push_back(lb.fluid[Nx-1][j][Nz].ux/u_max);
			OutputCSV(Velocity_numeric, "Velocity Numeric.csv");

			#if defined Poiseuille_Flow
				Velocity_analytic.push_back(Poiseuille_Analytic_Velocity(step*dt, j*dy)/u_max);
			#elif defined Couette_Flow
				Velocity_analytic.push_back(Couette_Analytic_Velocity(step*dt, j*dy)/u_max);
			#endif			

			OutputCSV(Velocity_analytic, "Velocity Analytic.csv");

			//Viscous stress
			#ifdef MORE_DETAILS
				#if defined Couette_Flow
					Viscous_numeric.push_back(lb.fluid[Nx-1][j][Nz].sigma[0][1]/(nu*rho0*u_max/(Ny+1)));
				#elif defined Poiseuille_Flow
					Viscous_numeric.push_back(lb.fluid[Nx-1][j][Nz].sigma[0][1]/(0.5*(Ny+1)*dp_dx));
				#endif
				
				OutputCSV(Viscous_numeric, "Viscous stress Numeric.csv");
			
				#if defined Poiseuille_Flow
					Viscous_analytic.push_back(Poiseuille_Analytic_Viscous(step*dt, j*dy)/(0.5*(Ny+1)*dp_dx));
				#elif defined Couette_Flow
					Viscous_analytic.push_back(Couette_Analytic_Viscous(step*dt, j*dy)/(nu*rho0*u_max/(Ny+1)));
				#endif
				OutputCSV(Viscous_analytic, "Viscous stress Analytic.csv");
			#endif

			#ifdef ENERGY
				//Temperature and heat flux 
				#if defined Couette_Flow
					Temperature_numeric.push_back((lb.fluid[Nx-1][j][Nz].T - Tb)/dT);
					HeatFlux_numeric.push_back(lb.fluid[Nx-1][j][Nz].q[1]);
				#elif defined Poiseuille_Flow
					Temperature_numeric.push_back((lb.fluid[Nx-1][j][Nz].T/Tb-1)*pow(10,6));
					HeatFlux_numeric.push_back( lb.fluid[Nx-1][j][Nz].q[1]/((Tt - Tb)/(Ny+1) - 8*Pr*pow(u_max,2)/(3*cp*(Ny+1))) );
				#endif
				

				OutputCSV(Temperature_numeric, "Temperature Numeric.csv");
				OutputCSV(HeatFlux_numeric, "Heat Flux Numeric.csv");

				#if defined Poiseuille_Flow
					Temperature_analytic.push_back((Poiseuille_Analytic_Temperature(step*dt, j*dy)/Tb-1)*pow(10,6));
					HeatFlux_analytic.push_back(Poiseuille_Analytic_HeatFlux(step*dt, j*dy)/((Tt - Tb)/(Ny+1) - 8*Pr*pow(u_max,2)/(3*cp*(Ny+1))));
				#elif defined Couette_Flow
					Temperature_analytic.push_back((Couette_Analytic_Temperature(step*dt, j*dy) - Tb)/dT);
					HeatFlux_analytic.push_back(Couette_Analytic_HeatFlux(step*dt, j*dy));
				#endif

				OutputCSV(Temperature_analytic, "Temperature Analytic.csv");
				OutputCSV(HeatFlux_analytic, "Heat Flux Analytic.csv");
			#endif
		}
	}
#endif
