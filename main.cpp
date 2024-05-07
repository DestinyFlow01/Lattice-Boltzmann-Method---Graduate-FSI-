#include<iostream>
#include<ctime>
#include<vector>
#include<string>
#include "LBM.h"
using namespace std;



int main() {
	// print_Logo();
	print_Info();
	//ofstream ofs;
	
	clock_t c_start = clock();

	//Defining variable for simulation : 
	#if defined Cylinder_Basic
		LBM lb = main_setup_Cylinder();
	#elif defined Poiseuille_Flow
		LBM lb = main_setup_Poiseuille();
	#elif defined Couette_Flow
		LBM lb = main_setup_Couette();
	#elif defined Rayleigh_Benard
		LBM lb = main_setup_RayleighBenard1();
	#elif defined Cylinder_Tandem_Circle
		LBM lb = main_setup_Tandem_Circle();
	#elif defined Cylinder_Tandem_Triangle
		LBM lb = main_setup_Tandem_Triangle();
	#elif defined Sphere_BB
		LBM lb = main_setup_Sphere();
	#elif defined Cylinder_IBM
		LBM lb = main_setup_Cylinder_IBM();
		
		//Markers : 
		vector <MARKER> marker(number_nodes);
		Generate_Cylinder_Marker(marker);
	#elif defined Cylinder_Inline_Oscillating
		LBM lb = main_setup_Cylinder_InlineOscillating();
		
		//Markers : 
		vector <MARKER> marker(number_nodes);
		Generate_Cylinder_Marker(marker);
	#elif defined Cylinder_Transverse_Oscillation
		LBM lb = main_setup_Cylinder_TransverseOscillating();
		
		//Markers : 
		vector <MARKER> marker(number_nodes);
		Generate_Cylinder_Marker(marker);
	#elif defined Cylinder_IBM_Moving_Arbie
		LBM lb = main_setup_Cylinder_IBM_Moving();
		
		//Markers : 
		vector <MARKER> marker(number_nodes);
		Generate_Cylinder_Marker(marker);
	#endif
	

	lb.Init();
	
	//Preliminary Result : 
	int step = 0, n = 0;
	OutputVTK(step,lb); //Fluid

	#ifdef IBM
		OutputVTK_Marker(step, marker);
	#endif
	
	
	int _NUMBER = floor(T_OUT/Nt);

	#if defined Cylinder_Basic
		vector <double> CL; //Lift Coefficient
		vector <double> CD; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now
	#elif defined Poiseuille_Flow
		//Error calculation 
		vector <double> error_Velocity; //Velocity error

		#if defined MORE_DETAILS
			vector <double> error_Viscous; //Viscous error
		#endif

		#if defined ENERGY
			vector<double> error_Temperature; //Temperature error
			vector<double> error_HeatFlux; //Heat Flux error
		#endif
		int total_data = 0; //Total data for the simulation right now
	#elif defined Couette_Flow
		vector <double> error_Velocity; //Velocity error

		#if defined MORE_DETAILS
			vector <double> error_Viscous; //Viscous error
		#endif

		#if defined ENERGY
			vector<double> error_Temperature; //Temperature error
			vector<double> error_HeatFlux; //Heat Flux error
		#endif
		int total_data = 0; //Total data for the simulation right now
	#elif defined Rayleigh_Benard
		#ifdef ENERGY
			vector<double> uy_max;
			vector<double> Nusselt_number;
		#endif
	#elif defined Cylinder_Tandem_Circle
		vector <double> CL[N]; //Lift coefficient
		vector <double> CD[N]; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now
	#elif defined Cylinder_Tandem_Triangle
		vector <double> CL[N]; //Lift coefficient
		vector <double> CD[N]; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now
	#elif defined Sphere_BB
		vector <double> CL; //Lift Coefficient
		vector <double> CD; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now
	#elif defined Cylinder_IBM
		vector <double> CL; //Lift coefficient
		vector <double> CD; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now

		vector<double> ds;
		for (int i = 0; i<number_nodes; i++) {
			int j = (i+1)%(number_nodes);
			double deltax = marker[i].x - marker[j].x;
			double deltay = marker[i].y - marker[j].y;
			double deltaz = marker[i].z - marker[j].z;

			ds.push_back(sqrt(pow(deltax,2) + pow(deltay,2) + pow(deltaz,2)));
		}

		//center of mass velocity :
		double V_CM[3] = {0,0,0};

	#elif defined Cylinder_Inline_Oscillating
		vector <double> CL; //Lift coefficient
		vector <double> CD; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now

		vector<double> ds;
		for (int i = 0; i<number_nodes; i++) {
			int j = (i+1)%(number_nodes);
			double deltax = marker[i].x - marker[j].x;
			double deltay = marker[i].y - marker[j].y;
			double deltaz = marker[i].z - marker[j].z;

			ds.push_back(sqrt(pow(deltax,2) + pow(deltay,2) + pow(deltaz,2)));
		}

		double V_CM[3] = {0,0,0};
	#elif defined Cylinder_Transverse_Oscillation
		vector <double> CL; //Lift coefficient
		vector <double> CD; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now

		vector<double> ds;
		for (int i = 0; i<number_nodes; i++) {
			int j = (i+1)%(number_nodes);
			double deltax = marker[i].x - marker[j].x;
			double deltay = marker[i].y - marker[j].y;
			double deltaz = marker[i].z - marker[j].z;

			ds.push_back(sqrt(pow(deltax,2) + pow(deltay,2) + pow(deltaz,2)));
		}

		double V_CM[3] = {0,0,0};

	#elif defined Cylinder_IBM_Moving_Arbie
		vector <double> CL; //Lift coefficient
		vector <double> CD; //Drag Coefficient
		int total_data = 0; //Total data for the simulation right now

		vector<double> ds;
		for (int i = 0; i<number_nodes; i++) {
			int j = (i+1)%(number_nodes);
			double deltax = marker[i].x - marker[j].x;
			double deltay = marker[i].y - marker[j].y;
			double deltaz = marker[i].z - marker[j].z;

			ds.push_back(sqrt(pow(deltax,2) + pow(deltay,2) + pow(deltaz,2)));
		}

		//center of mass velocity :
		double V_CM[3] = {0,0,0};
	#endif
	
	

	
	//Iterations 
	string name;
	int numstep = T_OUT/dt;
	cout<<"Total number of step = "<<numstep<<endl;
	for(step = 1; step<=numstep; step++) {
		
		#if defined IBM
			//IBM step : 
			IBLBM(lb, marker, ds, CL, CD, V_CM, step*dt);
			
			//Collision
			lb.Collision();
			
			//Streaming and BC (Since the initialization has been done)
			lb.Streaming();

			//Output partial solution 
			name = "Cl at Re = "+to_string(Re) + ".csv";
			OutputCSV(CL, name); //Output CL
			name = "Cd at Re = "+to_string(Re) + ".csv";
			OutputCSV(CD, name); //Output CD
		#else
			//Collision
			lb.Collision();
			
			//Streaming and BC (Since the initialization has been done)
			lb.Streaming();	
			
			//Macroscopic Property
			lb.MacroProp();
		#endif
		
		
		
		
		//Outputing file : 
		if(step%Nt == 0) {
			n++;
			OutputVTK(n,lb); // For fluid	

			#ifdef MORE_DETAILS
				OutputAVGVTK(n,lb); // For fluid
			#endif
			
			#ifdef IBM
				OutputVTK_Marker(step, marker);
			#endif
			
			cout<<"\nStep = "<<step<<endl;
			
			//Additional steps :
			#if defined Cylinder_Basic
				total_data++;
				
				#if defined MOMENTUM_EXCHANGE_ALGORITHM
					CL.push_back(MEA_CL(lb, center_x - D, center_x + D, center_y - D, center_y + D,1,1));
					cout<<"Lift coefficient at time step "<<step<<" = "<<CL[total_data-1]<<endl;

					CD.push_back(MEA_CD(lb, center_x - D, center_x + D, center_y - D, center_y + D,1,1));
					cout<<"Drag coefficient at time step "<<step<<" = "<<CD[total_data-1]<<endl;
				#endif
				

				//Output partial solution 
				name = "Value of Cl at Re = "+to_string(Re) + ".csv";
				OutputCSV(CL, name); //Output CL
				name = "Value of Cd at Re = "+to_string(Re) + ".csv";
				OutputCSV(CD, name); //Output CD
			#elif defined Poiseuille_Flow
				total_data++; Postprocessing(lb, step);

				//Error calculation and output partial solution : 
				error_Velocity.push_back(Velocity_Error_Calculation(lb, step));
				cout<<"Error for velocity "<<" = "<<error_Velocity.back()<<endl;
				OutputCSV(error_Velocity, "Velocity Error.csv");

				#if defined MORE_DETAILS
					error_Viscous.push_back(Viscous_Error_Calculation(lb, step));
					cout<<"Error for viscous stress = "<<error_Viscous.back()<<endl;
					OutputCSV(error_Viscous, "Viscous stress Error.csv"); //Output Viscous stress Error
				#endif

				#ifdef ENERGY
					error_Temperature.push_back(Temperature_Error_Calculation(lb, step));
					cout<<"Error for temperature = "<<error_Temperature.back()<<endl;
					OutputCSV(error_Temperature, "Temperature Error.csv"); //Output Temperature Error

					error_HeatFlux.push_back(HeatFlux_Error_Calculation(lb, step));
					cout<<"Error for heat flux = "<<error_HeatFlux.back()<<endl;
					OutputCSV(error_HeatFlux, "Heat Flux Error.csv"); //Output Heat Flux Error
				#endif
				cout<<endl;
			#elif defined Couette_Flow
				total_data++; Postprocessing(lb, step);

				//Error calculation and output partial solution : 
				error_Velocity.push_back(Velocity_Error_Calculation(lb, step));
				cout<<"Error for velocity "<<" = "<<error_Velocity.back()<<endl;
				OutputCSV(error_Velocity, "Velocity Error.csv");

				#if defined MORE_DETAILS
					error_Viscous.push_back(Viscous_Error_Calculation(lb, step));
					cout<<"Error for viscous stress = "<<error_Viscous.back()<<endl;
					OutputCSV(error_Viscous, "Viscous stress Error.csv"); //Output Viscous stress Error
				#endif

				#ifdef ENERGY
					error_Temperature.push_back(Temperature_Error_Calculation(lb, step));
					cout<<"Error for temperature = "<<error_Temperature.back()<<endl;
					OutputCSV(error_Temperature, "Temperature Error.csv"); //Output Temperature Error

					error_HeatFlux.push_back(HeatFlux_Error_Calculation(lb, step));
					cout<<"Error for heat flux = "<<error_HeatFlux.back()<<endl;
					OutputCSV(error_HeatFlux, "Heat Flux Error.csv"); //Output Heat Flux Error
				#endif

				cout<<endl;
			#elif defined Rayleigh_Benard
				#ifdef ENERGY
					uy_max.push_back(Find_Maximum_uy(lb));
					Nusselt_number.push_back(Find_Nusselt_number(lb));

					//Output CSV : 
					name = "Maximum velocity in y direction at Ra = "+ to_string(Ra) + ".csv";
					OutputCSV(uy_max, name); //Output maximum vertical velocity
					name = "Nusselt number at Ra = "+ to_string(Ra) + ".csv";
					OutputCSV(Nusselt_number, name); //Output Nusselt number

					cout<<"Maximum velocity in y direction = "<<uy_max.back()<<endl;
					cout<<"Nusselt number = "<<Nusselt_number.back()<<endl<<endl;
				#endif
				
			#elif defined Cylinder_Tandem_Circle
				total_data++;

				for(int n = 0; n<N; n++) {
					//Obtaining the Cl and Cd right now
					#if defined MOMENTUM_EXCHANGE_ALGORITHM
						CL[n].push_back(MEA_CL(lb, center_x - L/2, center_x + L/2, center_y - n*L - L/2, center_y - n*L + L/2,1,1));
						cout<<"Lift coefficient at time step "<<step<<" = "<<CL[n][total_data-1]<<endl;

						CD[n].push_back(MEA_CD(lb, center_x - L/2, center_x + L/2, center_y - n*L - L/2, center_y - n*L + L/2,1,1));
						cout<<"Drag coefficient at time step "<<step<<" = "<<CD[n][total_data-1]<<endl;
					#endif

					//Output partial solution 
					name = "Value of Cl at LpD = "+to_string(LpD)+" Re = "+to_string(Re) + " for circle " + to_string(n) + ".csv";
					OutputCSV(CL[n], name); //Output CL
					name = "Value of Cd at LpD = "+to_string(LpD)+" Re = "+to_string(Re) + " for circle " + to_string(n) + ".csv";
					OutputCSV(CD[n], name); //Output CD

				}
				cout<<"\n";
			#elif defined Cylinder_Tandem_Triangle
				total_data++;

				for(int n = 0; n<N; n++) {
					//Obtaining the Cl and Cd right now
					#if defined MOMENTUM_EXCHANGE_ALGORITHM
						CL[n].push_back(MEA_CL(lb, tip_x + n*L - D/4, tip_x + n*L + (D+L)/2, tip_y - D, tip_y + D,1,1));
						cout<<"Lift coefficient at time step "<<step<<" = "<<CL[n][total_data-1]<<endl;

						CD[n].push_back(MEA_CD(lb, tip_x + n*L - D/4, tip_x + n*L + (D+L)/2, tip_y - D, tip_y + D,1,1));
						cout<<"Drag coefficient at time step "<<step<<" = "<<CD[n][total_data-1]<<endl;
					#endif

					//Output partial solution 
					name = "output Tandem Triangle/Value of Cl at LpD = "+to_string(LpD) + " Re = "+ to_string(Re)+ " for triangle " + to_string(n) + ".csv";
					OutputCSV(CL[n], name); //Output CL
					name = "output Tandem Triangle/Value of Cd at LpD = "+to_string(LpD) + " Re = "+ to_string(Re)+ " for triangle " + to_string(n) + ".csv";
					OutputCSV(CD[n], name); //Output CD

				}
				cout<<"\n";
			#elif defined Sphere_BB
				total_data++;
				
				#if defined MOMENTUM_EXCHANGE_ALGORITHM
					CL.push_back(MEA_CL(lb, center_x - D, center_x + D, center_y - D, center_y + D, center_z - D, center_z + D));
					cout<<"Lift coefficient at time step "<<step<<" = "<<CL[total_data-1]<<endl;

					CD.push_back(MEA_CD(lb, center_x - D, center_x + D, center_y - D, center_y + D, center_z - D, center_z + D));
					cout<<"Drag coefficient at time step "<<step<<" = "<<CD[total_data-1]<<endl;
				#endif
				

				//Output partial solution 
				name = "Value of Cl at Re = "+to_string(Re) + ".csv";
				OutputCSV(CL, name); //Output CL
				name = "Value of Cd at Re = "+to_string(Re) + ".csv";
				OutputCSV(CD, name); //Output CD
			#elif defined Cylinder_IBM or defined Cylinder_Inline_Oscillating or defined Cylinder_Transverse_Oscillation or defined Cylinder_IBM_Moving_Arbie
				cout<<"Lift coefficient at time step "<<step<<" = "<<CL.back()<<endl;
				cout<<"Drag coefficient at time step "<<step<<" = "<<CD.back()<<endl;	
				// cout<<"V_CM[0] = "<<V_CM[0]<<endl;
				// cout<<"V_CM[1] = "<<V_CM[1]<<endl;
			#endif
			
			#if defined MORE_DETAILS
				lb.Sup_MacroProp(step);
			#endif
			
			// cout<<"Velocity = "<<lb.fluid[int(Nx/2)][int(Ny/2)][int(Nz/2)].ux<<", "<<lb.fluid[int(Nx/2)][int(Ny/2)][int(Nz/2)].uy<<", "<<lb.fluid[int(Nx/2)][int(Ny/2)][int(Nz/2)].uz<<endl<<endl;
		}

		
		if(step == numstep) {
			//Post processing
			#if defined Couette_Flow or defined Poiseuille_Flow
				Postprocessing(lb, step);
			#endif

			cout<<"Want to continue? (Y/N) "<<endl<<endl; 
			char ans; cin>>ans;
			if(ans == 'Y' or ans == 'y') {
				int add; 
				cout<<"How many timesteps that want to be added? "; cin>>add;
				numstep += add;
			}
		}
		
		
	}

	//--------------------------------------------------MAKING CSV ------------------------------------------
	#if defined Cylinder_Basic
		name = "Value of Cl at Re = "+to_string(Re) + ".csv";
		OutputCSV(CL, name); //Output CL
		name = "Value of Cd at Re = "+to_string(Re) + ".csv";
		OutputCSV(CD, name); //Output CD
	#elif defined Cylinder_Tandem_Circle
		for (int n = 0; n<N; n++) {
			name = "Value of Cl at LpD = "+to_string(LpD)+" Re = "+to_string(Re) + " for circle " + to_string(n) + ".csv";
			OutputCSV(CL[n], name); //Output CL
			name = "Value of Cd at LpD = "+to_string(LpD)+" Re = "+to_string(Re) + " for circle " + to_string(n) + ".csv";
			OutputCSV(CD[n], name); //Output CD
		}
	#elif defined Cylinder_Tandem_Triangle
		for (int n = 0; n<N; n++) {
			name = "output Tandem Triangle/LpD = "+to_string(LpD)+" Re = "+to_string(Re)+"/Value of Cl "+ " for triangle " + to_string(n) + ".csv";
			OutputCSV(CL[n], name); //Output CL
			name = "output Tandem Triangle/LpD = "+to_string(LpD)+" Re = "+to_string(Re)+"/Value of Cd "+ " for triangle " + to_string(n) + ".csv";
			OutputCSV(CD[n], name); //Output CD
		}
	#endif


	//-------------------------------------------------------------------------------------------------
	clock_t c_end = clock();
	long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
}
