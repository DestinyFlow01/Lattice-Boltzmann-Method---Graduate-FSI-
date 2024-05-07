#include<iostream>
#include<cmath>
#include<string>
using namespace std;



// -----------------------------------------DEFINITIONS OF SIMULATION PARAMETERS ----------------------------------
#ifndef DEFINES_H
	#define DEFINES_H
	
	//FLOW CASES (Basic Bounce Back)
	// #define Cylinder_Basic
	// #define Sphere_BB

	//Validation 1 (paper : Characteristics of the flow around four cylinders of various shapes)
	//#define Cylinder_Tandem_Circle
	
	//Validation 2 (paper : A Thermal Lattice Boltzmann Model for Flows with Viscous Heat Dissipation)
	// #define Poiseuille_Flow
	// #define Couette_Flow	

	//Validation 3 : Incorpera
	// #define Rayleigh_Benard	//For natural convection using Boussinesq approximation, energy must be activated

	//For paper 
	// #define Cylinder_Tandem_Triangle

	//Flow Cases (IBM)
	#define IBM
	
	//Method for IBM 
	#ifdef IBM 
		//IBM method
		// #define Explicit_Diffuse_Interface_Scheme
		#define Implicit_Diffuse_Interface_Scheme
	#endif

	//Cases for IBM 
	#ifdef IBM
		#define Cylinder_IBM
		// #define Cylinder_Inline_Oscillating
		// #define Cylinder_Transverse_Oscillation
		// #define Cylinder_IBM_Moving_Arbie 		//Suggestions from pak Rizqi Arbie
	#endif
	

	





	//More detail on fluid properties (pressure and shear stress) : 
	// #define MORE_DETAILS

	//Including energy equation : 
	// #define ENERGY

	//Dimensionalization of the quantity : 
	// #define Dimensional

	//Lift and drag calculation
	#ifdef IBM
		#define Volume_Force
	#else 
		#define MOMENTUM_EXCHANGE_ALGORITHM
	#endif

	//LBM Method : 
	#define Conventional_LBM


	static const double PI = 3.14159265;


	//Parameters for each cases
	//FLOW PARAMETERS
	#define rho0 1.0

	#define dt 1.0
	#define dx 1.0
	#define dy 1.0
	#define dz 1.0

	static const double cs2 = pow(dx,2)/(3*pow(dt,2));
	#define p0 rho0*cs2

	//Setup for each cases
	//Bounce back
	#if defined Cylinder_Basic
		#define Re 200.0
		static const double D = 20;

		static const int Nx = 100;//40*D;
		static const int Ny = 100;//10*D;
		static const int Nz = 1;

		static const double u_max = 0.1;
		static double center_x = dx*(Nx/5.0);
		static double center_y = dy*(Ny/2.0);

		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;

		#if defined Dimensional
			static const double D_dim = 1.;	//Diameter of the cylinder
			static const double V_dim = 1; //Freestream velocity
			static const double rho_dim = 1.225; //Air density

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream
			#endif
		#endif

	#elif defined Poiseuille_Flow
		#define Re 100.0

		static const int Nx = 5;
		static const int Ny = 100;
		static const int Nz = 5;

		static const double tau = 0.9;
		static const double nu = cs2*(tau - 0.5);

		static const double u_max = Re*nu/(Ny+1); //Maximum velocity

		//Forcing function 
		static double dp_dx = -8*rho0*nu*u_max/pow(Ny+1,2);

		#define Velocity_Ratio
		#define Temperature_Ratio

		#if defined Dimensional
			static const double D_dim = 1.;	//Width of the channel (m)
			static const double V_dim = 10; //Wall velocity (m/s)
			static const double rho_dim = 1000; //Air density (kg/m3)
			static const double c_dim = 340; //Speed of sound (m/s)

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream (Pa)
			#endif

			#ifdef ENERGY
				static const double Tb_dim = 300; 	//Bottom wall temperature (K)
				static const double cv_dim = 4200; 	//Heat capacity (J/kgK)
			#endif
		#endif

		#ifdef ENERGY
			static const double Pr = 0.7;

			static const double alpha = nu/Pr; //Thermal diffusivity
			static const double tau_g = 0.5 + alpha/(cs2*dx); 
			static const double k = cs2*tau_g*dx;	//Thermal conductivity
			static const double cp = k/(rho0*alpha); //specific heat capacity

			static const double Tt = 1;  //top wall temperature
			static const double Tb = 1;  //bottom wall temperature
			static const double j_dot_n = 0;
		#endif


	#elif defined Couette_Flow
		#define Re 300

		static const int Nx = 2;
		static const int Ny = 120;
		static const int Nz = 2;

		static const double tau = 0.9;
		static const double nu = cs2*(tau - 0.5);

		static const double u_wall = Re*nu/(Ny+1);  //Top wall velocity
		static const double u_max = u_wall;

		//For velocity boundary on top wall : 
		static double uwall_topx = u_wall;

		//For velocity boundary on bottom wall : 
		static double uwall_botx = 0;

		//Publishing the type of results of velocity and temperature
		#define Velocity_Ratio		
		#ifdef ENERGY
			#define Temperature_Ratio
		#endif
		//Prioritize the velocity and temperature ratio result, if not, prioritize dimensional value, if not output the simulation result

		#if defined Dimensional
			static const double D_dim = 1.;	//Width of the channel (m)
			static const double V_dim = 10; //Wall velocity (m/s)
			static const double rho_dim = 1000; //Air density (kg/m3)
			static const double c_dim = 340; //Speed of sound (m/s)

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream (Pa)
			#endif

			#ifdef ENERGY
				static const double Tb_dim = 300; 	//Bottom wall temperature (K)
				static const double cv_dim = 4200; 	//Heat capacity (J/kgK)
			#endif
		#endif

		#ifdef ENERGY
			static const double Pr = 0.7;	//Prandtl number
			static const double Br = 10;		//Brinkmann number 
			static const double Ec = Br/Pr;	//Eckart number 

			static const double alpha = nu/Pr; //Thermal diffusivity
			static const double tau_g = 0.5 + alpha/(cs2*dx); 
			static const double k = cs2*tau_g*dx;	//Thermal conductivity
			static const double cp = k/(rho0*alpha); //specific heat capacity

			static const double Tb = 0.5;	//Bottom wall temperature
			static const double dT = pow(u_wall,2)/(cp*Ec);	//Temperature difference between top and bottom wall
			static const double Tt = Tb + dT;	//Top wall temperature
			static const double j_dot_n = 0;
		#endif
		
	#elif defined Rayleigh_Benard
		#ifdef ENERGY
			static const int Nx = 200;		//w in incorpera
			static const int Ny = 50;		//H in incorepera
			static const int Nz = 5;		//L in incorepera

			static const double Tb = 2;	//Bottom wall temperature
			static const double dT = 1;	//Temperature difference
			static const double Tt = Tb - dT;	//Top wall temperature
			static const double j_dot_n = 0;	//For both top and bottom wall

			static const double Ra = 1710; //Rayleigh number to find nu
			static const double g_beta = pow(10,-5);
			static const double Pr = 0.71;	//Prandtl number
			static const double nu = sqrt(g_beta*dT*pow(Nz+1,3)*Pr/Ra);
			static const double tau = 0.5 + nu/cs2;	
			static const double alpha = nu/Pr; //Thermal diffusivity
			static const double tau_g = 0.5 + alpha/(cs2*dx); 
			static const double k = cs2*tau_g*dx;	//Thermal conductivity
			static const double cp = k/(rho0*alpha); //specific heat capacity
		#endif
	#elif defined Cylinder_Tandem_Circle
		#define Re 200.0

		static int D = 20;		//Diameter of circle
		static const int Nx = 800;
		static const int Ny = 400;
		static const int Nz = 1;

		//Geometry properties
		static double LpD = 1.5;	//Aspect ratio
		static double L = D*LpD; 	//Separation distance of each circle center
		static int N = 2;		//Number of circle
		static int center_x = Nx/5;			//Center of the upper circle x coor
		static int center_y = Ny/2 + L/2;	//Center of the upper circle y coor
				
		//Flow properties
		static double u_max = 0.05;
		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;
	#elif defined Cylinder_Tandem_Triangle
		#define Re 200.0

		static int D = 20;		//Length of equilateral triangle
		static const int Nx = 200;
		static const int Ny = 200;
		static const int Nz = 1;

		//Geometry properties
		static double LpD = 1.5;		//Aspect ratio
		static int N = 3;		//Number of triangle
		static int tip_x = Nx/5;	//Tip of the leftmost triangle x coor
		static int tip_y = Ny/2;	//Tip of the leftmost triangle y coor
		static double L = D*LpD; 	//Separation distance of each equilateral center of mass
		
		//Flow properties
		static double u_max = 0.05;
		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;

		#ifdef Dimensional
			static const double D_dim = 1.;	//Width of the channel (m)
			static const double V_dim = 10; //Wall velocity (m/s)
			static const double rho_dim = 1000; //Air density (kg/m3)
			static const double c_dim = 340; //Speed of sound (m/s)

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream (Pa)
			#endif
		#endif
	
	#elif defined Sphere_BB
		#define Re 200.0
		static const double D = 20;

		static const int Nx = 800;
		static const int Ny = 400;
		static const int Nz = 400;

		static const double u_max = 0.1;
		static double center_x = dx*(Nx/5.0);
		static double center_y = dy*(Ny/2.0);
		static double center_z = dz*(Nz/2.0);

		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;

		#if defined Dimensional
			static const double D_dim = 1.;	//Diameter of the cylinder
			static const double V_dim = 1; //Freestream velocity
			static const double rho_dim = 1.225; //Air density

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream
			#endif
		#endif

	//IBM
	#elif defined Cylinder_IBM
		//Flow properties
		#define Re 40.0
		static const double D = 20;

		static const int Nx = 40*D;
		static const int Ny = 20*D;
		static const int Nz = 1;

		static const double u_max = 0.1;
		static double center_x = dx*(Nx/5.0);
		static double center_y = dy*(Ny/2.0);

		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;

		//Dimensional quantity
		#if defined Dimensional
			static const double D_dim = 1.;	//Diameter of the cylinder
			static const double V_dim = 1; //Freestream velocity
			static const double rho_dim = 1.225; //Air density

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream
			#endif
		#endif

		//IBM properties
		static const int number_nodes = 200;
		static const int level_marker = 1;

		#ifdef Implicit_Diffuse_Interface_Scheme
			static int m_max = 20;
		#endif

		//For added mass : 
		static const double Vsolid = PI*D*D*0.25;

	#elif defined Cylinder_Inline_Oscillating
		#define Re 100.0
		#define KC 5.0
		static const double D = 20;

		static const int Nx = 30*D;
		static const int Ny = 20*D;
		static const int Nz = 1;

		//External movement : 
		static const double period = 7200;
		static const double frequency = 1/period;
		static const double u_max = KC*D*frequency;
		static const double amplitude = u_max/(2*PI*frequency);

		static double center_x = dx*(Nx/2.0);
		static double center_y = dy*(Ny/2.0);

		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;

		//For added mass : 
		static const double Vsolid = PI*D*D*0.25*(Nz+2);

		//Dimensional quantity
		#if defined Dimensional
			static const double D_dim = 1.;	//Diameter of the cylinder
			static const double V_dim = 1; //Freestream velocity
			static const double rho_dim = 1.225; //Air density

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream
			#endif
		#endif

		//IBM properties
		static const int number_nodes = 204;
		static const int level_marker = 1;

		#ifdef Implicit_Diffuse_Interface_Scheme
			static int m_max = 20;
		#endif
	#elif defined Cylinder_Transverse_Oscillation
		#define Re 185.0
		static const double D = 20;

		static const int Nx = 800;
		static const int Ny = 400;
		static const int Nz = 1;

		//External motion : 
		static const double u_max = 0.05;
		static const double ratio_Amd = 0.2;
		static const double amplitude = ratio_Amd*D;
		static const double Strouhal_number = 0.192908;
		static const double frequency_0 = Strouhal_number*u_max/D;
		static const double frequency_ratio = 1.2;
		static const double frequency_e = frequency_0*frequency_ratio;

		static double center_x = dx*(Nx/5.0);
		static double center_y = dy*(Ny/2.0);

		static const double nu = u_max*D/Re;
		static const double tau = 0.5*dt+nu/cs2;

		//For added mass : 
		static const double Vsolid = PI*D*D*0.25*(Nz+2);

		//Dimensional quantity
		#if defined Dimensional
			static const double D_dim = 1.;	//Diameter of the cylinder
			static const double V_dim = 1; //Freestream velocity
			static const double rho_dim = 1.225; //Air density

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream
			#endif
		#endif

		//IBM properties
		static const int number_nodes = 204;
		static const int level_marker = 1;

		#ifdef Implicit_Diffuse_Interface_Scheme
			static int m_max = 20;
		#endif
	#elif defined Cylinder_IBM_Moving_Arbie
		//Flow properties
		#define Re 40.0
		static const double D = 40;
		static const double target_tstar = 50;

		static const int Nx = (16 + target_tstar)*D;
		static const int Ny = 20*D;
		static const int Nz = 1;

		static const double u_max = 0.1;
		static double center_x = Nx - 8*D;
		static double center_y = dy*(Ny/2.0);

		static const double nu = u_max*D/Re; //problem here, no need to change the computational method for lift and drag (probably)
		static const double tau = 0.5*dt+nu/cs2;

		

		//Dimensional quantity
		#if defined Dimensional
			static const double D_dim = 1.;	//Diameter of the cylinder
			static const double V_dim = 1; //Freestream velocity
			static const double rho_dim = 1.225; //Air density

			#if defined MORE_DETAILS
				static const double pressure_dim = 101325;	//Pressure at freestream
			#endif
		#endif

		//IBM properties
		static const int number_nodes = 200;
		static const int level_marker = 1;

		#ifdef Implicit_Diffuse_Interface_Scheme
			static int m_max = 5;
		#endif

		//For added mass : 
		static const double Vsolid = PI*D*D*0.25;
	#endif


	


	//OUTPUT VALUE
	#if defined Cylinder_IBM_Moving_Arbie
		#define T_OUT D*target_tstar/u_max
		#define Nt 10
	#else
		#define T_OUT 30000
		#define Nt 10
	#endif
	
	
	//VELOCITY SET bb
	#define D2Q9
	// #define D3Q19
	
	//BOUNDARY CONDITION
	#define TYPE_F 0 // fluid domain
	//Momentum solver : 
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_OUT 2 // outflow boundary
	#define TYPE_IN 3 //Inflow boundary
    #define TYPE_P 4 // periodic boundary
	#define TYPE_L 5 // Slip boundary

	//Energy solver : 
	#define TYPE_T 6 //Dirichlet boundary condition (Temperature defined)
	#define TYPE_Q 7 //Neumann boundary condition (Heat Flux defined)
	
#endif
// ----------------------------------------------------------------------------------------------------------------------






//-------------------------------DEFINITION FOR LBM SIMULATION----------------------------------------------------------
#ifndef LBM_H
#define LBM_H
	//Velocity set 
	#if defined D2Q9 
		const int ndim = 2;
		const int npop = 9;
		//Velocity set 
		const int cx[npop] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
		const int cy[npop] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
		const int cz[npop] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
		//weight function 
		const double w[npop] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};
		//opposite for D2Q9 velocity set : 
		const int opposite[npop] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; //For No slip boundary condition
		////////////////////////{0, 1, 2, 3, 4, 5, 6, 7, 8};
		const int slip[npop] = {0, 1, 4, 3, 2, 7, 8, 5, 6}; //For slip boundary condition
	#elif defined D3Q19
		const int ndim = 3;
		const int npop = 19;
		//Velocity set 
		const int cx[npop] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0};
		const int cy[npop] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
		const int cz[npop] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
		//weight function 
		const double w[npop] = {1.0/3, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/36};
		//opposite for D2Q9 velocity set : 
		const int opposite[npop] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17}; //For No slip boundary condition
		////////////////////////{0, 1, 2, 3, 4, 5, 6, 7, 8};
		//const int slip[npop] = {0, 1, 4, 3, 2, 7, 8, 5, 6}; //For slip boundary condition
	#endif
	
	//lattice definition 
	class LATTICE{
		public : 
			short type = TYPE_F; //details in setup
			double density = 1.0;
			double ux = 0.0;
			double uy = 0.0;
			double uz = 0.0;
			double f[npop], f_star[npop];
			double F[3] = {0,0,0}; //External force

			#ifdef Implicit_Diffuse_Interface_Scheme
				double Fm[3] = {0,0,0};
			#endif

			#if defined MORE_DETAILS
				double pressure = 0.0; 
				double sigma[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};	//viscous stress
				double vorticity[3] = {0,0,0};

				//Time Averaged value : 
				double avg_velocity[3] = {0,0,0};
				double avg_pressure = 0;
				double avg_vorticity[3] = {0,0,0};
				double avg_ReynoldStress[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};

				//Time Rms value : 
				double rms_velocity[3] = {0,0,0};
				double rms_pressure = 0;
				double rms_vorticity[3] = {0,0,0};
				double rms_ReynoldStress[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
			#endif

			#if defined ENERGY
				double g[npop], g_star[npop];
				double T;	//temperature
				double q[3] = {0,0,0};	//Heat flux
				short type_Energy = TYPE_F;	//Type for energy equation
				double F_vis[npop];		//For viscous dissipation term 
			#endif
	};
	
	//Marker definition, will be used in the form of array
	class MARKER {
		//The definition of marker will be used in the form of level and the order will be CCW wrt theta = 0
		public : 
			//Variables : 
			double x, y, z = 1;  //Current position
			double u_boundary[3] = {0,0,0}; //Current marker velocity
			double u_boundary_des[3] = {0,0,0};	//Desired marker velocity
			double F_boundary[3] = {0,0,0}; //Lagrangian force of marker velocity
			double normal[3];

			#ifdef Implicit_Diffuse_Interface_Scheme
				double Fm_boundary[3] = {0,0,0};
			#endif

			//Finding normal : 
			void Normal_Marker(int n, double x1, double x2, double y1, double y2, double z1, double z2);
	};


	//LBM class definition 
	class LBM {
		private : 
			int Nx, Ny, Nz = 1;
			double tau;
		public : 
			LATTICE ***fluid; //Definition of lattice for fluid (IBM), accessible to all

			LBM(int Nx, int Ny, int Nz, double tau); //Parametric constructor for Lattice
			
			//Mutator
			void Init();
			

			void Collision(); 
			void Streaming(); 
			void MacroProp();
			
			#if defined Conventional_LBM
				void ConvLBM(int i, int j, int k);
				double Calc_feq(int i, int j, int k, int l);

				#ifdef ENERGY
					double Calc_geq(int i, int j, int k, int l, short type, vector<double> uwall);
				#endif
			#endif

			#ifdef MORE_DETAILS
				void Sup_MacroProp(int n);
			#endif
		
			void Determine_Wall_Velocity(int i, int j, int k, vector<double>& uwall);
			#ifdef ENERGY
				//For Dirichlet Boundary Condition
				double Determine_Wall_Temperature(int i, int j, int k);	
				
				//For Neumann Boundary Condition
				double Determine_Twall_Neumann(int i, int j, int k, vector<double> uwall, vector<int> cdotn_pos, vector<int> cdotn_notpos, double j_dot_n);
			#endif

			//Accessor
			int getNx() const;
			int getNy() const;
			int getNz() const;
			double getTAU() const;
	};
	
#endif
// ----------------------------------------------------------------------------------------------------------------------





//-----------------------------------STEPS IMMERSED BOUNDARY METHODS ---------------------------------------------
#ifndef IBM_H
#define IBM_H
	//Determine nearest lattice coordinate (ijk) from a marker
	void NearestLatticeCoordinate(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int& i, int& j, int& k, vector<MARKER>& marker, int n);

	//Kernel (Discretized Dirac delta) for Stencil = 4 :
	double kernel(double r);

	//Calculate average density for added mass contribution
	double Calc_AverageDensity(LBM& lb, vector<MARKER>& marker, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);


	//Steps for IBM (reference from Kang et al Disseration for IB-LBM): 
	void IBLBM(LBM& lb, vector<MARKER>& marker, vector<double>& ds, vector<double>& Cl, vector<double>& Cd, double V_CM[], double t);
#endif
//---------------------------------------------------------------------------------------------------------------





//DONE FOR IBM
// --------------------------------------SETUP OF THE PROBLEM -----------------------------------------------------------
#ifndef SETUP_H
	#define SETUP_H
	
	#if defined Cylinder_Basic
		LBM main_setup_Cylinder();
	#elif defined Poiseuille_Flow
		LBM main_setup_Poiseuille();
	#elif defined Couette_Flow
		LBM main_setup_Couette();
	#elif defined Rayleigh_Benard
		LBM main_setup_RayleighBenard1();
	#elif defined Cylinder_Tandem_Circle
		LBM main_setup_Tandem_Circle();
	#elif defined Cylinder_Tandem_Triangle
		LBM main_setup_Tandem_Triangle();
	#elif defined Sphere_BB
		LBM main_setup_Sphere();
	#elif defined Cylinder_IBM
		LBM main_setup_Cylinder_IBM();
	#elif defined Cylinder_Inline_Oscillating
		LBM main_setup_Cylinder_InlineOscillating();
	#elif defined Cylinder_Transverse_Oscillation
		LBM main_setup_Cylinder_TransverseOscillating();
	#elif defined Cylinder_IBM_Moving_Arbie
		LBM main_setup_Cylinder_IBM_Moving();
	#endif

#endif


//-----------------------------------------------------------------------------------------------------------------------










//DONE FOR GEOMETRY DEFINITION IN IBM
// --------------------------------------DEFINITION OF PROBLEM GEOMETRY -------------------------------------------------
#ifndef GEOM_H
	#define GEOM_H
	
	#if defined Cylinder_Basic 
		void Generate_Cylinder(LBM& lb, int& count);
	#elif defined Cylinder_Tandem_Circle
		void Generate_Circle(LBM& lb, int& count, int center_x, int center_y);
	#elif defined Cylinder_Tandem_Triangle
		void Generate_Triangle(LBM& lb, int& count, int tip_x, int tip_y);
	#elif defined Sphere_BB
		void Generate_Sphere(LBM& lb, int& count);
	#elif defined Cylinder_IBM
		void Generate_Cylinder_Marker(vector<MARKER>& marker);
	#elif defined Cylinder_Inline_Oscillating
		void Generate_Cylinder_Marker(vector<MARKER>& marker);
	#elif defined Cylinder_Transverse_Oscillation
		void Generate_Cylinder_Marker(vector<MARKER>& marker);
	#elif defined Cylinder_IBM_Moving_Arbie
		void Generate_Cylinder_Marker(vector<MARKER>& marker);
	#endif
	
#endif

// ----------------------------------------------------------------------------------------------------------------------









// -------------------------------------------OUTPUT HEADER--------------------------------------------------------------
#ifndef OUTPUT_H
	#define OUTPUT_H
	void print_Logo();
	void print_Info();
	void OutputVTK(int&, LBM&); //For fluid
	
	#ifdef MORE_DETAILS
		void OutputAVGVTK(int&, LBM&); //For fluid	
	#endif

	#if defined IBM
		void OutputVTK_Marker(int&, vector<MARKER>&); // For boundary marker
	#endif

	//Force calculation for non-IBM method  : 
	//Momentum Exchange Algorithm for stationary boundary : 
	#if defined MOMENTUM_EXCHANGE_ALGORITHM
		double MEA_CL(LBM &lb, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
		double MEA_CD(LBM &lb, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
	#endif
	
	
	//Analytical solution : 
	#if defined Couette_Flow
		double Couette_Analytic_Velocity(double t, double y);
		double Velocity_Error_Calculation(LBM& lb, int step);

		#if defined MORE_DETAILS
			double Couette_Analytic_Viscous(double t, double y);
			double Viscous_Error_Calculation(LBM& lb, int step);
		#endif

		#ifdef ENERGY
			//Temperature : 
			double Couette_Analytic_Temperature(double t, double y);
			double Temperature_Error_Calculation(LBM& lb, int step);

			//Heat Flux : 
			double Couette_Analytic_HeatFlux(double t, double y);
			double HeatFlux_Error_Calculation(LBM& lb, int step);
		#endif
	#elif defined Poiseuille_Flow
		double Poiseuille_Analytic_Velocity(double t, double y);
		double Velocity_Error_Calculation(LBM& lb, int step);

		#if defined MORE_DETAILS
			double Poiseuille_Analytic_Viscous(double t, double y);
			double Viscous_Error_Calculation(LBM& lb, int step);
		#endif

		#ifdef ENERGY
			double Poiseuille_Analytic_Temperature(double t, double y);
			double Temperature_Error_Calculation(LBM& lb, int step);

			double Poiseuille_Analytic_HeatFlux(double t, double y);
			double HeatFlux_Error_Calculation(LBM& lb, int step);
		#endif
	#elif defined Rayleigh_Benard
		#ifdef ENERGY
			double Find_Maximum_uy(LBM& lb);
			double Find_Nusselt_number(LBM& lb);
		#endif
			

	#endif

	//Output CSV for the data
	void OutputCSV(vector<double>& error, string name); 
	
	//Post processing purpose : 
	#if defined Couette_Flow or defined Poiseuille_Flow
		void Postprocessing (LBM& lb, int step);
	#endif
        
#endif

//-----------------------------------------------------------------------------------------------------------------------
