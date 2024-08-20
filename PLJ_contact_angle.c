/*
Russo Andrea && Anastasio Domenico.
*/

//We're including all the libraries needed to solve the problem.
//The first one solves the incompressible, variable-density Navier-Stokes Equation. 
#include "navier-stokes/centered.h"
//The second one lets us use a two-phase approach, usually used for two fluids separated by an interface. 
#include "two-phase.h"
//The third one lets us account for surface tension. 
#include "tension.h"
//The fourth one lets us know continuous information about the running solver. 
#include "profile6.h"
//The fifth one lets us implement momentum-conserving VOF advection of the velocity component, needed for our two-phase solver. 
#include "navier-stokes/conserving.h"
//The sixth one lets us implement the contact angle through a boundary condition as we'll see. 
#include "contact.h"


//Defining problem variables.
#define my_rho_l 997.
#define my_rho_g 9.97
#define my_mu_g 0.0000184
#define my_sigma 0.00603
#define my_g 9.81
#define my_H 0.0015
#define my_L 0.075
#define my_U 0.492
#define my_r_rho (my_rho_g)/(my_rho_l)
#define my_r_mu (my_mu_g)/(my_mu_l)

//We define all the adimensional coefficients that will be needed. 
#define my_RE (my_rho_l*my_U*my_H)/(2*my_mu_l) 
#define my_WE (my_rho_l*pow(my_U,2)*my_H)/(2*my_sigma) 
#define my_FR (pow(my_U,2))/(my_g*my_L)
#define my_eps (my_H)/(my_L)
#define my_p_ref my_rho_l*my_U*my_U

//This is the time at which we want the sinusoidal forcing at the inlet to start. 
#define my_t_start 0.2

#define my_u_p 1.5*my_U

//This is the time at which the simulation will end. 
#define my_tEnd 10.0

//These are the maximum and minimum level of refinement that will be used in the adapt_wavelet() function. 
#define my_maxlev 11
#define my_minlev 5

//That's the power to which we will elevate 2, and that will be our initial grid size. 
#define my_initlev 9
#define my_iter 5000
#define n_cells 512 

//These are the amplitudes and the frequency of the sinusoidal forcing that will be applied at the inlet. 
#define my_A 0.06
#define my_f 25

//This is the Strouhal number, which indicates the frequence of vortex shedding. 
#define my_Stf (my_f*my_H)/my_U

//We define the strings containing the file's name.
char Nomefile3[50];
char Nomefile4[60];
char Nomefile5[60];
char Nomefile6[60];
char Nomefile7[60];
char Nomefile11[60];
char Nomefile20[60];
char Nomefile21[60];
char Nomefile22[60];
char Nomefile23[60];
char Nomefile24[60];
char Nomefile30[60];

//We define the strings containing the directories' name. 
char Nomedir[60];
char Nomedir2[60];

//We initialize the field containing the height-function.
vector h[];
//We define the vectors containing all the viscosities and superficial tension coefficients for which we will iterate. 
double mu_l_Vec[3] = {0.00089, 0.00032, 0.001};
double sigma_Vec[3] = {0.0725, 0.02, 0.04};

//We initialize the contact angle variable.
double theta0;

//We initialize the counter variables that will be used to iterate.  
int k;
int z; 

//All of these conditions are applied on the left, which should be the inlet, since the gravity is applied from left to right. 
//We impose our boundary conditions. 
p[left] = neumann(0.);

//We will impose parabolic (Poisseuille) flow along the x-directions, but only inside the liquid phase, so it will be between 
//0.5*my_H and -0.5*my_H
u.n[left]  = fabs(y) <= 0.5*my_H ? dirichlet(-my_u_p*sq(y)/sq(0.5*my_H)+my_u_p) : dirichlet(0.);

//We will impose sinusoidal forcing at the inlet, but only inside the liquid phase. 
u.t[left]  = fabs(y) <= 0.5*my_H ? dirichlet(my_A*my_U*sin(2.*pi*my_f*t)*(t > my_t_start)) : dirichlet(0.);

//We will impose that the liquid phase is only present between 0.5*my_H and -0.5*my_H, thus changing the value of the tracers to 1. 
f[left]    = fabs(y) <= 0.5*my_H ? dirichlet(1.) : dirichlet(0.);

//On the right (which should be the lower edge), standard outflow conditions are imposed. 
u.n[right] = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
p[right]   = dirichlet(0.);
f[right] = neumann(0.);
u.t[right] = u.t[] > 0. ? neumann(0.) : 0.;

//For the right and left edges neumann conditions are imposed. 
u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
f[top] = neumann(0.);
p[top] = neumann(0.);
u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);
f[bottom] = neumann(0.);
p[bottom] = neumann(0.);

int main (int argc, char * argv[])
{
  
  //We define the initial grid dimension.
  init_grid (1 << my_initlev); 

  //We define the domain size. 
  size (my_L);

  //We move the origin to the middle (along y-direction) of our box, on the left. The origin before was in the top left corner.
  origin (0., -my_L/2);

  
  //These are the densities that the solver will use for both of our phases. 
  rho1 = my_rho_l, rho2 = my_rho_g;
  TOLERANCE = 1e-4;
 
 //We make a for loop to solve for different contact angles, different liquid viscosities, and different superficial 
 //tension coefficients. 
 for (k = 0; k < 1 ; k++ )
 for (z = 0; z < 1; z++ )
 for (theta0 = 40; theta0 <= 140; theta0 += 100)
 {
  //The field containing the tracers will need the superficial tension coefficient to build the interface.
  f.sigma = sigma_Vec[z];
  
  //This is the viscosity the field will use. 
  mu1 = mu_l_Vec[k];
  
  //We modify the height-function field using the contact_angle() macro, thus imposing the contact angle.
  h.t[left] = (fabs(y) >= my_H*0.5) ? contact_angle(theta0*pi/180.) : 0.;
  
  //We then have to pass the information to the tracers field. 
  f.height = h;
  
  //We give the two directories a name based on the Reynolds and Webers. 
  sprintf(Nomedir,"Re_%0.f_We_%0.f", my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]));
  sprintf(Nomedir2,"Re_%0.f_We_%0.f/Snapshot_theta_%0.f", my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]), theta0);

  //We make the directories. 
  if (mkdir(Nomedir, 0777) == -1)
        printf("Error");
    else
        printf("\nDirectory created\n");

  if (mkdir(Nomedir2, 0777) == -1)
        printf("Error");
    else
        printf("\nDirectory created\n");
      
  //We make the names for the different files that are going to be created, the first six are a zoom at the inlet for each variable.
  sprintf(Nomefile3,"Re_%0.f_We_%0.f/ux_jet_theta_%0.f_zoom.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile4,"Re_%0.f_We_%0.f/uy_jet_theta_%0.f_zoom.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile5,"Re_%0.f_We_%0.f/ux_theta_%0.f_zoom.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile6,"Re_%0.f_We_%0.f/uy_theta_%0.f_zoom.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile7,"Re_%0.f_We_%0.f/f_theta_%0.f_zoom.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile11,"Re_%0.f_We_%0.f/grid_%0.f.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  //We then name all the .mp4s that will contain the fields. 
  sprintf(Nomefile20,"Re_%0.f_We_%0.f/ux_jet_theta_%0.f.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile21,"Re_%0.f_We_%0.f/uy_jet_theta_%0.f.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile22,"Re_%0.f_We_%0.f/ux_theta_%0.f.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile23,"Re_%0.f_We_%0.f/uy_theta_%0.f.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf(Nomefile24,"Re_%0.f_We_%0.f/f_theta_%0.f.mp4",my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]),theta0);
  sprintf (Nomefile30, "Re_%0.f_We_%0.f/", my_H*my_rho_l*my_U/(mu_l_Vec[k]*2),(my_rho_l*pow(my_U,2)*my_H)/(2*sigma_Vec[z]));
  
  //We run the solver. 
  run();
}
}

//Event defining the initial conditions, t = 0 is the condition needed for the event to happen. 
event init (t = 0) {
  if (!restore (file = "restart")) {

  //We use a foreach() to iterate on all the cells in our fields. 
  foreach() {
   //We impose that the liquid phase is only present between 0.5*my_H and -0.5*my_H.
   f[] = fabs(y) <= 0.5*my_H ? 1. : 0.;
   //We impose parabolic (Poisseille) flow in the liquid phase, along the x direction. 
   u.x[] = fabs(y) <= 0.5*my_H ? (-my_u_p*sq(y)/sq(0.5*my_H)+my_u_p) : 0.;
   //We impose y-component equal to zero. 
   u.y[] = 0.; 
   p[] = 0;
            }
  
  //We use the boundary function to impose the BCs also on the fields defined in the CIs. 
  boundary({f,u.x,u.y,p});

  //We then save the level field in our field l[].
  int n1=pow(2,my_maxlev);
  scalar l[];
  foreach()
  l[] = level;
  boundary({l});

  //We then define these new fields, only containing the liquid phase's velocities and pressure. 
  scalar u_jet[],v_jet[],p_jet[],y_jet[];
  foreach() {
       u_jet[] = f[]*u.x[];
       v_jet[] = f[]*u.y[];
       p_jet[] = f[]*p[];
       y_jet[] = f[]*y;
            }    
 //We apply boundaries on our newly defined fields. 
boundary({u_jet,v_jet,p_jet,y_jet});

 //We put them all in a scalar list. 
scalar * list = {f,u_jet,v_jet,p_jet,y_jet};

//We the initiate profiling for our solver. 
char name[50];
sprintf (name, "1Dt_%g.dat",0*pow(my_tEnd,-1));
profile_equi(list, x,n_cells,fname = name);
profile(list, x,fname = name);
}
}

//This event ends the run.  
event end(t = my_tEnd) {
    printf("\nThe run is finished.\n");
}

//This event will impose the contact angle at every timestep. 
event contact (i++) {

//Just like in the main, we impose the contact_angle using the contact_angle() macro, which in turn modifies the height-function 
//field accordingly. 
h.t[left] = fabs(y) >= (0.5*my_H) ? contact_angle(theta0*pi/180.) : 0.;

//We print the theta variable, to know at which value we're iterating. 
printf("\n%f\n",theta0);

//We then pass the information to the tracers field. 
f.height = h;
}

//With this event we then introduce gravity, it is done at every timestep.  
event acceleration (i++) {
    
    //av[] is the vector that the navier-stokes solver will use to represent the acceleration. In this case, we just add 9.81 
    //every timestep. 
    face vector av = a;
    foreach_face(x)
    av.x[] += my_g;
}

//This event makes a logfile for each run. 
event my_logfile (i++) {

  //We define the file's name according to the Reynolds and the Webers. 
  FILE * p_file;
  char Nomefile10[30];
  sprintf(Nomefile10,"logfile_theta_%0.f_Re_%0.f.txt",theta0,my_H*my_rho_l*my_U/(mu_l_Vec[k]*2));
  
  //We create a pointer to the file. 
  p_file = fopen(Nomefile10, "a+");
  
  //We write inside the file all of the solver parameters. 
  fprintf(p_file, "%d %g %d %d %g %g %g %g %d %d %g %g \n",
	   i,dt, mgp.i,mgu.i,mgp.resb,mgu.resb,mgp.resa,mgu.resa,mgp.nrelax,mgu.nrelax,TOLERANCE,TOLERANCE/(dt*dt));
  fclose(p_file);
  
  //We print the timestep and the time at which we're at. 
fprintf(stdout, "i = %d , t = %g\n", i, t);

}

//We initialize the variables and the fields that will be used inside the next event. 
double du,dun;
scalar U_old[],U[];


//This event lets us check if there's convergence. It takes at each timestep the biggest variation of u.x[] and then puts it 
//in a file. 
event my_checkconv (i=0;i<=my_iter;i++) {
//We initialize the max variable and the udiff field. 
double max = 0.;
scalar udiff[];

//We calculate the velocity's module for each cell. 
foreach() {
U[] = sqrt(sq(u.x[])+sq(u.y[]));
}


//What's the reduction operator? It's an operator that lets us operate on various cells simultaneously (in parallel). It is needed 
//because if, for example, we're looking for the max variable it can happen that a variable is accessed and written in 
//simultaneously, invalidating the result of the operation. Reduction allows to choose the criteria with which we will choose which
//variable to actually write, if a cell is accessed simultaneously. For example, the + means that if two values are being written
//at the same time, what will be written inside the cell is the sum. In our case, with max, we will only pick the bigger value. 

//The first argument is the criteria, after the : is the variable on which the reduction should be done. 

foreach(reduction(max:max)) {

//We calculate the difference between the module now, and the module at the previous timestep. 
udiff[] = fabs(U[] -U_old[]);

//We put this value inside ds. 
double ds = fabs(U[] -U_old[]);

//We confront the value now obtained with the maximum value ever obtained. If it's bigger, then ds will be the new max. 
if (ds > max && fabs(y) < 0.5*my_H)
max = ds;
}
boundary((scalar*){udiff}); 

//Now, the maximum variation was found, and it is put inside du. 
du = max;

//We calculate the average variation. 
dun = normf(udiff).avg;

//We now put the velocity field inside U_old[], so that at the next timestep the value will have changed. 
foreach(reduction(max:max)) 
U_old[] = U[];

//We define the name of the file that will contain the maximum and the average variations for each timestep. 
char Nomefile[30];
sprintf(Nomefile, "CheckConv_Theta_%0.f",theta0);

//We define the name of the file that will contain the parameters for each timestep.
char Nomefile2[30];
sprintf(Nomefile2, "Check_param_Theta_%0.f",theta0);

//We create pointers to the files. 
FILE * fp_dvar = fopen(Nomefile, "a+");
FILE * fp_dpar = fopen(Nomefile2, "a+"); 

//Wr write inside the files.
fprintf (fp_dvar, "%d %g %g %g \n",i,t,du*pow(my_U,-1),dun);
fclose (fp_dvar);
fprintf (fp_dpar, "%g %g %g %g %g %g \n",my_r_rho,my_r_mu,my_RE,my_FR,my_WE,my_tEnd);
fclose (fp_dpar); 

}

//This event produces the outputs. 
event output (i++){

//It saves the grid in a field. 
 scalar l[];
  int n1=pow(2,my_maxlev);
  foreach()
  l[] = level;
  boundary({l});

//We redefine the fields only containing the liquid phase's velocities. 
scalar u_jett[], v_jett[];
foreach(){
  u_jett[] = f[]*u.x[];
  v_jett[] = f[]*u.y[];
}

//We produce an .mp4 containing the grid. 
output_ppm (l, file = Nomefile11);

//We produce all the .mp4s containing the fields, but there's a zoom at the inlet (to check for contact angle).
//The third option, linear = true, allows us to interpolate where the level isnt high enough. 
output_ppm(u_jett, file = Nomefile3, linear = true, box = {{0,-1.3*my_H},{my_L/10, 1.3*my_H}});
output_ppm(v_jett, file = Nomefile4, linear = true, box = {{0,-1.3*my_H},{my_L/10, 1.3*my_H}});
output_ppm(u.x, file = Nomefile5, linear = true, box = {{0,-1.3*my_H},{my_L/5, 2*my_H}});
output_ppm(u.y, file = Nomefile6, linear = true, box = {{0,-2*my_H},{my_L/5, 2*my_H}});
output_ppm(f, file = Nomefile7, box = {{0,-2*my_H},{my_L/10, 2*my_H}});

//We produce all the .mp4s containing the fields.
output_ppm(u_jett, file = Nomefile20, linear = true);
output_ppm(v_jett, file = Nomefile21, linear = true);
output_ppm(u.x, file = Nomefile22, linear = true);
output_ppm(u.y, file = Nomefile23, linear = true);
output_ppm(f, file = Nomefile24, linear = true);

//If the stationary regime is reached (t > my_t_start) every 6000 timestep it will take a snapshot. 
if ((i % 6000) == 0 && t > my_t_start){
//We initialize the strings that will be used to give a name to the Snapshots. 
char Nomefile31[100];
char Nomefile32[100];
char Nomefile33[100];
char Nomefile34[100];

//We put the name inside the strings. 
sprintf (Nomefile31, "Snapshot_theta_%0.f/t_%f.dat", theta0, t);
strcpy (Nomefile32,Nomefile30);
strcat (Nomefile32,Nomefile31);

sprintf (Nomefile33,"time_theta_%0.f", theta0);
strcpy (Nomefile34,Nomefile30);
strcat (Nomefile34,Nomefile33);

//We create a pointer to the file.
FILE * fp_snap = fopen (Nomefile32, "w");
//We output the field. 
output_field ({u.x, u.y, f, p}, fp_snap, n=n1, linear = true,box = {{0.0,-my_H},{my_L,my_H}});
fclose (fp_snap);

//We create the file containing the times at which the snapshots were taken. 
FILE * fp_time = fopen(Nomefile34, "a+");
fprintf (fp_time, "%f\n", t);
fclose (fp_time);

}
}
//This event is needed to  adapt the grid to the tracer's variation. Basically we increase the resolution until the error 
//is less than 1e-5. 
event adapt (i++) {
//To do this, we used the adapt_wavelet function. 
adapt_wavelet ({f}, (double[]){1e-5}, minlevel = my_minlev, maxlevel = my_maxlev);
adapt_wavelet ({u.x}, (double[]){1e-5}, minlevel = my_minlev, maxlevel = my_maxlev);
adapt_wavelet ({u.y}, (double[]){1e-5}, minlevel = my_minlev, maxlevel = my_maxlev);
//We refine our grid. 
refine(x > -1. && fabs(y) < my_H && level < my_maxlev);
}

