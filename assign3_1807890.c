/*********************************************************
* This program solves coupled equations which describes
* the diffusion-damped system with constant coefficients.
* A time splitting technique by Marchuk–Strang is used. 
* For two subproblems, A1 and A2 we solve A1 with time length
* tau=(1/2)*dt then solve the other one, A2, for time length
* dt and solve with A1 again for time length tau=(1/2)*dt.
* 
* The read_input function should not be modified, but the
* use of the function in main() may or may not be correct.
*
* To compile: gcc -Wall -Werror -std=c99 -o assign3 assign3.c -lm
*
* List of identified errors:
*--------+--------------------------------------------
*  Line  |     Brief description of a fix
* Number |
*--------+-------------------------------------------
PLEASE NOTE LINE NUMBERS REFER TO THE LINE NUMBERS OF THE
ORIGINAL BROKEN CODE, NOT THE CODE IN THIS FILE
*  24 ...... Include <math.h> to use sin and cos
*  29 ...... Change return type of main method to 'int'
*  43 ...... Change arguments of 'read_input' to pointers
*  45-49.... Fix formula for dx and decrease size of time step for stability
*  58 ...... Declare separate derivative variables for u and v (not strictly 
              an error fix but improves code readability)
* 61-68 .... Fix memory allocation
* 71-75 .... Fix booleans in the if statement checking allocation pointers
* 80 ....... Explicitely initialise ctime to 0
* 83-90 .... Change <= in for loop to < and add enforcement of boundary conditions to
              initial arrays
* 96-97 .... Change rotation factors to cos(dt) and sin(dt) (not 2dt)
* 113-114 .. Fix application of rotation matrix by ensuring signs are correct
* 99-124 ... Rewrite code to be able to calculate 2nd derivative at boundary cases
* 127-133 .. Fix copying of arrays for next time step
* 143-146 .. Fix freeing of memory

********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI  3.141592

void read_input(double *D, double *L, int *nx, double *t_F);

int main(void) {
  /******************************/
  /* Declarations of parameters */
  /******************************/
  /* Number of grid points */                  
  int nx;
  /* Length of domain */
  double L;
  /* Equation coefficient */
  double D;
  /* Length of time to run simulation. */
  double t_F;

  /* Read in from file; */
  read_input(&D, &L, &nx, &t_F);

  /* Grid spacing and small time step to ensure stability */
  double dx = L/(nx-1);
  double invdx2 = 1.0/(dx*dx);      
  double dt = 0.01;

  /************************************************/
  /* Solution Storage at Current / Next time step */
  /************************************************/
  double *uc, *un, *vc, *vn;
  /* Time splitting solutions */
  double *uts1, *uts2, *vts1, *vts2;
  /* Derivatives used in finite difference, one for u and one for v */
  double uderiv, vderiv;

  /* Allocate memory according to size of nx */
  uc = (double*) malloc(nx * sizeof(double));
  un = (double*) malloc(nx * sizeof(double));
  vc = (double*) malloc(nx * sizeof(double));
  vn = (double*) malloc(nx * sizeof(double));
  uts1 = (double*) malloc(nx * sizeof(double));
  uts2 = (double*) malloc(nx * sizeof(double));
  vts1 = (double*) malloc(nx * sizeof(double));
  vts2 = (double*) malloc(nx * sizeof(double));

  /* Check the allocation pointers */
  if (uc == NULL||un == NULL||vc == NULL||vn == NULL||uts1 == NULL||
  uts2 == NULL||vts1 == NULL||vts2 == NULL) {
    printf("Memory allocation failed\n");
    return 1;
  }
  
  /*Declare indexing variable, grid x-location variable and current-time variable*/
  int k;
  double x;
  double ctime = 0.0;

  /* Initialise arrays using the provided initial conditions */
  for(k = 0; k < nx; k++) {
    x = k*dx;
    uc[k]  = 1.0 + sin(2.0*PI*x/L);
    vc[k]  = 0.0;
    /* Set other arrays to 0 */
    uts1[k] = 0.0; uts2[k] = 0.0;
    vts1[k] = 0.0; vts2[k] = 0.0;
  }
  
  /*Enforce periodic boundary conditions for intial arrays*/
  uc[nx - 1] = uc[0];
  vc[nx - 1] = vc[0];

  /* Print initial output (for t=0) */
  for (k = 0; k < nx; k++) {
    x = k*dx;
    printf("%g %g %g %g\n",ctime,x,uc[k],vc[k]);
  }

  /* Declare rotation factors for time-splitting scheme. */
  double cfac = cos(dt);
  double sfac = sin(dt);
  
  /* Loop over timesteps */ 
  while (ctime < t_F){
    
    /* First substep for diffusion equation, A_1 */ 
    for (k = 0; k < nx; k++) {
      x = k*dx;
      /* Calculate 2nd derivative for u and v using the fact that the derivatives 
      are also periodic */
      if(k == 0){
        uderiv = (uc[nx - 1] + uc[1] - 2*uc[0])*invdx2;
        vderiv = (vc[nx - 1] + vc[1] - 2*vc[0])*invdx2;
      }
      else if(k == nx - 1){
        uderiv = (uc[nx - 2] + uc[0] - 2*uc[nx - 1])*invdx2;
        vderiv = (vc[nx - 2] + vc[0] - 2*vc[nx - 1])*invdx2;
      }
      else{
        uderiv = (uc[k-1] + uc[k+1] - 2*uc[k])*invdx2;
        vderiv = (vc[k-1] + vc[k+1] - 2*vc[k])*invdx2;
      }

      /*First diffusion step for u and v*/
      uts1[k] = uc[k] + (D * uderiv * 0.5*dt);
      vts1[k] = vc[k] + (D * vderiv * 0.5*dt);
    }

    /* Second substep for decay/growth terms, A_2 */
    for (k = 0; k < nx; k++) {
      x = k*dx;
      /* Apply rotation matrix to u and v, */
      uts2[k] = cfac*uts1[k] + sfac*vts1[k];
      vts2[k] = -sfac*uts1[k] + cfac*vts1[k];
    }

    /* Third substep for diffusion terms, A_1 */
    for (k = 0; k < nx; k++) {
      x = k*dx;
      /* Calculate 2nd derivative for u and v using the fact that the derivatives 
      are also periodic */
      if(k == 0){
        uderiv = (uts2[nx - 1] + uts2[1] - 2*uts2[0])*invdx2;
        vderiv = (vts2[nx - 1] + vts2[1] - 2*vts2[0])*invdx2;
      }
      else if(k == nx - 1){
        uderiv = (uts2[nx - 2] + uts2[0] - 2*uts2[nx - 1])*invdx2;
        vderiv = (vts2[nx - 2] + vts2[0] - 2*vts2[nx - 1])*invdx2;
      }
      else{
        uderiv = (uts2[k-1] + uts2[k+1] - 2*uts2[k])*invdx2;
        vderiv = (vts2[k-1] + vts2[k+1] - 2*vts2[k])*invdx2;
      }

      /*Third diffusion step for u and v*/
      un[k] = uts2[k] + (D * uderiv * 0.5*dt);
      vn[k] = vts2[k] + (D * vderiv * 0.5*dt);
    }

    /*Enforce periodic boundary conditions for the latest solution*/
    un[nx - 1] = un[0];
    vn[nx - 1] = vn[0];
    
    /* Copy next values at timestep to u, v arrays. */
    double *tmp;
    tmp = uc;
    uc = un;
    un = vc;
    vc = vn;
    vn = tmp;

    /*Increment time*/
    ctime += dt;

    /* Print output */
    for (k = 0; k < nx; k++) {
      x = k*dx;
      printf("%g %g %g %g\n",ctime,x,uc[k],vc[k]);
    }
  }
  
  /* Free allocated memory */
  free(uc); free(un);
  free(vc); free(vn);
  free(uts1); free(uts2);
  free(vts1); free(vts2);

  return 0;
}

// The lines below don't contain any bugs! Don't modify them
void read_input(double *D, double *L, int *nx, double *t_F) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(4!=fscanf(infile,"%lf %lf %d %lf",D,L,nx,t_F)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}