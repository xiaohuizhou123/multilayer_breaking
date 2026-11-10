/**
# Field scale wave breaking (multilayer solver)
*/

// #include "grid/multigrid.h"
// #include "view.h"
#include "layered/hydro.h"
//#include "layered/hydro_test.h"
#include "layered/nh.h"
#include "layered/remap.h"
//#include "remap_test.h"
#include "layered/perfs.h"
#include "input.h" // Used for older input method (matrices of quantities)
// Updated 2024/12/20: the newest Basilisk output_matrix is already compatible with MPI
#include "output_mpi.h" // Antoon's function for MPI compatible matrix output

/**
Definition of some controlling parameters. */

#define g_ 9.8
double kp_ = 2.*pi/10.; // Peak wave number, not used for initialization
                        // just for scaling the water depth
double h_ = 10; // depth of the water
double gpe_base = 0; // gauge of potential energy
/* double ETAE = 0.1; // refinement criteria for eta */
/* int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function */
/* int MINLEVEL = 6; // min level */
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;

#include "./spectrum.h" // Used for new input method (the spectrum info)

/** Function for writing fields at time t.
*/
int writefields (double t, const char *suffix) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50];
  vector u_temp;
  scalar w_temp, h_temp;
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_%s_t%g_l%d", suffix, t, j);
    sprintf (filename2, "field/uy_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename3, "field/uz_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename4, "field/h_%s_t%g_l%d", suffix, t, j);  
    if (j==0) {
      // The first layer is named u instead of u0
      sprintf (s, "u");
      u_temp = lookup_vector (s);
      sprintf (s, "w");
      w_temp = lookup_field (s);
      sprintf (s, "h");
      h_temp = lookup_field (s);
    }
    else {
      sprintf (s, "u%d", j);
      u_temp = lookup_vector (s);
      sprintf (s, "w%d", j);
      w_temp = lookup_field (s);
      sprintf (s, "h%d", j);
      h_temp = lookup_field (s);
    }
    FILE * fux = fopen (filename1, "w");
    output_matrix_mpi (u_temp.x, fux, N, linear = true);
    fclose (fux);
    FILE * fuy = fopen (filename2, "w");
    output_matrix_mpi (u_temp.y, fuy, N, linear = true);
    fclose (fuy);
    FILE * fuz = fopen (filename3, "w");
    output_matrix_mpi (w_temp, fuz, N, linear = true);
    fclose (fuz);
    FILE * fh = fopen (filename4, "w");
    output_matrix_mpi (h_temp, fh, N, linear = true);
    fclose (fh);    
  }
  return 0;
}

/* Extractor: restore a dump written by field_PM and write surface + full-depth
   matrices (same formats as the movie/writefields events).
   Usage: ./extract_surface DUMPFILE [NLAYER] [LEVEL] [L0] [Lp] */
char * dumpfile = "restart";

int main(int argc, char * argv[])
{
  if (argc > 1)
    dumpfile = argv[1];
  if (argc > 2)
    NLAYER = atoi(argv[2]); // # of layers
  if (argc > 3)
    LEVEL_data = atoi(argv[3]); // Horizontal resolution
  if (argc > 4)
    L0 = atof(argv[4]); // Box size
  if (argc > 5)
    kp_ = 2.*pi/atof(argv[5]); // Peak wavelength
  else
    kp_ = 2.*pi/(L0/5.);
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data;
  nl = NLAYER;
  G = g_;
  h_ = 2.*pi/kp_; // water depth = peak wavelength, must match the run
#if dimension == 2
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#else
  gpe_base = -0.5*sq(h_)*L0*g_;
#endif
  CFL_H = 1;
  fprintf (stderr, "Read in parameters!\n");
  run();
}

event init (i = 0)
{
  /* Read t from the dump header directly: restore()'s event catch-up in this
     snapshot leaves t/iter inconsistent, so we set them explicitly. */
  double tdump = 0.;
  FILE * fh = fopen (dumpfile, "r");
  if (fh) {
    if (fread (&tdump, sizeof(double), 1, fh) < 1)
      tdump = 0.;
    fclose (fh);
  }
  if (!restore (dumpfile)) {
    fprintf (stderr, "restore('%s') failed!\n", dumpfile);
    exit (1);
  }
  // t and iter are read-only globals here; we terminate immediately so
  // only the output labels need the correct time (tdump).
  geometric_beta (1./3., true); // same grid mapping as field_PM restart branch
  char s[80];
  char filename1[80], filename2[80], filename3[80];
  vector u_temp;
  sprintf (s, "u%d", nl-1); // top layer velocity
  u_temp = lookup_vector (s);
  sprintf (filename1, "surface/eta_matrix_%g", tdump);
  sprintf (filename2, "surface/ux_matrix_%g", tdump);
  sprintf (filename3, "surface/uy_matrix_%g", tdump);
  FILE * feta = fopen (filename1, "w");
  output_matrix_mpi (eta, feta, N, linear = true);
  fclose (feta);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear = true);
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi (u_temp.y, fuy, N, linear = true);
  fclose (fuy);
  char *suffix = "matrix";
  writefields (tdump, suffix);
  fprintf (stderr, "Extracted surface + full fields at t = %g from %s\n", tdump, dumpfile);
  return 1; // extraction only: end the run now
}

