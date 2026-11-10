/**
# Field scale wave breaking with wind forcing (multilayer solver)

Same as field_PM.c but with a wind stress applied on the top layer:
aligned or misaligned wind-wave cases at high wind speeds (U10 = 30-60
m/s). Two extra command line arguments:

  ... NLAYER LEVEL TEND nu RAND L0 Lp Htheta U10 THETA_W

with U10 the wind speed at 10 m (m/s, 0 = no wind, recovers the
decaying case) and THETA_W the wind direction relative to the dominant
wave direction +x (degrees).
*/

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "input.h"
#include "output_mpi.h"

/**
Definition of some controlling parameters. */

#define g_ 9.8
double kp_ = 2.*pi/10.;
double h_ = 10;
double gpe_base = 0;
double TEND = 50.;
int NLAYER = 10;
int LEVEL_data = 7;

/**
Wind parameters. The surface stress uses the relative-wind form

$$\vec\tau = \frac{\rho_a}{\rho_w} C_d \|\vec U_{10} - \vec u_{top}\|
(\vec U_{10} - \vec u_{top})$$

with the Large & Pond (1981) drag coefficient, capped at CD_MAX to
account for the saturation observed at hurricane wind speeds (Powell
et al. 2003; Donelan et al. 2004). */

double U10 = 0.;          // wind speed at 10 m (m/s), 0 = no wind
double THETA_W = 0.;      // wind direction relative to +x (degrees)
double CD_MAX = 2.5e-3;   // saturation value of the drag coefficient

/**
## Typhoon-passage mode (tier-1 realism)

WIND_MODE 0: constant U10/THETA_W (the matrix experiments).
WIND_MODE 1: the patch sits at a fixed point while a translating
Holland (1980) vortex passes it.  The storm moves at VT along
direction TRK; its center passes at closest-approach distance DIMP.
Wind at the patch = gradient tangential wind + VTFRAC*VT translation
asymmetry.  Wind time = t + T0, so T0 < 0 starts the run before
closest approach.  All parameters via argv 11..19. */

double WIND_MODE = 0.;
double VT   = 8.;       // storm translation speed (m/s)
double RM   = 40.e3;    // radius of maximum wind (m)
double BH   = 1.5;      // Holland shape parameter
double DPC  = 60.e2;    // central pressure deficit (Pa)
double CORF = 5.e-5;    // Coriolis parameter (s^-1)
double DIMP = 20.e3;    // closest-approach distance of the track (m)
double TRK  = 90.;      // storm translation direction (deg from +x)
double T0W  = 0.;       // time offset: wind time = t + T0W
double VTFRAC = 0.7;    // translation asymmetry fraction
#define RHOA 1.15       // air density (kg/m^3)
#define RHO_RATIO 1.25e-3 // rho_air / rho_water

#include "./spectrum.h"

/** Wind stress given the top-layer velocity. */
/** Instantaneous wind vector sampled from the passing vortex (mode 1)
or the constant settings (mode 0). */

static void wind_now (double time, double * u10, double * thdeg) {
  if (WIND_MODE < 0.5) { *u10 = U10; *thdeg = THETA_W; return; }
  double ts = time + T0W;                    // storm-relative time
  double tr = TRK*pi/180.;
  double ax = cos(tr), ay = sin(tr);         // along-track unit vector
  double px = -ay, py = ax;                  // +90 deg of it
  double xc = ts*VT*ax + DIMP*px;            // storm center relative to patch
  double yc = ts*VT*ay + DIMP*py;
  double r = sqrt (sq(xc) + sq(yc)) + 1.;
  double rp = pow (RM/r, BH);
  double dpdr = DPC*exp(-rp)*rp*BH/r;        // Holland radial pressure gradient
  double vgr = -CORF*r/2. + sqrt (sq(CORF*r/2.) + r/RHOA*dpdr);
  double wx = vgr*(-yc/r) + VTFRAC*VT*ax;    // CCW tangential + translation
  double wy = vgr*( xc/r) + VTFRAC*VT*ay;
  *u10 = sqrt (sq(wx) + sq(wy));
  *thdeg = atan2 (wy, wx)*180./pi;
}

static void wind_tau (double ux, double uy, double * tx, double * ty) {
  double u10, thd;
  wind_now (t, &u10, &thd);
  double cd = (0.49 + 0.065*u10)*1.e-3;
  if (cd > CD_MAX) cd = CD_MAX;
  double th = thd*pi/180.;
  double ur = u10*cos(th) - ux, vr = u10*sin(th) - uy;
  double um = sqrt (sq(ur) + sq(vr));
  *tx = RHO_RATIO*cd*um*ur;
  *ty = RHO_RATIO*cd*um*vr;
}

/** Function for writing fields at time t. */
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

int main(int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]);
  if (argc > 2)
    LEVEL_data = atoi(argv[2]);
  if (argc > 3)
    TEND = atof(argv[3]);
  if (argc > 4)
    nu = atof(argv[4]);
  else
    nu = 0.;
  if (argc > 5)
    RANDOM = atoi(argv[5]);
  if (argc > 6)
    L0 = atof(argv[6]);
  else
    L0 = 50.;
  if (argc > 7)
    kp_ = 2.*pi/atof(argv[7]);
  else
    kp_ = 2.*pi/(L0/5.);
  if (argc > 8)
    theta_H = atof(argv[8]);
  else
    theta_H = 0.5;
  if (argc > 9)
    U10 = atof(argv[9]);      // wind speed at 10 m (0 = no wind)
  if (argc > 10)
    THETA_W = atof(argv[10]); // wind direction (degrees from +x)
  // Typhoon-passage mode: WIND_MODE VT RM_km BH DPC_hPa CORF_1e-5 DIMP_km TRK_deg T0
  if (argc > 11) {
    WIND_MODE = atof(argv[11]);
    if (argc > 12) VT = atof(argv[12]);
    if (argc > 13) RM = atof(argv[13])*1.e3;
    if (argc > 14) BH = atof(argv[14]);
    if (argc > 15) DPC = atof(argv[15])*1.e2;
    if (argc > 16) CORF = atof(argv[16])*1.e-5;
    if (argc > 17) DIMP = atof(argv[17])*1.e3;
    if (argc > 18) TRK = atof(argv[18]);
    if (argc > 19) T0W = atof(argv[19]);
  }
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data;
  nl = NLAYER;
  G = g_;
  h_ = 2.*pi/kp_;
#if dimension == 2
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#else
  gpe_base = -0.5*sq(h_)*L0*g_;
#endif
  CFL_H = 1;
  fprintf (stderr, "Read in parameters! U10 = %g, THETA_W = %g, WIND_MODE = %g\n",
           U10, THETA_W, WIND_MODE);
  if (WIND_MODE > 0.5)
    fprintf (stderr, "Holland passage: VT=%g RM=%gkm BH=%g DPC=%ghPa f=%g "
             "DIMP=%gkm TRK=%g T0=%g\n",
             VT, RM/1.e3, BH, DPC/1.e2, CORF, DIMP/1.e3, TRK, T0W);
  run();
}

/**
## Initialization

Identical to field_PM.c (spectrum input or restore). */

event init (i = 0)
{
  if (!restore ("restart")) {
    power_input();
    dkx_ = kx_[1] - kx_[0];
    dky_ = ky_[1] - ky_[0];
    fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_);
    geometric_beta (1./3., true);
    foreach() {
      zb[] = -h_;
      eta[] = wave(x, y);
      double H = wave(x, y) - zb[];
      foreach_layer() {
	      h[] = H/nl;
      }
    }
    vertical_remapping (h, tracers);
    foreach() {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        u.x[] = u_x(x, y, z);
        u.y[] = u_y(x, y, z);
        w[] = u_z(x, y, z);
        z += h[]/2.;
      }
    }
    fprintf (stderr,"Done initialization!\n");
    dump("initial");
    char *suffix = "matrix";
    writefields (t, suffix);
  }
  else {
    geometric_beta (1./3., true);
    dtmax = 0.01;
    dt = dtnext (dtmax);
    char *suffix = "matrix";
    writefields (t, suffix);
  }
}

/** Horizontal diffusion. */

event viscous_term (i++)
  horizontal_diffusion ((scalar *){u}, nu, dt);

/**
## Wind stress on the top layer

Added to the face acceleration field ha between face_fields and
pressure (the standard multilayer acceleration hook, as in
coriolis.h). Since ha is d(hu)/dt, adding tau gives exactly the
momentum-flux boundary condition: du_top = tau*dt/h_top, applied
consistently to both fluxes and velocities. */

event acceleration (i++, last)
{
  if (U10 > 0. || WIND_MODE > 0.5) {
    // x-faces: normal = x, transverse = y
    foreach_face (x)
      foreach_layer()
        if (point.l == nl - 1) {
          double hl = h[-1] > dry ? h[-1] : 0.;
          double hr = h[]   > dry ? h[]   : 0.;
          if (hl + hr > dry) {
            double ux = (hl*u.x[-1] + hr*u.x[])/(hl + hr);
            double uy = 0.5*(u.y[-1] + u.y[]);
            double tx, ty;
            wind_tau (ux, uy, &tx, &ty);
            ha.x[] += tx;
          }
        }
    // y-faces: normal = y, transverse = x
    foreach_face (y)
      foreach_layer()
        if (point.l == nl - 1) {
          double hl = h[-1] > dry ? h[-1] : 0.;
          double hr = h[]   > dry ? h[]   : 0.;
          if (hl + hr > dry) {
            double vx = 0.5*(u.x[-1] + u.x[]);
            double vy = (hl*u.y[-1] + hr*u.y[])/(hl + hr);
            double tx, ty;
            wind_tau (vx, vy, &tx, &ty);
            ha.y[] += ty;
          }
        }
  }
}

/**
## Outputs */

event energy_before_remap (i+=10, last)
{
  if (i==10) {
    fprintf(stderr, "energy output before remap!\n");
    fflush(stderr);
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer () {
      double norm2 = sq(w[]);
      foreach_dimension()
        norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  static FILE * fp = fopen("energy_before_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe - gpe_base);
  fflush (fp);
}

/** Wind power input per unit area, P_in = int tau.u_top dA / L0^2. */

event energy_wind (i+=10, last)
{
  if (U10 > 0. || WIND_MODE > 0.5) {
    double pin = 0.;
    foreach (reduction(+:pin)) {
      foreach_layer()
        if (point.l == nl - 1) {
          double tx, ty;
          wind_tau (u.x[], u.y[], &tx, &ty);
          pin += (tx*u.x[] + ty*u.y[])*dv();
        }
    }
    static FILE * fp = fopen("wind_input.dat","w");
    if (pid() == 0) {
      fprintf (fp, "%g %g\n", t, pin/sq(L0));
      fflush (fp);
    }
  }
}

/** Log the imposed 10-m wind vector (mode 1: the passing typhoon). */

event wind_log (i += 10, last)
{
  if (WIND_MODE > 0.5 && pid() == 0) {
    double u10, thd;
    wind_now (t, &u10, &thd);
    static FILE * fp = fopen("wind10_time.dat","w");
    fprintf (fp, "%g %g %g\n", t, u10, thd);
    fflush (fp);
  }
}

/**
## Momentum budget

Total water-column momentum split into the depth-mean current and the
wave/orbital part (total minus column depth-mean), along x and y.
Together with the imposed wind stress (whose integral input rate is
exactly d(M)/dt since momentum is only exchanged at the surface), this
gives the partition of the momentum input between current growth and
wave-field growth (radiation stress storage) - the diagnostic for how
much along-wind stress the water column actually feels under wind-wave
misalignment. */

event momentum_budget (i += 10, last)
{
  scalar ubx[], uby[];
  foreach() {
    double H = 0., sx = 0., sy = 0.;
    foreach_layer() { H += h[]; sx += h[]*u.x[]; sy += h[]*u.y[]; }
    if (H > dry) { ubx[] = sx/H; uby[] = sy/H; }
    else { ubx[] = 0.; uby[] = 0.; }
  }
  double mxm = 0., mxw = 0., mym = 0., myw = 0.;
  foreach (reduction(+:mxm) reduction(+:mxw)
	   reduction(+:mym) reduction(+:myw)) {
    foreach_layer() {
      double dvh = h[]*dv();
      mxm += dvh*ubx[];
      mym += dvh*uby[];
      mxw += dvh*(u.x[] - ubx[]);
      myw += dvh*(u.y[] - uby[]);
    }
  }
  static FILE * fp = fopen ("momentum.dat", "w");
  if (pid() == 0) {
    fprintf (fp, "%g %g %g %g %g\n", t, mxm, mxw, mym, myw);
    fflush (fp);
  }
}

/** Horizontally-averaged velocity and thickness of each layer (drift
and shear profiles), every t = 10. */

event drift_profile (t += 10; t <= TEND)
{
  FILE * fp = pid() == 0 ? fopen ("drift_profile.dat", "a") : NULL;
  for (int l = 0; l < nl; l++) {
    double sx = 0., sy = 0., sh = 0., a = 0.;
    foreach (reduction(+:sx) reduction(+:sy) reduction(+:sh) reduction(+:a)) {
      foreach_layer()
        if (point.l == l) {
          sx += u.x[]*dv(); sy += u.y[]*dv(); sh += h[]*dv(); a += dv();
        }
    }
    if (fp)
      fprintf (fp, "%g %d %g %g %g\n", t, l, sx/a, sy/a, sh/a);
  }
  if (fp)
    fclose (fp);
}

/**
Per-layer second moments of the horizontal velocity (deviation from the
layer mean), every t = 1: <u>, <v>, <u'u'>, <v'v'>, <u'v'>, <h>.
Together with the energy budget this gives the vertical structure of
the wave + coherent-turbulence kinetic energy injected by breaking
(wave orbital energy decays like e^{kz}; near-surface excess above
that baseline is the breaking injection), the quantity to be
parameterized as a surface TKE flux / dissipation profile in GCM
turbulence closures. */

event layer_variance (t += 1; t <= TEND)
{
  FILE * fp = pid() == 0 ? fopen ("layer_variance.dat", "a") : NULL;
  for (int l = 0; l < nl; l++) {
    double sx = 0., sy = 0., sxx = 0., syy = 0., sxy = 0., a = 0., sh = 0.;
    foreach (reduction(+:sx) reduction(+:sy) reduction(+:sxx)
	     reduction(+:syy) reduction(+:sxy) reduction(+:a)
	     reduction(+:sh)) {
      foreach_layer()
        if (point.l == l) {
          double dva = dv();
          sx += u.x[]*dva;
          sy += u.y[]*dva;
          sxx += u.x[]*u.x[]*dva;
          syy += u.y[]*u.y[]*dva;
          sxy += u.x[]*u.y[]*dva;
          a += dva;
          sh += h[]*dva;
        }
    }
    double ux = sx/a, uy = sy/a;
    if (fp)
      fprintf (fp, "%g %d %g %g %g %g %g %g\n", t, l,
               ux, uy, sxx/a - ux*ux, syy/a - uy*uy, sxy/a - ux*uy, sh/a);
  }
  if (fp)
    fclose (fp);
}

/**
Surface matrices (eta, ux, uy) every t = 1 for postprocessing. */

event movie (t += 1; t <= TEND)
{
  char s[80];
  char filename1[50], filename2[50], filename3[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  sprintf (filename3, "surface/uy_matrix_%g", t);
  FILE * feta = fopen (filename1, "w");
  output_matrix_mpi (eta, feta, N, linear=true);
  fclose (feta);
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear=true);
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi (u_temp.y, fuy, N, linear=true);
  fclose (fuy);
}

event regulardump (t = 0; t += 20; t < TEND) {
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}

event endrun (t = TEND) {
  dump();
}
