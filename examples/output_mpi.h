/*
This function is an extension of output_matrix() inside "output.h" which allows
to write 2D fields in a gnuplot-complatible format when running in MPI by
performing a MPI_Reduce.
Updated 2024/12/20: the newest Basilisk output_matrix is already compatible with MPI so this file will not be needed.
*/

// trace
// void output_matrix_mpi (scalar f,
// 		    FILE * fp = stdout,
// 		    int n = N,
// 		    bool linear = false,
// 		    const char * file = NULL,
// 		    coord box[2] = {{X0, Y0}, {X0 + L0, Y0 + L0}})
// {
//   coord cn = {n}, p;
//   double delta = (box[1].x - box[0].x)/n;
//   cn.y = (int)((box[1].y - box[0].y)/delta);
    
//   double ** ppm = (double **) matrix_new (cn.x, cn.y, sizeof(double));
//   double * ppm0 = &ppm[0][0];
//   unsigned int len = cn.x*cn.y;
//   for (int i = 0; i < len; i++)
//     ppm0[i] = - HUGE;

// #if _MPI
//   foreach_region (p, box, cn, reduction(max:ppm0[:len]))
// #else
//   foreach_region (p, box, cn, cpu)
// #endif
//   {
//     int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;
//     int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;
//     double ** alias = ppm; // so that qcc considers ppm a local variable
//     alias[i][j] = linear ? interpolate_linear (point, f, p.x, p.y, p.z) : f[];
//   }
  
//   if (pid() == 0) {
//     if (file) {
//       fp = fopen (file, "wb");
//       if (!fp) {
// 	perror (file);
// 	exit (1);
//       }
//     }
//     float fn = cn.y;
//     fwrite (&fn, sizeof(float), 1, fp);
//     coord delta = {(box[1].x - box[0].x)/cn.x, (box[1].y - box[0].y)/cn.y};
//     for (int j = 0; j < cn.y; j++) {
//       float yp = box[0].y + delta.y*(j + 0.5);
//       fwrite (&yp, sizeof(float), 1, fp);
//     }
//     for (int i = 0; i < cn.x; i++) {
//       float xp = box[0].x + delta.x*(i + 0.5);
//       fwrite (&xp, sizeof(float), 1, fp);
//       for (int j = 0; j < cn.y; j++) {
// 	float z = ppm[i][j];
// 	fwrite (&z, sizeof(float), 1, fp);
//       }
//     }
//     if (file)
//       fclose (fp);
//     else
//       fflush (fp);
//   }
    
//   matrix_free (ppm);
// }

/* Note: the signature below was rewritten into the "struct parameter" style
   used by the (2021-12) Basilisk snapshot provided in basilisk.tar.gz. The
   previous signature used the newer (2023+) default-argument style which that
   qcc does not support. Calls such as
     output_matrix_mpi (f, fp, n, linear = true)
   are unchanged: qcc translates them into a compound literal of the struct. */

struct OutputMatrixMpi {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

void output_matrix_mpi (struct OutputMatrixMpi p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));

  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      field[i][j] = interpolate(p.f, xp, yp, p.linear);
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < p.n; i++){
      float xp = Delta*i + X0 + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
@endif

  matrix_free (field);
}