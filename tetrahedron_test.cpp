
// tetrahedron_test.cpp (by M. Q. Rieck, updated: 10/17/2021)

// Note: This is test code for the results in my "tetrahedron and toroids" paper.

// Note: Recommend redirecting the output to a file, and scrolling through that file.

// Note: Can use three command line integer parameters to specify proportion A:B:C.

// Note: This C++ program uses passing-by-reference. It can be easily converted to a C
// program by altering this aspect of function call, and by changing the includes.

#include <cstdio>
#include <cmath>

#define M 1100                // how many (alpha, beta, gamma) points (M^3)?
#define N 100                 // how fine to subdivide the interval [0, pi]
#define O 0                   // set higher to avoid low "tilt planes"
#define pi M_PI               // pi = 3.141592654..., of course
#define ACUTE_TEST            // only appropriate for acute base triangle ABC
//#define COSINES_TEST        // include the "cosines test" when using an acute triangle
//#define SHOW_CUTOFFS        // show when alpha = A, beta = B or gamma = C

// The tau's are "tilt angles" for three planes, each containing one of the sidelines of
// the triangle ABC. Dihedral angle formulas are used to find the "view angles", alpha,
// beta and gamma, at the point of intersection of the three planes.
bool tilt_to_view_angles(double tau1, double tau2, double tau3, double cosA, double cosB,
  double cosC, double& alpha, double& beta, double& gamma, int& rejected) {
    double cos_tau1, cos_tau2, cos_tau3, sin_tau1, sin_tau2, sin_tau3;
    double cos_delta1, cos_delta2, cos_delta3, sin_delta1, sin_delta2, sin_delta3;
    cos_tau1 = cos(tau1), cos_tau2 = cos(tau2), cos_tau3 = cos(tau3);
    sin_tau1 = sin(tau1), sin_tau2 = sin(tau2), sin_tau3 = sin(tau3);
    cos_delta1 = sin_tau2 * sin_tau3 * cosA - cos_tau2 * cos_tau3;
    cos_delta2 = sin_tau3 * sin_tau1 * cosB - cos_tau3 * cos_tau1;
    cos_delta3 = sin_tau1 * sin_tau2 * cosC - cos_tau1 * cos_tau2;
    sin_delta1 = sqrt(1 - cos_delta1*cos_delta1);
    sin_delta2 = sqrt(1 - cos_delta2*cos_delta2);
    sin_delta3 = sqrt(1 - cos_delta3*cos_delta3);
    alpha = acos((cos_delta1 + cos_delta2 * cos_delta3) / (sin_delta2 * sin_delta3));
    beta  = acos((cos_delta2 + cos_delta3 * cos_delta1) / (sin_delta3 * sin_delta1));
    gamma = acos((cos_delta3 + cos_delta1 * cos_delta2) / (sin_delta1 * sin_delta2));
    if (alpha < 0 || alpha > pi || beta < 0 || beta > pi || gamma < 0 || gamma > pi ||
      alpha > beta+gamma || beta > gamma+alpha || gamma > alpha+beta || alpha+beta+
        gamma > 2*pi) { rejected++; return false; } else return true;
}

void clear_array(int a[N][N][N]) {
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
        a[i][j][k] = 0;
}

inline int ind(double angle) {
  int i = (int) (N*angle/pi);
  if (i < 0) i = 0;
  if (i >= N) i = N-1;
  return i;
}

void show_array(int a[N][N][N], int i0, int j0, int k0) {
  printf("\n");
  for (int i=0; i<N-1; i++) {
    printf("α = %1.3f:\n", (i+.5)*pi/N);
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
#ifdef SHOW_CUTOFFS
  if (i==i0 || j==j0 || k==k0) printf("#");
  else  switch (a[i][j][k]) {
          case 0:  printf("."); break;  // a "prohibited" cell that is empty
          case 1:  printf("x"); break;  // a "prohibited" cell containing a data pt.
          case 2:  printf(" "); break;  // an "allowable" cell that is empty
          case 3:  printf("o"); break;  // an "allowable" cell containing a data pt.
        }
#else
        switch (a[i][j][k]) {
          case 0:  printf("."); break;  // a "prohibited" cell that is empty
          case 1:  printf("x"); break;  // a "prohibited" cell containing a data pt.
          case 2:  printf(" "); break;  // an "allowable" cell that is empty
          case 3:  printf("o"); break;  // an "allowable" cell containing a data pt.
        }
#endif
      }
      printf("\n");
    }
    printf("\n");
    for (int k=0; k<N; k++) printf("_");
    printf("\n\n");
  }
  printf("\n");
}

int main(int argc, char **argv) {

  int states[N][N][N], state, total, count0, count1, count2, count3, rejected = 0, i0, j0, k0;
  double A, B, C, cosA, cosB, cosC, alpha, beta, gamma, cos_alpha, cos_beta, cos_gamma, den;
  // Set the angles for the base triangle ABC
  // Can use three command line integer parameters to specify the proportion A : B : C
  if (argc == 4) {
    den = atoi(argv[1]) + atoi(argv[2]) + atoi(argv[3]);
    A = atoi(argv[1])*pi / den;
    B = atoi(argv[2])*pi / den;
    C = atoi(argv[3])*pi / den;
  } else {
  // Or else specify the proportions here
  //   Acute triangles:
    A =  8;  B =  6;  C =  5;
//  A =  4;  B =  6;  C =  9;
//  A =  5;  B =  7;  C =  7;
//  A =  9;  B =  9;  C =  1;
//  A =  2;  B =  9;  C =  8;
//  A =  5;  B =  9;  C =  5;
  //  Obtuse triangles:
//  A = 12;  B =  3;  C =  4;
//  A =  4;  B =  3;  C = 12;
//  A =  5;  B = 10;  C =  4;
    den = A + B + C;
    A = A*pi / den; B = B*pi / den; C = C*pi / den;
  }
  cosA = cos(A); cosB = cos(B); cosC = cos(C);
#ifdef SHOW_CUTOFFS
  i0 = ind(A); j0 = ind(B); k0 = ind(C);
  if (i0 < 0) i0 = 0; if (i0 >= N) i0 = N-1;
  if (j0 < 0) j0 = 0; if (j0 >= N) j0 = N-1;
  if (k0 < 0) k0 = 0; if (k0 >= N) k0 = N-1;
#endif

//  clear_array(states);
  printf("\n\nThe base triangle angles: A = %.4f , B = %.4f , C = %.4f\n\n", A, B, C);
  printf("The following plots show slices of the cube [0,π] x [0,π] x [0,π] whose coordinates are α, β and γ. A system\n");
  printf("of inequalities defines an \"allowable\" portion of this cube. The slices are divided into cells. Each cell is\n");
  printf("designated to be \"allowable\" or \"unallowable,\" based on the system of inequalities. However, a cell that has\n");
  printf("been designated to be \"unallowable\" might actually contain some of the allowable portion of the cube together\n");
  printf("with some of the unallowable portion of the cube, in which case calling the cell \"unallowable\" is an unfortunate\n");
  printf("mistake. This can only happen at the boundary of the allowable portion of the cube.\n\n"); 

  printf("The allowable portion of the cube bounds all of the points (α, β, γ) for which α, β and γ can be the angles at\n");
  printf("a point P = (x, y, z) that extends the triangle ABC to form a tetrahedron ABCP. If a cell contains such a point\n");
  printf("(α, β, γ), then we call it \"occupied;\" otherwise the cell is \"unoccupied.\" (A basic understanding of the problem\n");
  printf("in the paper is presumed here.)\n\n");

  printf("Each cell is represented by a character. A space character represents an unoccupied allowable cell, an \'o\'\n");
  printf("represents an occupied allowable cell, a dot represents an unoccupied unallowable cell, and an \'x\' represents\n");
  printf("an occupied unallowable cell. This latter case is possible since an \"unallowable\" cell might contain an allowable\n");
  printf("portion of the cube (when it contains part of the boundary). Pound signs show where α = A, β = B or γ = C.\n");

  // Use 3D array to record possible (alpha, beta, gamma) triples for given triangle
  for (int i=O; i<M-O; i++)
    for (int j=O; j<M-O; j++)
      for (int k=O; k<M-O; k++)
        if (tilt_to_view_angles(i*pi/M, j*pi/M, k*pi/M, cosA, cosB, cosC, alpha,
          beta, gamma, rejected)) states[ind(alpha)][ind(beta)][ind(gamma)] = 1;

  // Also use array to record which cells in the array are within system of bounds
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++) {
        alpha = (i+.5)*pi/N;
        beta  = (j+.5)*pi/N;
        gamma = (k+.5)*pi/N;
        cos_alpha = cos(alpha);
        cos_beta  = cos(beta);
        cos_gamma = cos(gamma);
        if (
          alpha +  beta + gamma < 2*pi &&
          alpha <  beta + gamma &&
           beta < gamma + alpha &&
          gamma < alpha +  beta &&
          A + beta + gamma  < 2*pi &&
          alpha + B + gamma < 2*pi &&
          alpha + beta + C  < 2*pi &&
           beta + gamma - alpha < 2*(B+C) &&
          gamma + alpha -  beta < 2*(C+A) &&
          alpha +  beta - gamma < 2*(A+B)
#ifdef ACUTE_TEST
&&
          (alpha >= A || beta  < B || beta  < C + alpha) &&
          (alpha >= A || gamma < C || gamma < B + alpha) &&
          (beta  >= B || gamma < C || gamma < A + beta ) &&
          (beta  >= B || alpha < A || alpha < C + beta ) &&
          (gamma >= C || alpha < A || alpha < B + gamma) &&
          (gamma >= C || beta  < B || beta  < A + gamma)
#ifdef COSINES_TEST
&&
          (alpha >= A || cosC * cos_beta  + cosB * cos_gamma > 0) &&
          (beta  >= B || cosA * cos_gamma + cosC * cos_alpha > 0) &&
          (gamma >= C || cosB * cos_alpha + cosA * cos_beta  > 0)
#endif
#endif
        ) states[i][j][k] += 2;
  }

  // Show slices of the array, indicating the nature of each cell.
  show_array(states, i0, j0, k0);

  // Compute and display statistices for the given triangle ABC.
  total = count0 = count1 = count2 = count3 = 0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++) {
        switch (states[i][j][k]) {
          case 0: count0++; break;
          case 1: count1++; break;
          case 2: count2++; break;
          case 3: count3++;
        }
        total++;
  }
  printf("Number of   occupied   allowable cells:    %8d\n",   count3);
  printf("Number of unoccupied   allowable cells:    %8d\n",   count2);
  printf("Number of   occupied unallowable cells:    %8d\n",   count1);
  printf("Number of unoccupied unallowable cells:    %8d\n",   count0);
  printf("Total number of cells in the array:        %8d\n",   total);
  printf("Number of rejected calls for a data point: %8d\n\n", rejected);
}
