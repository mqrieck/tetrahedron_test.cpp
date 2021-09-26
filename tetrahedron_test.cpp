
// tetrahedron_test.cpp (by M. Q. Rieck, updated: 9/26/2021)

// Note: This is test code for the results in my "tetrahedron and toroids" paper.

// Note: Recommend redirecting the output to a file, and scrolling through that file.

// Note: This C++ program uses passing-by-reference. It can be easily converted to a C
// program by altering this aspect of function call, and by changing the includes.

#include <cstdio>
#include <cmath>

#define M 1100                // how many (alpha, beta, gamma) points (M^3)?
#define N 100                 // how fine to subdivide the interval [0, pi]
#define O 1                   // set higher to avoid low "tilt planes"
#define pi M_PI               // pi = 3.141592654..., of course
#define ACUTE_TEST            // only appropriate for acute base triangle ABC
//#define COSINES_TEST        // include the "cosines test" when using an acute triangle
#define SHOW_CUTOFFS          // show when alpha = A, beta = B or gamma = C

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

void show_array(int a[N][N][N]) {
  printf("\n\n\n");
  for (int i=0; i<N-1; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
        switch (a[i][j][k]) {
          case 0:  printf("."); break;  // a "prohibited" cell that is empty
          case 1:  printf("x"); break;  // a "prohibited" cell containing a data pt.
          case 2:  printf(" "); break;  // an "allowable" cell that is empty
          case 3:  printf("o"); break;  // an "allowable" cell containing a data pt.
          default: printf("#");
        }
      }
      printf("\n");
    }
    printf("\n");
    for (int k=0; k<N; k++) printf("_");
    printf("\n\n");
  }
  printf("\n");
}

int main() {

  int states[N][N][N], state, total, count0, count1, count2, count3, rejected = 0;
  double A, B, C, cosA, cosB, cosC, alpha, beta, gamma, cos_alpha, cos_beta, cos_gamma;

  //  Set the angles for the base triangle ABC
  //  Acute triangles:
   A =  8*pi/19;  B =  6*pi/19;  C =  5*pi/19;
// A =  4*pi/19;  B =  6*pi/19;  C =  9*pi/19;
// A =  5*pi/19;  B =  7*pi/19;  C =  7*pi/19;
// A =  9*pi/19;  B =  9*pi/19;  C =  1*pi/19;
// A =  2*pi/19;  B =  9*pi/19;  C =  8*pi/19;
// A =  5*pi/19;  B =  9*pi/19;  C =  5*pi/19;
  // Obtuse triangles:
// A =  12*pi/19; B =  3*pi/19;  C =  4*pi/19;
// A =  4*pi/19;  B =  3*pi/19;  C = 12*pi/19;
// A =  5*pi/19;  B = 10*pi/19;  C =  4*pi/19;

  cosA = cos(A); cosB = cos(B); cosC = cos(C);
//  clear_array(states);

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
          (alpha > A || beta  < B || beta  < C + alpha) &&
          (alpha > A || gamma < C || gamma < B + alpha) &&
          (beta  > B || gamma < C || gamma < A + beta ) &&
          (beta  > B || alpha < A || alpha < C + beta ) &&
          (gamma > C || alpha < A || alpha < B + gamma) &&
          (gamma > C || beta  < B || beta  < A + gamma)
#ifdef COSINES_TEST
&&
          (alpha > A || cosC * cos_beta  + cosB * cos_gamma > 0) &&
          (beta  > B || cosA * cos_gamma + cosC * cos_alpha > 0) &&
          (gamma > C || cosB * cos_alpha + cosA * cos_beta  > 0)
#endif
#endif
        ) states[i][j][k] += 2;
#ifdef SHOW_CUTOFFS
        if (fabs(alpha-A) < .015)    states[i][j][k] = 10;
        if (fabs( beta-B) < .015)    states[i][j][k] = 10;
        if (fabs(gamma-C) < .015)    states[i][j][k] = 10;
//      if (fabs(alpha-pi+A) < .015) states[i][j][k] = 10;
//      if (fabs( beta-pi+B) < .015) states[i][j][k] = 10;
//      if (fabs(gamma-pi+C) < .015) states[i][j][k] = 10;
#endif
  }

  // Show slices of the array, indicating the nature of each cell.
  show_array(states);

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
  printf("Number of   occupied   \"allowable\" cells:    %8d\n", count3);
  printf("Number of unoccupied   \"allowable\" cells:    %8d\n", count2);
  printf("Number of   occupied \"unallowable\" cells:    %8d\n", count1);
  printf("Number of unoccupied \"unallowable\" cells:    %8d\n", count0);
  printf("Total number of cells in the array:          %8d\n", total);
  printf("Number of rejected calls for a data point:   %8d\n", rejected);
  printf("(Note: near the boundary, an \"unallowable\" cell might actually ");
  printf("have an \"allowable\" portion.)\n\n");
}
