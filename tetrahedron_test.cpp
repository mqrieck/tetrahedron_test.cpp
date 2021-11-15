
// tetrahedron_test.cpp (by M. Q. Rieck, updated: 11/14/2021)

// Note: This is test code for the results in my "tetrahedron and toroids" paper, and more.

// Note: Recommend redirecting the output to a file, and scrolling through that file.

// Note: You can use three command line integer parameters to specify angle proportion A:B:C.

// Note: This C++ program uses passing-by-reference. It can be easily converted to a C
// program by altering this aspect of function call, and by changing the includes.

#include <cstdio>
#include <cstdlib>
#include <cmath>

#define M 1000                  // how many (alpha, beta, gamma) points (M^3)?
#define N 80                    // how fine to subdivide the interval [0, pi]
#define O 0                     // set this higher to avoid low "tilt planes"
#define pi M_PI                 // pi = 3.141592654..., of course
//#define EXTRA_RULES_1         // some extra tests that could be superfluous
//#define EXTRA_RULES_2         // some additional such tests
#define ACUTE_TESTING           // only appropriate for an acute base triangle ABC
//#define MAX_RULES             // some testing based on toroid analysis
#define EASY_COSINE_RULES       // more testing based of toroid analysis
#define GRUNERT_DISCR_RULE      // a test based on Grunert's system discriminant
#define REFINED                 // more refined testing for cell acceptance/rejection
//#define SHOW_CUTOFFS          // show when alpha = A, beta = B or gamma = C

using namespace std;

// Obtain "view angles" at P based on triangle ABC and three "tilt angles", the tau's.
// Each tilt angle is the dihedral angle between the ABC side and another side of the
// tetrahedron ABCP. Dihedral angle formulas are used to find the view angles, alpha,
// beta and gamma, at the point P.
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

inline int ind(double angle) {
  int i = (int) (N*angle/pi);
  if (i < 0) i = 0;
  if (i >= N) i = N-1;
  return i;
}

void show_array(int a[N][N][N], int i0, int j0, int k0) {
  printf("\n");
  for (int i=0; i<N-2; i++) {
    printf("α = %1.3f:\n", (i+.5)*pi/N);
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
#ifdef SHOW_CUTOFFS
        if (i==i0 || j==j0 || k==k0) printf("#"); else  
#endif
        switch (a[i][j][k]) {
          case 0:  printf("."); break;  // a "prohibited" cell that is empty
          case 1:  printf("x"); break;  // a "prohibited" cell containing a data pt.
          case 2:  printf(" "); break;  // an "allowable" cell that is empty
          case 3:  printf("o"); break;  // an "allowable" cell containing a data pt.
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

int main(int argc, char **argv) {
  int states[N][N][N], state, total, count0, count1, count2, count3, rejected = 0, i0, j0, k0, i, j, k, x, y,
    choice, delta_i, delta_j, delta_k;
  double A, B, C, cosA, cosB, cosC, alpha, beta, gamma, cos_alpha, cos_beta, cos_gamma, den, tol = 0.05;
  double x10, x20, x30, y10, y20, y30, x1, x2, x3, y1, y2, y3, cos_turn, sin_turn, c1, c2, c3, C0, C1, C2, C3,
    eta_sq, L, R, E, D;
  bool accept;
  // Set the angles for the base triangle ABC
  // Can use three command line integer parameters to specify the proportion A : B : C
  if (argc == 4) {
    den = atoi(argv[1]) + atoi(argv[2]) + atoi(argv[3]);
    A = atoi(argv[1])*pi / den;
    B = atoi(argv[2])*pi / den;
    C = atoi(argv[3])*pi / den;
  } else A = B = C = pi/3;
  cosA = cos(A); cosB = cos(B); cosC = cos(C);
#ifdef SHOW_CUTOFFS
  i0 = ind(A); j0 = ind(B); k0 = ind(C);
  if (i0 < 0) i0 = 0; if (i0 >= N) i0 = N-1;
  if (j0 < 0) j0 = 0; if (j0 >= N) j0 = N-1;
  if (k0 < 0) k0 = 0; if (k0 >= N) k0 = N-1;
#endif
// Assume control point A = (x1, y1) = (1, 0) intially
  x10 = 1; y10 = 0;
  x20 = cos(2*C); y20 =  sin(2*C);
  x30 = cos(2*B); y30 = -sin(2*B);
// But now turn all control points to achieve my standard orientation
  cos_turn = cos(2*(B-C)/3);
  sin_turn = sin(2*(B-C)/3);
  x1 = x10 * cos_turn - y10 * sin_turn;
  y1 = x10 * sin_turn + y10 * cos_turn;
  x2 = x20 * cos_turn - y20 * sin_turn;
  y2 = x20 * sin_turn + y20 * cos_turn;
  x3 = x30 * cos_turn - y30 * sin_turn;
  y3 = x30 * sin_turn + y30 * cos_turn;
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
#ifdef REFINED
        accept = false;
        for (delta_i=0; delta_i<2; delta_i++)
          for (delta_j=0; delta_j<2; delta_j++)
            for (delta_k=0; delta_k<2; delta_k++) {
              alpha = (i+delta_i)*pi/N;
              beta  = (j+delta_j)*pi/N;
              gamma = (k+delta_k)*pi/N;
#else
              alpha = (i+0.5)*pi/N;
              beta  = (j+0.5)*pi/N;
              gamma = (k+0.5)*pi/N;
#endif
        cos_alpha = cos(alpha);
        cos_beta  = cos(beta);
        cos_gamma = cos(gamma);
        // The following is taken from my "Grunert" paper, for the discriminant D
        c1 = cos_alpha; c2 = cos_beta; c3 = cos_gamma;
        C0 = c1*c2*c3; C1 = c1*c1; C2 = c2*c2; C3 = c3*c3;
        eta_sq = 1 - C1 - C2 - C3 + 2*C0;
        L = 2 * ( (1-x1)*y1*(C1-1) + (1-x2)*y2*(C2-1) + (1-x3)*y3*(C3-1) + (  y1+y2+y3)*(1-C0) );
        R = 2 * ( (1+x1)*x1*(C1-1) + (1+x2)*x2*(C2-1) + (1+x3)*x3*(C3-1) + (1+x1+x2+x3)*(1-C0) );
        E = L*L + (R+eta_sq)*(R+eta_sq);
        D = E*E + 18*E*eta_sq*eta_sq + 8*(R+eta_sq)*((R+eta_sq)*(R+eta_sq)-3*L*L)*eta_sq
          - 27*eta_sq*eta_sq*eta_sq*eta_sq;
        if (
          alpha +  beta + gamma < 2*pi &&
          alpha <  beta + gamma &&
          beta  < gamma + alpha &&
          gamma < alpha +  beta 
#ifdef EXTRA_RULES_1
          &&
          A + beta + gamma  < 2*pi &&
          alpha + B + gamma < 2*pi &&
          alpha + beta + C  < 2*pi
#endif
#ifdef EXTRA_RULES_2
          &&
          beta  + gamma - alpha < 2*(B+C) &&
          gamma + alpha -  beta < 2*(C+A) &&
          alpha +  beta - gamma < 2*(A+B)
#endif
#ifdef ACUTE_TESTING
#ifdef MAX_RULES
          &&
          (alpha >= A || beta  < B || beta  < C + alpha) &&
          (alpha >= A || gamma < C || gamma < B + alpha) &&
          (beta  >= B || gamma < C || gamma < A + beta ) &&
          (beta  >= B || alpha < A || alpha < C + beta ) &&
          (gamma >= C || alpha < A || alpha < B + gamma) &&
          (gamma >= C || beta  < B || beta  < A + gamma)
#endif
#ifdef EASY_COSINE_RULES
          &&
          (alpha >= A || cosC * cos_beta  + cosB * cos_gamma > 0) &&
          (beta  >= B || cosA * cos_gamma + cosC * cos_alpha > 0) &&
          (gamma >= C || cosB * cos_alpha + cosA * cos_beta  > 0)
#endif
#ifdef GRUNERT_DISCR_RULE
          && // if outside CSDC then cannot be inside exactly two basic toroids!
          ( D < 0 || (
            (alpha >= A || beta <  B || gamma <  C) &&
            (alpha <  A || beta >= B || gamma <  C) &&
            (alpha <  A || beta <  B || gamma >= C) ) )
#endif
#endif
#ifdef REFINED
              ) { accept = true; break; }
        }
        if (accept) states[i][j][k] += 2;
#else
              ) states[i][j][k] += 2;
#endif
  }
  // Show slices of the array, indicating the nature of each cell.
  show_array(states, i0, j0, k0);
  // Compute and display statistices for the given triangle ABC.
  total = count0 = count1 = count2 = count3 = 0;
  for (int i=0; i<N-2; i++)
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
