
// tetrahedron_test.cpp (by M. Q. Rieck, updated: 3/31/2022)

// Note: This is test code for the results in my "tetrahedron and toroids" paper, and beyond.

// Note: Recommend redirecting the output to a file, and scrolling through that file.

// Note: You can use three command line integer parameters to specify angle proportion A:B:C.

// Note: This C++ program uses passing-by-reference. It can be easily converted to a C
// program by altering this aspect of function call, and by changing the includes.

// Note: For faster results, reduce M, N and/or REF_NUM (or comment out REFINED). Using the
// settings M = 700, N = 50 and REF_NUM = 10, and working with the equilateral triangle case,
// it should take a few minutes to produce the following results:
//
//    Number of   occupied   allowable cells:    24648
//    Number of unoccupied   allowable cells:        4
//    Number of   occupied unallowable cells:        2
//    Number of unoccupied unallowable cells:    95346
//
// This means that out of 120000 cells, there were 6 incorrect results. This is an error
// rate of 6 / 120000 = 0.00005. By increasing M and REF_NUM, this number can be reduced.
// Similar results can be obtained for other acute triangles.

#define M 700                   // how many (alpha, beta, gamma) points (M^3)?
#define N 50                    // how fine to subdivide the interval [0, pi]
#define O 0                     // set this higher to avoid low "tilt planes"
#define pi M_PI                 // pi = 3.141592654..., of course
#define TOL1 0.05               // tolerance for some inequalities
#define TOL2 0                  // tolerance for some other inequalities
//#define RESTRICT_DATA         // reject potentially troublesome data points
#define BASIC_RULES_1           // Some basic linear rules
#define BASIC_RULES_2           // Some more basic linear rules
#define BASIC_RULES_3           // Even more basic linear rules
#define ACUTE_TESTING           // only appropriate for an acute base triangle ABC
#define MAX_RULES               // some testing based on toroid analysis
#define EASY_COSINE_RULES       // more testing based of toroid analysis
//#define GRUNERT_DISCR_RULE_1  // a test based on Grunert's system discriminant
#define GRUNERT_DISCR_RULE_2    // a possibly more restictive version of that
#define REFINED                 // more refined testing for cell acceptance/rejection
#define REF_NUM 10              // how much refinement?
//#define SHOW_ARRAY		// display the slices
//#define SHOW_CUTOFFS          // show when alpha = A, beta = B or gamma = C

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ncurses.h>

using namespace std;

// Obtain the "view angles" at P based on triangle ABC and three "tilt angles" (the tau's).
// Each tilt angle is the dihedral angle between the ABC side and another side of the
// tetrahedron ABCP. Dihedral angle formulas are used to find the view angles at P, i.e.,
// the angles alpha = <BPC, beta = <CPA, and gamma = <APB.
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
#ifdef RESTRICT_DATA
    if (abs(sin_delta1) < TOL2 || abs(sin_delta2) < TOL2 || abs(sin_delta3) < TOL2) 
      { rejected++; return false; }
#endif
    alpha = acos((cos_delta1 + cos_delta2 * cos_delta3) / (sin_delta2 * sin_delta3));
    beta  = acos((cos_delta2 + cos_delta3 * cos_delta1) / (sin_delta3 * sin_delta1));
    gamma = acos((cos_delta3 + cos_delta1 * cos_delta2) / (sin_delta1 * sin_delta2));
#ifdef RESTRICT_DATA
    if (alpha < 0 || alpha > pi || beta < 0 || beta > pi || gamma < 0 || gamma > pi ||
      alpha > beta+gamma || beta > gamma+alpha || gamma > alpha+beta || alpha+beta+
        gamma > 2*pi) { rejected++; return false; } 
#endif
    return true;
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

bool test_bounds(double A, double B, double C, double cosA, double cosB, double cosC, double x1, double x2, double x3,
 double y1, double y2, double y3, double alpha, double beta, double gamma) {
  double cos_alpha, cos_beta, cos_gamma, c1, c2, c3, C0, C1, C2, C3, H, L, R, E, D;
  bool accept;
  cos_alpha = cos(alpha);
  cos_beta  = cos(beta);
  cos_gamma = cos(gamma);
  // The following formulas come from my "Grunert" paper and from my paper with Bo Wang
  c1 = cos_alpha; c2 = cos_beta; c3 = cos_gamma;
  C0 = c1*c2*c3; C1 = c1*c1; C2 = c2*c2; C3 = c3*c3;
  H = 1 - C1 - C2 - C3 + 2*C0;
  L = 2 * ( (1-x1)*y1*(C1-1) + (1-x2)*y2*(C2-1) + (1-x3)*y3*(C3-1) + (  y1+y2+y3)*(1-C0) );
  R = 2 * ( (1+x1)*x1*(C1-1) + (1+x2)*x2*(C2-1) + (1+x3)*x3*(C3-1) + (1+x1+x2+x3)*(1-C0) );
  E = L*L + (R+H)*(R+H);
  D = E*E + 18*E*H*H + 8*(R+H)*((R+H)*(R+H)-3*L*L)*H - 27*H*H*H*H;
  accept = (
#ifdef BASIC_RULES_1
    alpha +  beta + gamma < 2*pi &&
    alpha <  beta + gamma &&
    beta  < gamma + alpha &&
    gamma < alpha +  beta
#else
    H > 0
#endif
#ifdef BASIC_RULES_2
    &&
    A + beta + gamma  < 2*pi &&
    alpha + B + gamma < 2*pi &&
    alpha + beta + C  < 2*pi
#endif
#ifdef BASIC_RULES_3
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
#ifdef GRUNERT_DISCR_RULE_1
    && // if outside the CSDC then cannot be inside exactly two of the basic toroids
    ( D < 0 || (
      (alpha >= A || beta <  B || gamma <  C) &&
      (alpha <  A || beta >= B || gamma <  C) &&
      (alpha <  A || beta <  B || gamma >= C) ) )
#endif
#ifdef GRUNERT_DISCR_RULE_2
    && // if outside the CSDC then cannot be inside exactly two of the basic toroids
       // nor inside all three basic toroids and outside their supplementary toroids
    ( D < 0 || ! (
      (alpha <  A && beta >= B && gamma >= C)  ||
      (alpha >= A && beta <  B && gamma >= C)  ||
      (alpha >= A && beta >= B && gamma  < C)  ||
      (alpha >= A && beta >= B && gamma >= C && alpha < pi-A && beta < pi-B && gamma < pi-C) ) )
#endif
#endif
  );
  return accept;
}

int main(int argc, char **argv) {
  int states[N][N][N], state, total, count0, count1, count2, count3, rejected = 0, i0, j0, k0, i, j, k,
    x, y, choice, delta_i, delta_j, delta_k;
  double A, B, C, cosA, cosB, cosC, x1, x2, x3, y1, y2, y3, x10, x20, x30, y10, y20, y30,
    alpha, beta, gamma, cos_alpha, cos_beta, cos_gamma, den, cos_turn, sin_turn;
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
// Set control points (at least intially)
  x10 = 1; y10 = 0;
  x20 = cos(2*C); y20 =  sin(2*C);
  x30 = cos(2*B); y30 = -sin(2*B);
// Turn all control points to achieve my standard orientation
  cos_turn = cos(2*(B-C)/3);
  sin_turn = sin(2*(B-C)/3);
  x1 = x10 * cos_turn - y10 * sin_turn;
  y1 = x10 * sin_turn + y10 * cos_turn;
  x2 = x20 * cos_turn - y20 * sin_turn;
  y2 = x20 * sin_turn + y20 * cos_turn;
  x3 = x30 * cos_turn - y30 * sin_turn;
  y3 = x30 * sin_turn + y30 * cos_turn;
  printf("\n\nThe base triangle angles: A = %.4fπ , B = %.4fπ , C = %.4fπ.\n\n", A/pi, B/pi, C/pi);
  printf("The following plots show slices of the cube [0,π] x [0,π] x [0,π] whose coordinates are α, β and γ. A system\n");
  printf("of inequalities defines an \"allowable\" portion of this cube. The slices are divided into cells. Each cell is\n");
  printf("designated to be \"allowable\" or \"unallowable,\" based on the system of inequalities. However, a cell that has\n");
  printf("been designated to be \"unallowable\" might actually contain some of the allowable portion of the cube together\n");
  printf("with some of the unallowable portion of the cube, in which case calling the cell \"unallowable\" is an unfortunate\n");
  printf("mistake. This can only happen at the boundary of the allowable portion of the cube.\n\n"); 
  printf("The allowable portion of the cube bounds all of the points (α, β, γ) for which α, β and γ can be the angles at\n");
  printf("a point P = (x, y, z) that extends the triangle ABC to form a tetrahedron ABCP. If a cell contains such a point\n");
  printf("(α, β, γ), then we call it \"occupied;\" otherwise the cell is \"unoccupied.\"\n\n");
#ifdef SHOW_ARRAY
  printf("Each cell is represented by a character. A space character represents an unoccupied allowable cell, an \'o\'\n");
  printf("represents an occupied allowable cell, a dot represents an unoccupied unallowable cell, and an \'x\' represents\n");
  printf("an occupied unallowable cell. This latter case is possible since an \"unallowable\" cell might contain an allowable\n");
  printf("portion of the cube (when it contains part of the boundary). When enabled, pound signs show where α = A, β = B or γ = C.\n\n");
#endif
  printf("PLEASE WAIT (patience is a virtue) while data is being generated .... \n\n\n\n\n");
  // Use 3D array to record possible (alpha, beta, gamma) triples for given triangle
  for (int i=O; i<M-O; i++)
    for (int j=O; j<M-O; j++)
      for (int k=O; k<M-O; k++)
        if ( tilt_to_view_angles(i*pi/M, j*pi/M, k*pi/M, cosA, cosB, cosC, alpha, beta, gamma, rejected) ) {
          states[ind(alpha)][ind(beta)][ind(gamma)] = 1;
          if ( test_bounds(A, B, C, cosA, cosB, cosC, x1, x2, x3, y1, y2, y3, alpha, beta, gamma) )
            states[ind(alpha)][ind(beta)][ind(gamma)] = 3;
        }
  // Also use array to record which cells in the array are within system of bounds
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
        if (states[i][j][k] < 2) {
#ifdef REFINED
          accept = false;
          for (delta_i=1; !accept && delta_i<=REF_NUM; delta_i++)
            for (delta_j=1; !accept && delta_j<=REF_NUM; delta_j++)
              for (delta_k=1; delta_k<=REF_NUM; delta_k++) {
                alpha = (i+(double)delta_i/(1+REF_NUM))*pi/N;
                beta  = (j+(double)delta_j/(1+REF_NUM))*pi/N;
                gamma = (k+(double)delta_k/(1+REF_NUM))*pi/N;
                if ( test_bounds(A, B, C, cosA, cosB, cosC, x1, x2, x3, y1, y2, y3, alpha, beta, gamma) )
                  { accept = true; break; }
              }
          if (accept) states[i][j][k] += 2;
#else
          alpha = (i+0.5)*pi/N;
          beta  = (j+0.5)*pi/N;
          gamma = (k+0.5)*pi/N;
          if ( test_bounds(A, B, C, cosA, cosB, cosC, x1, x2, x3, y1, y2, y3,
            alpha, beta, gamma) ) states[i][j][k] += 2;
#endif
        }
  // Show slices of the array, indicating the nature of each cell.
#ifdef SHOW_ARRAY
  show_array(states, i0, j0, k0);
#endif
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
