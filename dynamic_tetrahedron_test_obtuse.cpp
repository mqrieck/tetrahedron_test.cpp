// dynamic_tetrahedon_test_obtuse.cpp (by M. Q. Rieck, updated: 5/29/2022)

// Note: This is test code for the obtuse base triangle case, which uses far better bounds
// than those in my dynamic_tetrahedron_test.cpp, at least for the obtuse case. Further
// work is needed however, but I am definitely getting close now.

// Note: You need the ncurses library, and should then be able to compile at the command
// line using something like this:
//
//     g++ -c dynamic_tetrahedron_test_obtuse.cpp
//     g++ dynamic_tetrahedron_test_obtuse.o -lncurses -o test
//     ./test
//
// (You will need a sufficiently large virtual terminal screen and a small enough font.)

// Note: You must supply three command line integer parameters to specify angle proportion
// A:B:C. These must of course be positive, but to impose the obtuse case, it is also
// required that A > B + C, and so A will always be the obtuse angle.

// Note: This C++ program uses passing-by-reference. It can be easily converted to a C
// program by altering this aspect of function call, and by changing the includes.

// Note: For faster results, reduce M, N and/or REF_NUM (or comment out REFINED).

//#define M 3000                // how many (alpha, beta, gamma) points (M^3)?
//#define N 100                 // how fine to subdivide the interval [0, PI]
#define M 500                   // how many (alpha, beta, gamma) points (M^3)?
#define N 100                   // how fine to subdivide the interval [0, PI]
#define O 0                     // set this higher to avoid low "tilt planes"
#define PI M_PI                 // PI = 3.141592654..., of course
#define TOL1 0.05               // tolerance for some inequalities
#define TOL2 0                  // tolerance for some other inequalities
#define RESTRICT_DATA           // reject potentially troublesome data points
#define EXTRA_LOW_PLANES        // use more low elevation tilt planes
#define REFINED                 // more refined testing for cell acceptance/rejection
#define REF_NUM 9               // how much refinement?
//#define SHOW_EXTRA            // display a couple significant regions
//#define SHOW_MORE_DISCR
//#define SHOW_SPECIAL_PTS      // display special points
#define STARTX 2                // horizontal start of displayed character grid
#define STARTY 2                // vertical start of displayed character grid

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ncurses.h>
#ifdef COMPLEX_GRUNERT_DISCR
#include <complex>
#endif

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
    if (alpha < 0 || alpha > PI || beta < 0 || beta > PI || gamma < 0 || gamma > PI ||
      alpha > beta+gamma || beta > gamma+alpha || gamma > alpha+beta || alpha+beta+
        gamma > 2*PI) { rejected++; return false; }
#endif
    return true;
}

inline double f(double x) {
  return PI*(1-cos(x))/2;
}

inline int ind(double angle) {
  int i = (int) (N*angle/PI);
  if (i < 0) i = 0;
  if (i >= N) i = N-1;
  return i;
}

// Get cosine values for a point on the danger cylinder
inline void get_cosines(double phi1, double phi2, double phi3, double theta, double Z, double& c1, double& c2, double& c3) {
  double r1, r2, r3;
  r1 = sqrt(2-2*cos(theta-phi1)+Z);
  r2 = sqrt(2-2*cos(theta-phi2)+Z);
  r3 = sqrt(2-2*cos(theta-phi3)+Z);
  c1 = (1 + cos(phi2-phi3) - cos(theta-phi2) - cos(theta-phi3) + Z) / (r2*r3);
  c2 = (1 + cos(phi3-phi1) - cos(theta-phi3) - cos(theta-phi1) + Z) / (r3*r1);
  c3 = (1 + cos(phi1-phi2) - cos(theta-phi1) - cos(theta-phi2) + Z) / (r1*r2);
}

int main(int argc, char **argv) {
  int states[N][N][N], state, total, count0, count1, count2, count3, rejected = 0, i0, j0, k0, i, j, k, x, y,
    choice, delta_i, delta_j, delta_k, i12, i13, i22, i23, i32, i33, i42, i43, i52, i53, i62, i63, i72, i73, i82, i83,
    i92, i93, i102, i103, itemp;
  double A, B, C, cosA, cosB, cosC, alpha, beta, gamma, cos_alpha, cos_beta, cos_gamma, den,
    x10, x20, x30, y10, y20, y30, x1, x2, x3, y1, y2, y3, cos_turn, sin_turn, c1, c2, c3, C0, C1, C2, C3,
    H, L, R, E, D, t0, t1, t2, t3, sr, t10, t20, t30, G1, G2, G3, phi1, phi2, phi3, theta0, theta1, theta2,
    Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, c11, c12, c13, c21, c22, c23, c31, c32, c33, c41, c42, c43, c51, c52, c53,
    c61, c62, c63, c71, c72, c73, c81, c82, c83, c91, c92, c93, c101, c102, c103, sqt, temp;
  char ch, chars[N][N][N];
  bool accept, all_done, flags1[N][N][N], flags2[N][N][N], got1, got2, got3, got4, got5, got6, got7, got8, got9, got10;
#ifdef COMPLEX_GRUNERT_DISCR
  complex<double> two = 2, four = 4, eighteen = 18, twenty_seven = 27, zeta1, zeta2, zeta3, zeta1conj, zeta2conj,
    zeta3conj, zeta_prod, zeta_prod_sqr, zeta_prod_sqr_conj, xi, xi_conj, xi_norm, xi_cubed, xi_cubed_conj, D_complex;
#endif
  // Set the angles for the base triangle ABC
  // Can use three command line integer parameters to specify the proportion A : B : C
  if (argc == 4) {
    den = atoi(argv[1]) + atoi(argv[2]) + atoi(argv[3]);
    A = atoi(argv[1])*PI / den;
    B = atoi(argv[2])*PI / den;
    C = atoi(argv[3])*PI / den;
  } else {
    printf("\nYou must supply three command line integers. (See the program header.)\n\n");
    exit(0);
  }
  if (A < 0 || B < 0 || C < 0 || A < B+C) {
    printf("\nPlease read the program header for the parameter requirements.\n\n");
    exit(0);
  }
  cosA = cos(A); cosB = cos(B); cosC = cos(C);
  i0 = ind(A); j0 = ind(B); k0 = ind(C);
  if (i0 < 0) i0 = 0; if (i0 >= N) i0 = N-1;
  if (j0 < 0) j0 = 0; if (j0 >= N) j0 = N-1;
  if (k0 < 0) k0 = 0; if (k0 >= N) k0 = N-1;
// Set control points (at least initially)
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
  phi1 = atan2(y1,x1); phi2 = atan2(y2,x2); phi3 = atan2(y3,x3);
  printf("\n\nThe base triangle angles: A = %.4fπ , B = %.4fπ , C = %.4fπ.\n\n", A/PI, B/PI, C/PI);
  printf("The following plots show slices of the cube [0,π] x [0,π] x [0,π] whose coordinates are α, β and γ. A system\n");
  printf("of inequalities defines an \"allowable\" portion of this cube. The slices are divided into cells. Each cell is\n");
  printf("designated to be \"allowable\" or \"unallowable,\" based on the system of inequalities. However, a cell that has\n");
  printf("been designated to be \"unallowable\" might actually contain some of the allowable portion of the cube together\n");
  printf("with some of the unallowable portion of the cube, in which case calling the cell \"unallowable\" is an unfortunate\n");
  printf("mistake. Such mistakes can only happen at the boundary of the allowable portion of the cube.\n\n");
  printf("The allowable portion of the cube bounds all of the points (α, β, γ) for which α, β and γ can be the angles at\n");
  printf("a point P = (x, y, z) that extends the triangle ABC to form a tetrahedron ABCP. If a cell contains such a point\n");
  printf("(α, β, γ), then we call it \"occupied;\" otherwise the cell is \"unoccupied.\" (A basic understanding of the problem\n");
  printf("in the research paper that this program supports is presumed here.)\n\n");
  printf("Each cell is represented by a character. A space character represents an unoccupied allowable cell, an \'o\'\n");
  printf("represents an occupied allowable cell, a dot represents an unoccupied unallowable cell, and an \'x\' represents\n");
  printf("an occupied unallowable cell. This latter case is possible since an \"unallowable\" cell might contain an allowable\n");
  printf("portion of the cube (when it contains part of the boundary).\n\n");
  printf("PLEASE WAIT (patience is a virtue) while data is being generated (press the enter/return key if stuck) ....\n\n\n\n\n");
  // Use a 3D array to record possible (alpha, beta, gamma) triples for given triangle
  for (i=O; i<M-O; i++)
    for (j=O; j<M-O; j++)
      for (k=O; k<M-O; k++)
#ifdef EXTRA_LOW_PLANES
        if (tilt_to_view_angles(f(i*PI/M), f(j*PI/M), f(k*PI/M), cosA, cosB, cosC, alpha,
          beta, gamma, rejected)) states[ind(alpha)][ind(beta)][ind(gamma)] = 1;
#else
        if (tilt_to_view_angles(i*PI/M, j*PI/M, k*PI/M, cosA, cosB, cosC, alpha,
          beta, gamma, rejected)) states[ind(alpha)][ind(beta)][ind(gamma)] = 1;
#endif
  // Also use an array to record which cells in the array are within system of bounds
  for (i=0; i<N; i++) {
    alpha = (i+0.5)*PI/N;
    c1 = cos_alpha = cos(alpha);
    C1 = c1*c1;
    got1 = got2 = got3 = got4 = got5 = got6 = got7 = got8 = got9 = got10 = false;
// Self-crossing points are needed for the D = 0 curve when B+C <= alpha < A
    if (fabs(cosC*c1/cosA) <= 1) {
      c12 = cosB; c13 = -cosC*c1/cosA; i12 = ind(B); i13 = ind(acos(c13));
      states[i][i12][i13] = 10 + states[i][i12][i13] % 10;
      got1 = true;
    }
    if (fabs(cosB*c1/cosA) <= 1) {
      c22 = -cosB*c1/cosA; c23 = cosC; i22 = ind(acos(c22)); i23 = ind(C);
      states[i][i22][i23] = 20 + states[i][i22][i23] % 10;
      got2 = true;
    }
    if (alpha >= B+C && alpha < A) {
      theta0 = acos((C1-cos(phi2-phi3))/(1-C1));
      theta1 = phi1 + theta0;
      theta2 = phi1 - theta0;
      sqt = (2*cos(theta1)-1)*cos(theta1/2+phi1)*cos(theta1/2+phi2)*cos(theta1/2+phi3)/cos(theta1/2);
      if (sqt >= 0) {
        sqt = sqrt(sqt);
        Z3 = cos(theta1-phi1) + cos(theta1-phi2) + cos(theta1-phi3) - 1 + (
             12 + 16*cos(theta1) + 4*cos(2*theta1) + 8*cos(theta1-phi2) + 2*cos(2*theta1-phi2) + 12*cos(phi2) + 8*cos(theta1+phi2)
             + 2*cos(2*theta1+phi2) + 8*cos(theta1-phi3) + 2*cos(2*theta1-phi3) + 4*cos(theta1+phi1) + cos(2*theta1+phi1) + 6*cos(phi2-phi3)
             + 4*cos(theta1+phi2-phi3) + cos(2*theta1+phi2-phi3) + 12*cos(phi3) + 8*cos(theta1+phi3) + 2*cos(2*theta1+phi3)
             + 4*cos(theta1-phi2+phi3) + cos(2*theta1-phi2+phi3) + 6*cos(phi1) + 4*cos(theta1-phi1) + cos(2*theta1-phi1)
        ) * sqt / (8*(1+cos(phi2))*(1+cos(phi3))*(1+cos(theta1))*cos(theta1/2));
        if (Z3 >= 0) {
          get_cosines(phi1, phi2, phi3, theta1, Z3, c31, c32, c33);
          if (fabs(c1+c31) < fabs(c1-c31)) {c31 = -c31; c32 = -c32;}
          if (c32 < 0) {c32 = -c32; c33 = -c33;}
          i32 = ind(acos(c32)); i33 = ind(acos(c33));
          states[i][i32][i33] = 30 + states[i][i32][i33] % 10;
          got3 = true;
//printf("3.  %d  %f  %f %f %f\n", states[i][i32][i33], c1, c31, c32, c33);
        }
      }
      sqt = (2*cos(theta2)-1)*cos(theta2/2+phi1)*cos(theta2/2+phi2)*cos(theta2/2+phi3)/cos(theta2/2);
      if (sqt >= 0) {
        sqt = sqrt(sqt);
        Z4 = cos(theta2-phi1) + cos(theta2-phi2) + cos(theta2-phi3) - 1 + (
             12 + 16*cos(theta2) + 4*cos(2*theta2) + 8*cos(theta2-phi2) + 2*cos(2*theta2-phi2) + 12*cos(phi2) + 8*cos(theta2+phi2)
             + 2*cos(2*theta2+phi2) + 8*cos(theta2-phi3) + 2*cos(2*theta2-phi3) + 4*cos(theta2+phi1) + cos(2*theta2+phi1) + 6*cos(phi2-phi3)
             + 4*cos(theta2+phi2-phi3) + cos(2*theta2+phi2-phi3) + 12*cos(phi3) + 8*cos(theta2+phi3) + 2*cos(2*theta2+phi3)
             + 4*cos(theta2-phi2+phi3) + cos(2*theta2-phi2+phi3) + 6*cos(phi1) + 4*cos(theta2-phi1) + cos(2*theta2-phi1)
        ) * sqt / (8*(1+cos(phi2))*(1+cos(phi3))*(1+cos(theta2))*cos(theta2/2));
        if (Z4 >= 0) {
          get_cosines(phi1, phi2, phi3, theta2, Z4, c41, c42, c43);
          if (fabs(c1+c41) < fabs(c1-c41)) {c41 = -c41; c42 = -c42;}
          if (c43 < 0) {c42 = -c42; c43 = -c43;}
          i42 = ind(acos(c42)); i43 = ind(acos(c43));
          states[i][i42][i43] = 40 + states[i][i42][i43] % 10;
          got4 = true;
//printf("4.  %d  %f  %f %f %f\n", states[i][i42][i43], c1, c41, c42, c43);
        }
      }
// Other singular points on the D = 0 curve
      sqt = 4*C1*(1-cos(phi2-phi3))*(C1-cos(phi2-phi3)+(1-C1)*cos(phi2+phi3));
      if (sqt >= 0) {
        sqt = sqrt(sqt) / 2;
        Z5 = (2*C1 - 1 - cos(phi2-phi3) + (C1-1)*(cos(phi2)+cos(phi3)) + sqt) / (1 - C1);
        if (Z5 >= 0) {
          get_cosines(phi1, phi2, phi3, PI, Z5, c51, c52, c53);
          if (fabs(c1+c51) < fabs(c1-c51)) {c51 = -c51; c52 = -c52;}
          if (c53 < 0) {c52 = -c52; c53 = -c53;}
          i52 = ind(acos(c52)); i53 = ind(acos(c53));
          states[i][i52][i53] = 50 + states[i][i52][i53] % 10;
          got5 = true;
//printf("5.  %d  %f  %f %f %f\n", states[i][i52][i53], c1, c51, c52, c53);
        }
        Z6 = (2*C1 - 1 - cos(phi2-phi3) + (C1-1)*(cos(phi2)+cos(phi3)) - sqt) / (1 - C1);
        if (Z6 >= 0) {
          get_cosines(phi1, phi2, phi3, PI, Z6, c61, c62, c63);
          if (fabs(c1+c61) < fabs(c1-c61)) {c61 = -c61; c62 = -c62;}
          if (c62 < 0) {c62 = -c62; c63 = -c63;}
          i62 = ind(acos(c62)); i63 = ind(acos(c63));
          states[i][i62][i63] = 60 + states[i][i62][i63] % 10;
          got6 = true;
//printf("6.  %d  %f  %f %f %f\n", states[i][i62][i63], c1, c61, c62, c63);
        }
      }
      sqt = 2*(C1 + 2*C1*C1 + C1*cos(2*(phi2-phi3)) - 2*C1*(1+C1)*cos(phi2-phi3) - 2*C1*(1-C1)*sin(PI/6-phi2-phi3)
            + C1*(1-C1)*(sin(PI/6-2*phi2) + sin(PI/6-2*phi3)));
      if (sqt >= 0) {
        sqt = sqrt(sqt) / 2;
        Z7 = (2*C1 - 1 - cos(phi2-phi3) + (1-C1)*(sin(PI/6+phi2)+sin(PI/6+phi3)) + sqt) / (1 - C1);
        if (Z7 >= 0) {
          get_cosines(phi1, phi2, phi3, PI/3, Z7, c71, c72, c73);
          if (fabs(c1+c71) < fabs(c1-c71)) {c71 = -c71; c72 = -c72;}
          if (c72 < 0) {c72 = -c72; c73 = -c73;}
          i72 = ind(acos(c72)); i73 = ind(acos(c73));
          states[i][i72][i73] = 70 + states[i][i72][i73] % 10;
          got7 = true;
//printf("7.  %d  %f  %f %f %f\n", states[i][i72][i73], c1, c71, c72, c73);
        }
        Z8 = (2*C1 - 1 - cos(phi2-phi3) + (1-C1)*(sin(PI/6+phi2)+sin(PI/6+phi3)) - sqt) / (1 - C1);
        if (Z8 >= 0) {
          get_cosines(phi1, phi2, phi3, PI/3, Z8, c81, c82, c83);
          if (fabs(c1+c81) < fabs(c1-c81)) {c81 = -c81; c82 = -c82;}
          if (c82 < 0) {c82 = -c82; c83 = -c83;}
          i82 = ind(acos(c82)); i83 = ind(acos(c83));
          states[i][i82][i83] = 80 + states[i][i82][i83] % 10;
          got8 = true;
//printf("8.  %d  %f  %f %f %f\n", states[i][i82][i83], c1, c81, c82, c83);
        }
      }
      sqt = 2*(C1 + 2*C1*C1 + C1*cos(2*(phi2-phi3)) - 2*C1*(1+C1)*cos(phi2-phi3) - 2*C1*(1-C1)*sin(PI/6+phi2+phi3)
            + C1*(1-C1)*(sin(PI/6+2*phi2) + sin(PI/6+2*phi3)));
      if (sqt >= 0) {
        sqt = sqrt(sqt) / 2;
        Z9 = (2*C1 - 1 - cos(phi2-phi3) + (1-C1)*(sin(PI/6-phi2)+sin(PI/6-phi3)) + sqt) / (1 - C1);
        if (Z9 >= 0) {
          get_cosines(phi1, phi2, phi3, -PI/3, Z9, c91, c92, c93);
          if (fabs(c1+c91) < fabs(c1-c91)) {c91 = -c91; c92 = -c92;}
          if (c93 < 0) {c92 = -c92; c93 = -c93;}
          i92 = ind(acos(c92)); i93 = ind(acos(c93));
          states[i][i92][i93] = 90 + states[i][i92][i93] % 10;
          got9 = true;
//printf("9.  %d  %f  %f %f %f\n", states[i][i92][i93], c1, c91, c92, c93);
        }
        Z10 = (2*C1 - 1 - cos(phi2-phi3) + (1-C1)*(sin(PI/6-phi2)+sin(PI/6-phi3)) - sqt) / (1 - C1);
        if (Z10 >= 0) {
          get_cosines(phi1, phi2, phi3, -PI/3, Z10, c101, c102, c103);
          if (fabs(c1+c101) < fabs(c1-c101)) {c101 = -c101; c102 = -c102;}
          if (c103 < 0) {c102 = -c102; c103 = -c103;}
          i102 = ind(acos(c102)); i103 = ind(acos(c103));
          states[i][i102][i103] = 100 + states[i][i102][i103] % 10;
          got10 = true;
//printf("10. %d  %f  %f %f %f\n", states[i][i102][i103], c1, c101, c102, c103);
        }
      }
      if (got7 && got8 && c73 > c83) {
        temp = c72; c72 = c82; c82 = temp;
        temp = c73; c73 = c83; c83 = temp;
        itemp = i72; i72 = i82; i82 = itemp;
        itemp = i73; i73 = i83; i83 = itemp;
        states[i][i72][i73] -= 10;
        states[i][i82][i83] += 10;
      }
      if (got9 && got10 && c92 > c102) {
        temp = c92; c92 = c102; c102 = temp;
        temp = c93; c93 = c103; c103 = temp;
        itemp = i92; i92 = i102; i102 = itemp;
        itemp = i93; i93 = i103; i103 = itemp;
        states[i][i92][i93] -= 10;
        states[i][i102][i103] += 10;
      }
    }
    for (j=0; j<N; j++)
      for (k=0; k<N; k++) {
#ifdef REFINED
        accept = false;
        for (delta_i=1; !accept && delta_i<=REF_NUM; delta_i++)
          for (delta_j=1; !accept && delta_j<=REF_NUM; delta_j++)
            for (delta_k=1; delta_k<=REF_NUM; delta_k++) {
              alpha = (i+(double)delta_i/(1+REF_NUM))*PI/N;
              beta  = (j+(double)delta_j/(1+REF_NUM))*PI/N;
              gamma = (k+(double)delta_k/(1+REF_NUM))*PI/N;
#else
              beta  = (j+0.5)*PI/N;
              gamma = (k+0.5)*PI/N;
#endif
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
#ifdef SHOW_MORE_DISCR
              flags1[i][j][k] = 0;
#else
              flags1[i][j][k] = (H < 0);
#endif
              flags2[i][j][k] = (D < 0);

// Test the constraints in my papers concerning the angles A, B, C, alpha, beta and gamma
              if (
// sphere-based constraints
                alpha +  beta + gamma < 2*PI &&
                alpha <  beta + gamma &&
                beta  < gamma + alpha &&
                gamma < alpha +  beta
                &&
                A + beta + gamma  < 2*PI &&
                alpha + B + gamma < 2*PI &&
                alpha + beta + C  < 2*PI
                &&
                beta  + gamma - alpha < 2*(B+C) &&
                gamma + alpha -  beta < 2*(C+A) &&
                alpha +  beta - gamma < 2*(A+B)
// small alpha constraints
		&& (
                  !(alpha <= B+C) || (
                  (cosC * cos_beta  + cosB * cos_gamma > 0) &&
                  beta < B || gamma < C || (D < 0 &&
                    (alpha + B - C) * beta + (alpha - B + C) * gamma < (alpha + B + C) * alpha)
                ))
// trick border case
                && (
                  !(alpha < A) || got3 || got4 || got7 || got9 || (
                    (cosC * cos_beta  + cosB * cos_gamma > 0) &&
                    beta < B || gamma < C || (D < 0 &&
                      (alpha + B - C) * beta + (alpha - B + C) * gamma < (alpha + B + C) * alpha)
                ))

// midrange alpha constraints
                && (
                  !(alpha > B+C && alpha < A) || (
                  (cosC * cos_beta  + cosB * cos_gamma > 0) &&
                  (beta  >= B || gamma < C || gamma < A + beta ) &&
                  (beta  >= B || alpha < A || alpha < C + beta ) &&
                  (gamma >= C || alpha < A || alpha < B + gamma) &&
                  (gamma >= C || beta  < B || beta  < A + gamma) &&
                  (alpha < B + C || alpha > A || (beta < C + alpha && gamma < B + alpha)) &&
                  (alpha + B - C) * beta + (alpha - B + C) * gamma < (alpha + B + C) * alpha &&
                  (!got3 || c2 >= c32 || c3 >= c33) &&
                  (!got4 || c2 >= c42 || c3 >= c43) &&
                  (!got7 || c2 >= c72 || c3 >= c73) &&
                  (!got9 || c2 >= c92 || c3 >= c93) &&
                  (!got3 || !got4 || c2 > c32 || c2 < c42 || c3 < c33 || c3 > c43 || D >= 0) &&
                  (!got7 || !got3 || c2 > c72 || c2 < c32 || c3 < c73 || c3 > c33 || D <= 0) &&
                  (!got4 || !got9 || c2 > c42 || c2 < c92 || c3 < c43 || c3 > c93 || D <= 0)
                ))
// large alpha constraints
                && (
                  !(alpha >= A) || (beta >= B && gamma >= C)
                )
#ifdef REFINED
              ) { accept = true; break; }
        }
        if (accept && states[i][j][k] >= 0) states[i][j][k] += 2;
#else
              && states[i][j][k] >= 0) states[i][j][k] += 2;
#endif
        switch(states[i][j][k]) {
          case   0:  chars[i][j][k] = '.'; break;
          case   1:  chars[i][j][k] = 'x'; break;
          case   2:  chars[i][j][k] = ' '; break;
          case   3:  chars[i][j][k] = 'o'; break;
          case  10: case  11: case  12: case  13:  chars[i][j][k] = '1'; break;
          case  20: case  21: case  22: case  23:  chars[i][j][k] = '2'; break;
          case  30: case  31: case  32: case  33:  chars[i][j][k] = '3'; break;
          case  40: case  41: case  42: case  43:  chars[i][j][k] = '4'; break;
          case  50: case  51: case  52: case  53:  chars[i][j][k] = '5'; break;
          case  60: case  61: case  62: case  63:  chars[i][j][k] = '6'; break;
          case  70: case  71: case  72: case  73:  chars[i][j][k] = '7'; break;
          case  80: case  81: case  82: case  83:  chars[i][j][k] = '8'; break;
          case  90: case  91: case  92: case  93:  chars[i][j][k] = '9'; break;
          case 100: case 101: case 102: case 103:  chars[i][j][k] = '0'; break;
          default:  chars[i][j][k] = '?';
        }
    }
  }
  initscr();
  start_color();
  init_pair(1,  COLOR_WHITE,   COLOR_BLACK);
  init_pair(2,  COLOR_BLUE,    COLOR_WHITE);
  init_pair(3,  COLOR_RED,     COLOR_WHITE);
  init_pair(4,  COLOR_MAGENTA, COLOR_WHITE);
  init_pair(5,  COLOR_BLACK,   COLOR_WHITE);
  init_pair(6,  COLOR_BLUE,    COLOR_CYAN);
  init_pair(7,  COLOR_RED,     COLOR_CYAN);
  init_pair(8,  COLOR_MAGENTA, COLOR_CYAN);
  init_pair(9,  COLOR_BLUE,    COLOR_GREEN);
  init_pair(10, COLOR_RED,     COLOR_GREEN);
  init_pair(11, COLOR_MAGENTA, COLOR_GREEN);
  init_pair(12, COLOR_BLUE,    COLOR_YELLOW);
  init_pair(13, COLOR_RED,     COLOR_YELLOW);
  init_pair(14, COLOR_MAGENTA, COLOR_YELLOW);
  init_pair(15, COLOR_WHITE,   COLOR_RED);
  curs_set(0);
  cbreak();
  noecho();
  keypad(stdscr,TRUE);
  attron(COLOR_PAIR(1));
  printf("\033[2J");    // Clears the screen in a Unix-like OS
  nodelay(stdscr,TRUE);
  attron(COLOR_PAIR(1));
  mvprintw(STARTY   , STARTX+2*N+3, "A = %.4f", A);
  mvprintw(STARTY+ 1, STARTX+2*N+3, "B = %.4f", B);
  mvprintw(STARTY+ 2, STARTX+2*N+3, "C = %.4f", C);
  mvprintw(STARTY+ 6, STARTX+2*N+3, "Use up and down arrow keys to view different slices of the bounding region.");
  mvprintw(STARTY+ 7, STARTX+2*N+3, "Dots are outside the bounding region; spaces are inside but without data");
  mvprintw(STARTY+ 8, STARTX+2*N+3, "points to bound; o's are inside and contain data points; x's are subtle.");
  mvprintw(STARTY+10, STARTX+2*N+3, "Use up and down arrow keys to view different slices of bounding region.");
  mvprintw(STARTY+11, STARTX+2*N+3, "You can also select to display one of these special lines/curves");
  mvprintw(STARTY+12, STARTX+2*N+3, "(imperfectly rendered):");
  mvprintw(STARTY+14, STARTX+2*N+3, "1. beta = B");
  mvprintw(STARTY+15, STARTX+2*N+3, "2. gamma = C");
  mvprintw(STARTY+16, STARTX+2*N+3, "3. beta + gamma = 2 PI - alpha");
  mvprintw(STARTY+17, STARTX+2*N+3, "4. beta + gamma = alpha");
  mvprintw(STARTY+18, STARTX+2*N+3, "5. beta - gamma = alpha");
  mvprintw(STARTY+19, STARTX+2*N+3, "6. gamma - beta = alpha");
  mvprintw(STARTY+20, STARTX+2*N+3, "7. beta + gamma = 2 PI - A");
  mvprintw(STARTY+21, STARTX+2*N+3, "8. beta = 2 PI - alpha - C");
  mvprintw(STARTY+22, STARTX+2*N+3, "9. gamma = 2 PI - alpha - B");
  mvprintw(STARTY+23, STARTX+2*N+3, "a. beta = C + alpha");
  mvprintw(STARTY+24, STARTX+2*N+3, "b. gamma = B + alpha");
  mvprintw(STARTY+25, STARTX+2*N+3, "c. beta - gamma = A");
  mvprintw(STARTY+26, STARTX+2*N+3, "d. gamma - beta = A");
  mvprintw(STARTY+27, STARTX+2*N+3, "e. cos C cos beta + cos B cos gamma = 0");
  mvprintw(STARTY+28, STARTX+2*N+3, "f. cos A cos gamma + cos C cos alpha = 0");
  mvprintw(STARTY+29, STARTX+2*N+3, "g. cos B cos alpha + cos A cos beta = 0");
  mvprintw(STARTY+30, STARTX+2*N+3, "h. (alpha+B-C) beta + (alpha+C-B) gamma = alpha(alpha+B+C)");
  mvprintw(STARTY+32, STARTX+2*N+3, "Press the escape key to quit.");
  all_done = false;
  choice = i = 0;
  do {
    int ch = getch();
    if (ch == KEY_DOWN) { i = i>0   ? i-1 : i; choice = 0; }
    if (ch == KEY_UP)   { i = i<N-2 ? i+1 : i; choice = 0; }
    if (ch == '1') choice =  1;
    if (ch == '2') choice =  2;
    if (ch == '3') choice =  3;
    if (ch == '4') choice =  4;
    if (ch == '5') choice =  5;
    if (ch == '6') choice =  6;
    if (ch == '7') choice =  7;
    if (ch == '8') choice =  8;
    if (ch == '9') choice =  9;
    if (ch == 'a') choice = 10;
    if (ch == 'b') choice = 11;
    if (ch == 'c') choice = 12;
    if (ch == 'd') choice = 13;
    if (ch == 'e') choice = 14;
    if (ch == 'f') choice = 15;
    if (ch == 'g') choice = 16;
    if (ch == 'h') choice = 17;
    if (ch == 27) all_done = true;
    alpha = (i+.5)*PI/N;
    attron(COLOR_PAIR(1));
    mvprintw(STARTY+ 4, STARTX+2*N+3, "alpha = %.4f", alpha);
    if (alpha == A) mvprintw(STARTY+ 4, STARTX+2*N+20, "= A");
    if (alpha  < A) mvprintw(STARTY+ 4, STARTX+2*N+20, "< A");
    if (alpha  > A) mvprintw(STARTY+ 4, STARTX+2*N+20, "> A");
    for(j=0, x=STARTX; j < N; j++, x+=2)
      for(k=0, y=STARTY; k < N; k++, y++) {
        if (states[i][j][k] > 9) attron(COLOR_PAIR(15));
        else if (chars[i][j][k] == '.') {
#ifdef SHOW_EXTRA
          if (flags1[i][j][k]) attron(COLOR_PAIR(7)); else
          if (flags2[i][j][k]) attron(COLOR_PAIR(10)); else
          attron(COLOR_PAIR(13));
#else
          attron(COLOR_PAIR(3));
#endif
        } else if (chars[i][j][k] == 'x') {
#ifdef SHOW_EXTRA
          if (flags1[i][j][k]) attron(COLOR_PAIR(8)); else
          if (flags2[i][j][k]) attron(COLOR_PAIR(11)); else
          attron(COLOR_PAIR(14));
#else
          attron(COLOR_PAIR(4));
#endif
        } else {
#ifdef SHOW_EXTRA
          if (flags1[i][j][k]) attron(COLOR_PAIR(6)); else
          if (flags2[i][j][k]) attron(COLOR_PAIR(9)); else
          attron(COLOR_PAIR(12));
#else
          attron(COLOR_PAIR(2));
#endif
        }
        mvprintw(y,x  ,"%c",chars[i][j][k]);
        mvprintw(y,x+1,"%c",chars[i][j][k]);
      }

#ifdef SHOW_SPECIAL_PTS
    attron(COLOR_PAIR(5));
    x = STARTX + 2*ind(B);
    y = STARTY + ind(acos(-cos(alpha)*cosC/cosA));
    mvprintw(y,x  ,"%c",'>');
    mvprintw(y,x+1,"%c",'<');
    x = STARTX + 2*ind(acos(-cos(alpha)*cosB/cosA));
    y = STARTY + ind(C);
    mvprintw(y,x  ,"%c",'>');
    mvprintw(y,x+1,"%c",'<');
#endif
    if (choice == 1) {
      attron(COLOR_PAIR(5));
      x = STARTX + 2*j0;
      for(k=0, y=STARTY; k < N; k++, y++) {
         mvprintw(y,x  ,"%c",'%');
         mvprintw(y,x+1,"%c",'%');
      }
    }
    if (choice == 2) {
      attron(COLOR_PAIR(5));
      y = STARTY + k0;
      for(j=0, x=STARTX; j < N; j++, x+=2) {
        mvprintw(y,x  ,"%c",'%');
        mvprintw(y,x+1,"%c",'%');
      }
    }
    if (choice >= 3) {
      attron(COLOR_PAIR(5));
      for(j=0, x=STARTX; j < N; j++, x+=2) {
        beta = (j+.5)*PI/N;
        for(k=0, y=STARTY; k < N; k++, y++) {
          gamma = (k+.5)*PI/N;
          switch(choice) {
             case 3: if (fabs(alpha + beta + gamma - 2*PI) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  4: if (fabs(beta + gamma - alpha) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  5: if (fabs(gamma + alpha - beta) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  6: if (fabs(alpha + beta - gamma) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  7: if (fabs(A + beta + gamma - 2*PI) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  8: if (fabs(alpha + B + gamma - 2*PI) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  9: if (fabs(alpha + beta + C - 2*PI) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 10: if (fabs(alpha - beta + C) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 11: if (fabs(alpha + B - gamma) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 12: if (fabs(A - beta + gamma) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 13: if (fabs(A + beta - gamma) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 14: if (fabs(cos(C)*cos(beta) + cos(B)*cos(gamma)) < TOL1/2) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 15: if (fabs(cos(A)*cos(gamma) + cos(C)*cos(alpha)) < TOL1/2) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 16: if (fabs(cos(B)*cos(alpha) + cos(A)*cos(beta)) < TOL1/2) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 17: if (fabs((alpha+B-C)*beta + (alpha+C-B)*gamma - alpha*(alpha+B+C)) < TOL1) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
          }
        }
      }
    }
    refresh();
  } while (!all_done);
  attroff(COLOR_PAIR(1));
  endwin();
}