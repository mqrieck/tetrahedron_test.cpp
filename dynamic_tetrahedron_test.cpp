
// dynamic_tetrahedon_test.cpp (by M. Q. Rieck, updated: 11/9/2021)

// Note: This is test code for the results in my "tetrahedron and toroids" paper, and more.

// Note: Need the ncurses library, and should be able to compile at the command line
// using something like this:
//
//     g++ -c extended_dynamic_test.cpp
//     g++ extended_dymanic_test.o -lncurses -o test
//     ./test
//
// (Will need sufficiently large virtual terminal screen and small enough font.)

// Note: Can use three command line integer parameters to specify angle proportion A:B:C.

// Note: This C++ program uses passing-by-reference. It can be easily converted to a C
// program by altering this aspect of function call, and by changing the includes.

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ncurses.h>

#define M 1600                  // how many (alpha, beta, gamma) points (M^3)?
#define N 100                   // how fine to subdivide the interval [0, pi]
#define O 0                     // set higher to avoid low "tilt planes"
#define pi M_PI                 // pi = 3.141592654..., of course
#define ACUTE_TEST              // only appropriate for acute base triangle ABC
#define COSINES_TEST            // include the "cosines test" when using an acute triangle
#define NEW_TESTS		// use the new tests based on Grunert's system discriminant
#define REFINED                 // more refined testing for accepting a cell
//#define SHOW_EXTRA            // display a couple significant regions
#define STARTX 2                // horizontal start of displayed character grid
#define STARTY 2                // vertical start of displayed character grid

using namespace std;

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

inline int ind(double angle) {
  int i = (int) (N*angle/pi);
  if (i < 0) i = 0;
  if (i >= N) i = N-1;
  return i;
}

int main(int argc, char **argv) {
  int states[N][N][N], state, total, count0, count1, count2, count3, rejected = 0, i0, j0, k0, i, j, k, x, y,
    choice, delta_i, delta_j, delta_k;
  double A, B, C, cosA, cosB, cosC, alpha, beta, gamma, cos_alpha, cos_beta, cos_gamma, den, tol = 0.05,
    x10, x20, x30, y10, y20, y30, x1, x2, x3, y1, y2, y3, cos_turn, sin_turn, c1, c2, c3, C0, C1, C2, C3,
    eta_sq, L, R, E, D;
  char ch, chars[N][N][N];
  bool accept, all_done, flags1[N][N][N], flags2[N][N][N];
  // Set the angles for the base triangle ABC
  // Can use three command line integer parameters to specify the proportion A : B : C
  if (argc == 4) {
    den = atoi(argv[1]) + atoi(argv[2]) + atoi(argv[3]);
    A = atoi(argv[1])*pi / den;
    B = atoi(argv[2])*pi / den;
    C = atoi(argv[3])*pi / den;
  } else A = B = C = pi/3;
  cosA = cos(A); cosB = cos(B); cosC = cos(C);
  i0 = ind(A); j0 = ind(B); k0 = ind(C);
  if (i0 < 0) i0 = 0; if (i0 >= N) i0 = N-1;
  if (j0 < 0) j0 = 0; if (j0 >= N) j0 = N-1;
  if (k0 < 0) k0 = 0; if (k0 >= N) k0 = N-1;
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
  printf("in the research paper that this program supports is presumed here.)\n\n");
  printf("Each cell is represented by a character. A space character represents an unoccupied allowable cell, an \'o\'\n");
  printf("represents an occupied allowable cell, a dot represents an unoccupied unallowable cell, and an \'x\' represents\n");
  printf("an occupied unallowable cell. This latter case is possible since an \"unallowable\" cell might contain an allowable\n");
  printf("portion of the cube (when it contains part of the boundary).\n\n");
  printf("PLEASE WAIT while data is being generated (press the enter key if stuck) ....\n\n\n\n\n");
  // Use 3D array to record possible (alpha, beta, gamma) triples for given triangle
  for (i=O; i<M-O; i++)
    for (j=O; j<M-O; j++)
      for (k=O; k<M-O; k++)
        if (tilt_to_view_angles(i*pi/M, j*pi/M, k*pi/M, cosA, cosB, cosC, alpha,
          beta, gamma, rejected)) states[ind(alpha)][ind(beta)][ind(gamma)] = 1;
  // Also use array to record which cells in the array are within system of bounds
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      for (k=0; k<N; k++) {
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
              flags1[i][j][k] = (eta_sq < 0);
              flags2[i][j][k] = (D < 0);
              if (
                alpha +  beta + gamma < 2*pi &&
                alpha <  beta + gamma &&
                beta  < gamma + alpha &&
                gamma < alpha +  beta &&
                A + beta + gamma  < 2*pi &&
                alpha + B + gamma < 2*pi &&
                alpha + beta + C  < 2*pi &&
                beta  + gamma - alpha < 2*(B+C) &&
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
#ifdef NEW_TESTS
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
        switch(states[i][j][k]) {
          case 0: chars[i][j][k] = '.'; break;
          case 1: chars[i][j][k] = 'x'; break;
          case 2: chars[i][j][k] = ' '; break;
          case 3: chars[i][j][k] = 'o'; break;
        }
  }
  initscr();
  start_color();
  init_pair(1,  COLOR_WHITE,   COLOR_BLACK);
  init_pair(2,  COLOR_BLUE,    COLOR_WHITE);
  init_pair(3,  COLOR_RED,     COLOR_WHITE);
  init_pair(4,  COLOR_GREEN,   COLOR_WHITE);
  init_pair(5,  COLOR_BLACK,   COLOR_WHITE);
  init_pair(6,  COLOR_BLUE,    COLOR_CYAN);
  init_pair(7,  COLOR_RED,     COLOR_CYAN);
  init_pair(8,  COLOR_GREEN,   COLOR_CYAN);
  init_pair(9,  COLOR_BLUE,    COLOR_MAGENTA);
  init_pair(10, COLOR_RED,     COLOR_MAGENTA);
  init_pair(11, COLOR_GREEN,   COLOR_MAGENTA);
  init_pair(12, COLOR_BLUE,    COLOR_YELLOW);
  init_pair(13, COLOR_RED,     COLOR_YELLOW);
  init_pair(14, COLOR_GREEN,   COLOR_YELLOW);
  curs_set(0);
  cbreak();
  noecho();
  keypad(stdscr,TRUE);
  attron(COLOR_PAIR(1));
  printf("\033[2J"); // Clears the screen! Uncomment for Unix(-like). Comment out for VC++.
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
  mvprintw(STARTY+16, STARTX+2*N+3, "3. beta + gamma = 2 pi - alpha");
  mvprintw(STARTY+17, STARTX+2*N+3, "4. beta + gamma = alpha");
  mvprintw(STARTY+18, STARTX+2*N+3, "5. beta - gamma = alpha");
  mvprintw(STARTY+19, STARTX+2*N+3, "6. gamma - beta = alpha");
  mvprintw(STARTY+20, STARTX+2*N+3, "7. beta + gamma = 2 pi - A");
  mvprintw(STARTY+21, STARTX+2*N+3, "8. beta = 2 pi - alpha - C");
  mvprintw(STARTY+22, STARTX+2*N+3, "9. gamma = 2 pi - alpha - B");
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
    alpha = (i+.5)*pi/N;
    attron(COLOR_PAIR(1));
    mvprintw(STARTY+ 4, STARTX+2*N+3, "alpha = %.4f", alpha);
    if (alpha == A) mvprintw(STARTY+ 4, STARTX+2*N+20, "= A");
    if (alpha  < A) mvprintw(STARTY+ 4, STARTX+2*N+20, "< A");
    if (alpha  > A) mvprintw(STARTY+ 4, STARTX+2*N+20, "> A");
    for(j=0, x=STARTX; j < N; j++, x+=2)
      for(k=0, y=STARTY; k < N; k++, y++) {
        if (chars[i][j][k] == '.') {
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
        beta = (j+.5)*pi/N;
        for(k=0, y=STARTY; k < N; k++, y++) {
          gamma = (k+.5)*pi/N;
          switch(choice) {
             case 3: if (fabs(alpha + beta + gamma - 2*pi) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  4: if (fabs(beta + gamma - alpha) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  5: if (fabs(gamma + alpha - beta) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  6: if (fabs(alpha + beta - gamma) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  7: if (fabs(A + beta + gamma - 2*pi) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  8: if (fabs(alpha + B + gamma - 2*pi) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case  9: if (fabs(alpha + beta + C - 2*pi) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 10: if (fabs(alpha - beta + C) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 11: if (fabs(alpha + B - gamma) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 12: if (fabs(A - beta + gamma) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 13: if (fabs(A + beta - gamma) < tol) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 14: if (fabs(cos(C)*cos(beta) + cos(B)*cos(gamma)) < tol/2) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 15: if (fabs(cos(A)*cos(gamma) + cos(C)*cos(alpha)) < tol/2) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 16: if (fabs(cos(B)*cos(alpha) + cos(A)*cos(beta)) < tol/2) {
                mvprintw(y,x  ,"%c",'%'); mvprintw(y,x+1,"%c",'%'); }
                break;
             case 17: if (fabs((alpha+B-C)*beta + (alpha+C-B)*gamma - alpha*(alpha+B+C)) < tol) {
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
