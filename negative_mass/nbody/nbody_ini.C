#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

typedef double  real;

const int ndim = 3;
const int nposp = 50;
const int nnegp = 100;
const real massp = 1.0;
const real massn = -0.0001;
const real span = 50;
const real vmax = 0.0001;

float random(int i){return (float) rand()/RAND_MAX * i}

int main()
{
  srand ( time(NULL) );

  int ntot = nposp + nnegp;

  cout << ntot << endl;
  cout << "0." << endl;

  for (int i = 1; i <= ntot; i++){
    real mass = massp;
    if (i > nposp) {mass = massn;}
    cout << mass << "  ";
    // randomized positions in x, y
    for (int j = 1; j<=ndim; j++){
      real sign = 1.0;
      if (random(2) > 1){ sign = -1.0;}
      cout << sign*random(1) * span << "  ";
      }
    if (ndim == 2) {cout << "0.0" << "  ";}
    // random velocities: vx, vy
    for (int j = 1; j<=ndim; j++){
      real sign = 1.0;
      if (random(2) > 1){ sign = -1.0;}
      cout << sign*random(1) * vmax << "  ";
    }
    if (ndim == 2) {cout << "0.0" << "  ";}
    cout << endl;
  }
}


