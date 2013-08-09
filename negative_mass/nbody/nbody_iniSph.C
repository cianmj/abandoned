#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

typedef double real;

const int ndim = 3;
const int nposp = 200;
const int nnegp = 600;
const real massp = 1.0;
const real massn = -0.0001;
const real span = 50;
const real radius_pos = 20;
const real vmax = 0.0001;

float random(int i){return (float) rand()/RAND_MAX * i}

int main()
{
  srand ( time(NULL) );

  int ntot = nposp + nnegp;
  real pi = 4.*atan(1.);

  cout << ntot << endl;
  cout << "0." << endl;

  for (int i = 1; i <= ntot; i++){
    real mass = massp;
    real radius = (random(1) * radius_pos;
    if (i > nposp) {
      mass = massn; 
      radius = radius_pos + random(1) * (span-radius_pos);
    }
    cout << mass << "  ";
    // randomized positions in x, y, z
    real theta = random(1) * pi;
    real phi = random(1) * 2.*pi;
    cout << radius*sin(theta)*cos(phi) << "  "
	 << radius*sin(theta)*sin(phi) << "  "
	 << radius*cos(theta) << " ";
     // random velocities: vx, vy
    for (int j = 1; j<=ndim; j++){
      real sign = 1.0;
      if (random(2) > 1){ sign = -1.0;}
      cout << random(1) * vmax << "  ";
    }
    if (ndim == 2) {cout << "0.0" << "  ";}
    cout << endl;
  }
}
