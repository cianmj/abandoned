/*
 *  nbody.C:  an N-body integrator with a shared but variable time step
 *                (the same for all particles but changing in time), using
 *                the Hermite integration scheme.
 *                        
 *  External data format:
 *                      n
 *                      t
 *                      m1 r1_x r1_y r1_z v1_x v1_y v1_z
 *                      m2 r2_x r2_y r2_z v2_x v2_y v2_z
 *                      ...
 *                      mn rn_x rn_y rn_z vn_x vn_y vn_z
 *
 *_____________________________________________________________________________
 */

#include  <iostream>
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
#include  <unistd.h>                       // for getopt()
using namespace std;

typedef double  real;                      // "real" as a general name for the
                                           // standard floating-point data type

const int NDIM = 3;                        // number of spatial dimensions
const real boundary = 100;                 // boundary of periodic box

void correct_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt);
void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM],
            int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag);
void evolve_step(const real mass[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time);
void get_acc_jerk_pot_coll(const real mass[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time);
void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM], int n);
void predict_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt);
void put_snapshot(const real mass[], const real pos[][NDIM],
                  const real vel[][NDIM], int n, real t);
bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag);
void write_diagnostics(const real mass[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag);
void periodic_boundary(real pos[][NDIM], int n);

/*
 *
 */

int main(int argc, char *argv[])
{
    real  dt_param = 0.01;  // control parameter to determine time step size
    real  dt_dia = 1;         // time interval between diagnostics output
    real  dt_out = 1;         // time interval between output of snapshots
    real  dt_tot = 400;        // duration of the integration
    bool  init_out = false;    // if true: snapshot output with start at t = 0
                               //          with an echo of the input snapshot
    bool  x_flag = false;      // if true: extra debugging diagnostics output

    if (! read_options(argc, argv, dt_param, dt_dia, dt_out, dt_tot, init_out,
                       x_flag))
        return 1;                // halt criterion detected by read_options()

    int n;                       // N, number of particles in the N-body system
    cin >> n;

    real t;                      // time
    cin >> t;

    real * mass = new real[n];                  // masses for all particles
    real (* pos)[NDIM] = new real[n][NDIM];     // positions for all particles
    real (* vel)[NDIM] = new real[n][NDIM];     // velocities for all particles

    get_snapshot(mass, pos, vel, n);

    evolve(mass, pos, vel, n, t, dt_param, dt_dia, dt_out, dt_tot, init_out,
           x_flag);

    delete[] mass;
    delete[] pos;
    delete[] vel;
}

/*
 *
 */

bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag)
{
    int c;
    while ((c = getopt(argc, argv, "hd:e:o:t:ix")) != -1)
        switch(c){
            case 'h': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << endl;
                      return false;         // execution should stop after help
            case 'd': dt_param = atof(optarg);
                      break;
            case 'e': dt_dia = atof(optarg);
                      break;
            case 'i': i_flag = true;
                      break;
            case 'o': dt_out = atof(optarg);
                      break;
            case 't': dt_tot = atof(optarg);
                      break;
            case 'x': x_flag = true;
                      break;
            case '?': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << endl;
                      return false;        // execution should stop after error
            }

    return true;                         // ready to continue program execution
}

/*
 *
 */

void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM], int n)
{
    for (int i = 0; i < n ; i++){
        cin >> mass[i];                       // mass of particle i
        for (int k = 0; k < NDIM; k++)
            cin >> pos[i][k];                 // position of particle i
        for (int k = 0; k < NDIM; k++)
            cin >> vel[i][k];                 // velocity of particle i
    }
}
    
/*
 *
 */

void put_snapshot(const real mass[], const real pos[][NDIM],
                  const real vel[][NDIM], int n, real t)
{
    cout.precision(16);                       // full double precision

    cout << n << endl;                        // N, total particle number
    cout << t << endl;                        // current time
    int itag = 0;
    for (int i = 0; i < n ; i++){
      if (itag == 0 && mass[i] < 0){
	itag = 1;
        cout << endl;
        cout << endl;
      }
        cout << mass[i];                      // mass of particle i
        for (int k = 0; k < NDIM; k++)
            cout << ' ' << pos[i][k];         // position of particle i
        for (int k = 0; k < NDIM; k++)
            cout << ' ' << vel[i][k];         // velocity of particle i
        cout << endl;
    }
}
    
/*
 *
 */

void write_diagnostics(const real mass[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag)
{
    real ekin = 0;                       // kinetic energy of the n-body system
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            ekin += 0.5 * mass[i] * vel[i][k] * vel[i][k];

    real etot = ekin + epot;             // total energy of the n-body system

    if (init_flag)                       // at first pass, pass the initial
        einit = etot;                    // energy back to the calling function

    cerr << "at time t = " << t << " , after " << nsteps
         << " steps :\n  E_kin = " << ekin
         << " , E_pot = " << epot
         << " , E_tot = " << etot << endl;
    cerr << "                "
         << "absolute energy error: E_tot - E_init = "
         << etot - einit << endl;
    cerr << "                "
         << "relative energy error: (E_tot - E_init) / E_init = "
         << (etot - einit) / einit << endl;

    if (x_flag){
        cerr << "  for debugging purposes, here is the internal data "
             << "representation:\n";
        for (int i = 0; i < n ; i++){
            cerr << "    internal data for particle " << i+1 << " : " << endl;
            cerr << "      ";
            cerr << mass[i];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << pos[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << vel[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << acc[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << jerk[i][k];
            cerr << endl;
        }
    }
}
    
/*
 *  evolve/  --  integrates an N-body system, for a total duration dt_tot.
 *              Snapshots are sent to the standard output stream once every
 *              time interval dt_out.  Diagnostics are sent to the standard
 *              error stream once every time interval dt_dia.
 */

void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM],
            int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag)
{
    cerr << "Starting a Hermite integration for a " << n
         << "-body system,\n  from time t = " << t 
         << " with time step control parameter dt_param = " << dt_param
         << "  until time " << t + dt_tot 
         << " ,\n  with diagnostics output interval dt_dia = "
         << dt_dia << ",\n  and snapshot output interval dt_out = "
         << dt_out << "." << endl;

    real (* acc)[NDIM] = new real[n][NDIM];          // accelerations and jerks
    real (* jerk)[NDIM] = new real[n][NDIM];         // for all particles
    real epot;                        // potential energy of the n-body system
    real coll_time;                   // collision (close encounter) time scale

    get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time);

    /*    cout << mass[0] << " " << mass[1] << endl;
    cout << acc[0][0] << " " << acc[1][0] << endl;
    cout << jerk[0][0] << " " << jerk[1][0] << endl;
    cout << coll_time << endl;
    return;*/

    int nsteps = 0;               // number of integration time steps completed
    real einit;                   // initial total energy of the system

    write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps, einit,
                      true, x_flag);
    if (init_out)                                    // flag for initial output
        put_snapshot(mass, pos, vel, n, t);

    real t_dia = t + dt_dia;           // next time for diagnostics output
    real t_out = t + dt_out;           // next time for snapshot output
    real t_end = t + dt_tot;           // final time, to finish the integration

    while (true){
      // cerr << t << endl;
        while (t < t_dia && t < t_out && t < t_end){
	  real dt = dt_param * coll_time;
	  //	    real dt = 0.0001;
	  // cerr << dt << ' ' << coll_time << ' ' << endl;
            evolve_step(mass, pos, vel, acc, jerk, n, t, dt, epot, coll_time);
            nsteps++;
        }
        if (t >= t_dia){
            write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps,
                              einit, false, x_flag);
            t_dia += dt_dia;
        }
        if (t >= t_out){
            put_snapshot(mass, pos, vel, n, t);
            t_out += dt_out;
        }
        if (t >= t_end)
            break;
    }

    delete[] acc;
    delete[] jerk;
}

/*
 *  evolve_step  --  takes one integration step for an N-body system, using the
 *                   Hermite algorithm.
 */

void evolve_step(const real mass[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time)
{
    real (* old_pos)[NDIM] = new real[n][NDIM];
    real (* old_vel)[NDIM] = new real[n][NDIM];
    real (* old_acc)[NDIM] = new real[n][NDIM];
    real (* old_jerk)[NDIM] = new real[n][NDIM];

    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
          old_pos[i][k] = pos[i][k];
          old_vel[i][k] = vel[i][k];
          old_acc[i][k] = acc[i][k];
          old_jerk[i][k] = jerk[i][k];
        }

    predict_step(pos, vel, acc, jerk, n, dt);
    get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time);
    correct_step(pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk,
                 n, dt);
    t += dt;

    periodic_boundary(pos,n);

    delete[] old_pos;
    delete[] old_vel;
    delete[] old_acc;
    delete[] old_jerk;
}

/*
 *  predict_step  --  takes the first approximation of one Hermite integration
 *                    step.
 */

void predict_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2
                                      + jerk[i][k]*dt*dt*dt/6;
            vel[i][k] += acc[i][k]*dt + jerk[i][k]*dt*dt/2;
        }
}

/*
 *  correct_step  --  takes one iteration to improve the new values of position
 *                    and velocities
 */

void correct_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
                                      + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
            pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
                                      + (old_acc[i][k] - acc[i][k])*dt*dt/12;
        }
}

/*
 *  get_acc_jerk_pot_coll  --  calculates accelerations and jerks, and as side
 *                             effects also calculates potential energy and
 *                             the time scale coll_time for significant changes
 *                             in local configurations to occur.
 */

void get_acc_jerk_pot_coll(const real mass[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            acc[i][k] = jerk[i][k] = 0;
    epot = 0;
    const real VERY_LARGE_NUMBER = 1e300;
    real coll_time_q = VERY_LARGE_NUMBER;      // collision time to 4th power
    real coll_est_q;                           // collision time scale estimate
                                               // to 4th power (quartic)
    for (int i = 0; i < n ; i++){
        for (int j = i+1; j < n ; j++){            // rji[] is the vector from
            real rji[NDIM];                        // particle i to particle j
            real vji[NDIM];                        // vji[] = d rji[] / d t
            for (int k = 0; k < NDIM ; k++){
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k];
            }
            real r2 = 0;                           // | rji |^2
            real v2 = 0;                           // | vji |^2
            real rv_r2 = 0;                        // ( rij . vij ) / | rji |^2
            for (int k = 0; k < NDIM ; k++){
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k];
            }
            rv_r2 /= r2;
            real r = sqrt(r2);                     // | rji |
            real r3 = r * r2;                      // | rji |^3

// add the {i,j} contribution to the total potential energy for the system:

            epot -= mass[i] * mass[j] / r;

// add the {j (i)} contribution to the {i (j)} values of acceleration and jerk:

            real da[3];                            // main terms in pairwise
            real dj[3];                            // acceleration and jerk
            for (int k = 0; k < NDIM ; k++){
                da[k] = rji[k] / r3;                           // see equations
                dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;    // in the header
            }
            for (int k = 0; k < NDIM ; k++){
                acc[i][k] += mass[j] * da[k];                 // using symmetry
                acc[j][k] -= mass[i] * da[k];                 // find pairwise
                jerk[i][k] += mass[j] * dj[k];                // acceleration
                jerk[j][k] -= mass[i] * dj[k];                // and jerk
            }

// first collision time estimate, based on unaccelerated linear motion:

            coll_est_q = (r2*r2) / (v2*v2);
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;

// second collision time estimate, based on free fall:

            real da2 = 0;                                  // da2 becomes the 
            for (int k = 0; k < NDIM ; k++)                // square of the 
                da2 += da[k] * da[k];                      // pair-wise accel-
            double mij = abs(mass[i]) + abs(mass[j]);      // eration between
            da2 *= mij * mij;                              // particles i and j

            coll_est_q = r2/da2;
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;
        }                                     
    }                                               // from q for quartic back
    coll_time = sqrt(sqrt(coll_time_q));            // to linear collision time
}                                             


/*
 * Implement periodic boundary conditions
 *
 */
void periodic_boundary(real pos[][NDIM], int n)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
	  if (abs(pos[i][k]) > boundary)
	    pos[i][k] = -pos[i][k];
        }
}

