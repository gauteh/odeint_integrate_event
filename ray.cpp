/* ray.cpp: Trace ray with events 
 *
 * Author: Gaute Hope <eg@gaute.vetsj.com> / 2013-08-04
 *
 */

# include "ray.h"

# include <thread>
# include <chrono>
# include <stdexcept>

# include <sys/times.h>

# include <vector>
# include <array>
# include <string>
# include <iostream>
# include <utility>

# include <boost/numeric/odeint.hpp>
# include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
# include <boost/numeric/ublas/vector.hpp>
# include <boost/math/special_functions.hpp>

using namespace boost::numeric;
using namespace boost::numeric::odeint;
using namespace boost::numeric::odeint::detail;
using namespace boost::math;
using namespace std;
using namespace oie;

int main () {
  cout << "Creating test ray and solving.." << endl;
  Ray * r = new Ray ();

  r->solve ();

  cout << "Finished." << endl;
}

namespace oie {
  Ray::Ray () {
    ix.resize (3);
    ip.resize (3);

    ix[0] = 100;
    ix[1] = 100;
    ix[2] = 100;

    theta = 0;
    phi   = 0;

    // calculate slowness vector
    double veli = 3000;

    ip[0] = cos (phi) * sin (theta) / veli;
    ip[1] = sin (phi) * sin (theta) / veli;
    ip[2] = cos (theta)             / veli;
  }

  /* observer (of ray path) method structure (pass x and t) */
  struct Ray::pb_s_t {
    vector <state_type> & m_states;
    vector <double> & m_times;

    pb_s_t (vector<state_type> &states, vector<double> &times)
      : m_states (states), m_times (times) { }

    void operator() (const state_type &x, double t) {
      m_states.push_back (x);
      m_times.push_back (t);

      // print out
      cout << "t = " << t << ", x = " << str_vector (x) << endl;
    }
  };

  /* ray equation: step */
  void Ray::rayeq (const state_type &x, state_type &dxdt, const double t) {

    ublas::vector_range <state_type> dx (dxdt, ublas::range (0, 3));
    ublas::vector_range <state_type> dp (dxdt, ublas::range (3, 6));

    //cout << "step t = " << t << " -------------" << endl;
    ublas::vector <double> p(3);
    p[0] = x[3];
    p[1] = x[4];
    p[2] = x[5];

    // find velocity at x
    //double veli = vel->interp_vel (x[0], x[1], x[2], 1.6);
    double gradient = 0.001 * 3000;
    double veli = 3000 + gradient * x[2]; // linear increasing velocity
    //cout << "veli = " << veli << endl;

    // find velocity gradient at x
    //auto vg = vel->interp_grad (x[0], x[1], x[2], 0.0);
    ublas::vector <double> vg(3);
    vg[0] = 0;
    vg[1] = 0;
    vg[2] = gradient;
    //cout << "vg = " << str_vector (vg) << endl;

    //cout << "x = " << str_vector (x) << endl;
    //cout << "p = " << str_vector (p) << endl;

    // ray equation:
    // x is position
    // p is slowness
    //
    // equations from: Keers et. al., 1997
    //
    // dx = v(x)^2 * p
    // dp = - (1 / v(x)) * grad(v(x))

    dx = veli * veli    * p;
    dp = (-1.0 / veli)  * vg;

    cout << "ray: step t = " << t << ", dxdt = " << str_vector (dxdt) << endl;
  }

  void Ray::solve () { // {{{
    //cout << "w" << workerid << ": ray: solving, theta: " << theta << ", phi: " << phi << endl;
    //cout << "w" << workerid << ": ix = " << str_vector(ix) << ", ip = " << str_vector (ip) << endl;

    // result of ode stepper will be stored in x and t, guessing size
    // to reserve and avoid re-allocation of memory too often.
    x.reserve (ESTIMATE_STEPS);
    t.reserve (ESTIMATE_STEPS);

    // initial conditions of ray: creating a state_type from the
    // initial position and slowness vector.
    state_type xx (6);
    xx[0] = ix[0];
    xx[1] = ix[1];
    xx[2] = ix[2];
    xx[3] = ip[0];
    xx[4] = ip[1];
    xx[5] = ip[2];

    // solve ray:
    // - stop at out of bounds
    // - detect interface hits and:
    //   * continue transmitted ray
    //   * create   reflected ray

    // set up ode stepper algorithm (equivialent of ode45 in MATLAB)
    // this is suitable for a stiff system ( which I hope this is.. ),
    // and is used with an adaptive integrator. the controller parameters
    // specifies the error tolerances.
    typedef runge_kutta_dopri5 <state_type> stepper_type;
    //typedef euler < state_type > stepper_type; // for a very basic stepper

    // create controller for stepper, monitors error tolerance for
    // the adaptive integrator: respectively absolute and relative error.
    typedef boost::numeric::odeint::result_of::make_dense_output<
      runge_kutta_dopri5< state_type > >::type dense_stepper_type;

    //auto controlled_stepper = make_controlled <stepper_type> (1e-4, 1e-4);
    dense_stepper_type stepper = make_dense_output <stepper_type> (1e-5, 1e-5);
    //auto stepper = stepper_type (); // for a un-controlled stepper

    /*
    // test of un-controlled stepper
    state_type next (6);

    for (int i = 0; i < 3; i++) {
      double dt = 0.01;
      stepper.do_step (rl::Ray::rayeq, xx, dt * i, next, dt);

      cout << "next = " << str_vector (next) << endl;
      xx = next;
    }
    */

    // for testing stepper
    // stepper.do_step (rl::Ray::rayeq, xx, 0.0, 0.01);

    // start timing ode solver
    struct tms t0, t1;
    times (&t0);

    // create wrapper for member function as system
    auto rayeq_w = make_ode_wrapper (*this, &oie::Ray::rayeq);

    // solve ode: detect out of bounds and interface hits
    // based on integrate_adaptive for dense stepper
    state_type dxdt (6);

    const int max_steps = 1000;
    int       steps = 0;
    controlled_step_result res = fail;

    // time parameters
    double tt      = 0.0; // t start
    double t_end   = 10;
    double dt      = t_end - tt; // initial step size

    auto obs = pb_s_t (x, t);

    bool terminal = false;
    stepper.initialize (xx, tt, dt);

    double     prevt = tt;
    state_type prev (xx);
    double     curt  = 0.0;
    state_type cur (6);

    while (!terminal &&
           less_with_sign (stepper.current_time (),
           t_end,
           stepper.current_time_step ())
          && steps < max_steps)
    {
      while (!terminal &&
             less_eq_with_sign (stepper.current_time () + stepper.current_time_step (),
             t_end,
             stepper.current_time_step ()))
      {
        obs ( stepper.current_state (), stepper.current_time () );
        auto r = stepper.do_step (rayeq_w);
        steps++;

        // event detection:
        curt = stepper.current_time ();
        cur  = stepper.current_state ();

        // EVENT detection:

        // check for interface intersects (only from interface 1 and up, handle
        // surface specially)

        // surface
        cout << "checking if hit surface.." << endl;
        if (cur[2] <= 0.0) {
          cout << "ray: Event: hit surface" << endl;
          terminal    = true;

          // Determine intersect with desired accuracy:
          //auto ei = find_surface_intersect <dense_stepper_type> (stepper, r, prev);

          // Store event
          //obs ( ei.second, ei.first );

          continue;
        }

        // check whether we hit any general interfaces
        for (int i = 0; i < interfacesn; i++) {

          double intf = interfaces[i];
          cout << "checking interface i = " << i << ", at = " << intf << endl;

          if (intf > prev[2] && intf < cur[2] ||
              intf < prev[2] && intf > cur[2]) {
            cout << "ray: Event: hit interface (" << intf << ") between: " << str_vector (prev) << " and " << str_vector (cur) << endl;

            // Determine intersect with desired accuracy
            //auto ei = find_intersect <dense_stepper_type> (stepper, r, prev, cur, intf);

            // Handle event

            // Store event
            //obs ( ei.second, ei.first );

            if (terminal)
              continue;
          }
        }

        // check for out of bounds
        if (cur[0] < 0 || cur[0] >= (1000) ||
            cur[1] < 0 || cur[1] >= (1000) ||
            cur[2] < 0 || cur[2] >= (1000))
        {
          cout << "ray: event: state out of bounds." << endl;
          cout << str_vector (cur) << endl;
          terminal = true;
          continue;
        }

        // check for bad values
        if (isbad (stepper.current_time ()) ||
            isbad (cur[0]) || isbad (cur[1]) || isbad (cur[2]) ||
            isbad (cur[3]) || isbad (cur[4]) || isbad (cur[5]))
        {
          cout << "ray: event: bad in time or state." << endl;
          terminal = true;
          continue;
        }

        prevt = curt;
        prev  = cur;
      }

      // re-initialize with smaller time step (prev is still correct state and time)
      if (!terminal)
        stepper.initialize (stepper.current_state (), stepper.current_time (), t_end - stepper.current_time ());
    }

    // include last step
    if (!terminal) {
      obs ( stepper.current_state (), stepper.current_time () );
    }

    // stop time
    times (&t1);

    //cout << "w" << workerid << ": ray: done: steps: " << steps << ", utime: " << utime_diff(t0,t1) << " ticks, ";

    // shrink x and t to actual, final, size
    x.shrink_to_fit ();
    t.shrink_to_fit ();

    //cout << "x size: " << x.size () << ", t size: " << t.size () << endl;

    solved = true;
  } // }}}

  /* Utility functions: str_vector, write_ray {{{ */
  string Ray::str_vector (ublas::vector <double> v) {
    ostringstream ss;
    ss << "[ ";

    for (auto i = v.begin (); i < v.end (); i++) {
      ss << *i;
      if (i != v.end () -1)
        ss << ", ";
    }

    ss << " ]";

    return ss.str ();
  }

  // cross product
  template <typename VT>
    ublas::vector <double> Ray::cross_prod (VT &b, VT &c) {
      // http://en.wikipedia.org/wiki/Cross_product#Mnemonic
      ublas::vector <double> a(3);

      a[0] = b[1] * c[2] - b[2] * c[1];
      a[1] = b[2] * c[0] - b[0] * c[2];
      a[2] = b[0] * c[1] - b[1] * c[0];

      return a;
    }
}

