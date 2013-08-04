/* ray.h: Class for ray with tracer
 *
 * Author: Gaute Hope <eg@gaute.vetsj.com> / 2013-08-04
 *
 */

# pragma once


# include <array>
# include <vector>
# include <utility>

# include <boost/numeric/odeint.hpp>
# include <boost/numeric/ublas/vector.hpp>

# include "rl.h"

# define isbad(x) (std::isnan(x) || std::isinf(x))

using namespace std;
using namespace boost::numeric;

// state type of system: position coordinates and slowness vector
// combined into one vector of size 6.
typedef ublas::vector<double> state_type; // size = 6

// make state_type is_resizable
// http://headmyshoulder.github.com/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/state_types__algebras_and_operations.html#boost_numeric_odeint.odeint_in_detail.state_types__algebras_and_operations.algebras_and_operations.boost_ublas
namespace boost { namespace numeric { namespace odeint {
  template <>
    struct is_resizeable <state_type>
    {
      typedef boost::true_type type;
      const static bool value = type::value;

    };

} } }

// ODE wrapper to use a member function as ODE system: {{{
// http://www.boost.org/doc/libs/1_53_0_beta1/libs/numeric/odeint/examples/bind_member_functions.cpp
template< class Obj , class Mem >
class ode_wrapper
{
  Obj m_obj;
  Mem m_mem;

  public:

  ode_wrapper( Obj obj , Mem mem ) : m_obj( obj ) , m_mem( mem ) { }

  template< class State , class Deriv , class Time >
    void operator()( const State &x , Deriv &dxdt , Time t )
    {
      (m_obj.*m_mem)( x , dxdt , t );
    }
};

template< class Obj , class Mem >
ode_wrapper< Obj , Mem > make_ode_wrapper( Obj obj , Mem mem )
{
  return ode_wrapper< Obj , Mem >( obj , mem );
}

// }}}

namespace oie {
  class Ray {
    public:
      Ray ();

      int interfacesn = 3;
      double interfaces[3] = { 100, 200, 300 };

      /* initial conditions */
      double theta;
      double phi;
      ublas::vector <double> ix; // initial position vector
      ublas::vector <double> ip; // initial slowness vector

      bool solved = false;

      // ODE solver
      void solve ();

      /*
      template <class StepperType>
        std::pair<double, state_type> find_intersect (StepperType &s, std::pair<double, double> r,
                                   state_type x0, state_type x1, Interface * i);

      */

      // ray path
      struct pb_s_t;
      vector <state_type> x;
      vector <double>     t;

      // ray equation
      void rayeq (const state_type &x, state_type &dxdt, const double t);

      // cross product
      template <typename VT> static
        ublas::vector <double> cross_prod (VT &v0, VT &v1);


      // wildish guess of steps generated by ode solver as initial guess
      // of state and time vectors to prevent un-necessary re-allocation
      # define ESTIMATE_STEPS 60

      // output
      static string str_vector (ublas::vector <double>);
  };
}

