/* rl.h: various other stuff..
 *
 * Author: Gaute Hope <eg@gaute.vetsj.com> / 2013-08-04
 *
 */

# pragma once

# ifdef GIT_DESC
  # define VERSION GIT_DESC
# else
  # define VERSION "[??]"
# endif

// helper for timing
# define utime_diff(t0,t1) ((double)t1.tms_utime - (double)t0.tms_utime)


