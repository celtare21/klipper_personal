// Jerk-limiting computation of acceleration limits for groups of moves
//
// Copyright (C) 2019  Dmitry Butyugin <dmbutyugin@google.com>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // sqrt
#include "accelgroup.h"
#include "compiler.h" // unlikely

void
fill_accel_group(struct accel_group *ag, struct qmove *m, int accel_order
                 , double accel, double jerk, double min_jerk_limit_time)
{
    ag->accel_order = accel_order;
    ag->max_accel = accel;
    ag->max_jerk = jerk;
    ag->min_jerk_limit_time = min_jerk_limit_time;
    ag->min_accel = jerk * min_jerk_limit_time / 6.;
    if (unlikely(ag->min_accel > ag->max_accel)) {
        ag->min_accel = ag->max_accel;
    }
    ag->move = m;
}

void
limit_accel(struct accel_group *ag, double accel, double jerk)
{
    if (unlikely(accel < 0.)) accel = 0.;
    ag->max_accel = MIN(ag->max_accel, accel);
    ag->max_jerk = MIN(ag->max_jerk, jerk);
    double min_accel = ag->max_jerk * ag->min_jerk_limit_time / 6.;
    if (unlikely(ag->min_accel > min_accel)) {
        ag->min_accel = min_accel;
    }
    if (unlikely(ag->min_accel > ag->max_accel)) {
        ag->min_accel = ag->max_accel;
    }
}

void
set_max_start_v2(struct accel_group *ag, double start_v2)
{
    ag->max_start_v2 = start_v2;
    ag->max_start_v = sqrt(start_v2);
}

double
calc_max_v2(const struct accel_group* ag)
{
    double dist = ag->combined_d;
    // Check if accel is the limiting factor
    double start_v2 = ag->start_accel->max_start_v2;
    double max_accel_v2 = start_v2 + 2.0 * dist * ag->max_accel;
    // Compute maximum achievable speed with limited kinematic jerk using
    // max(jerk) == 6 * accel / accel_time, which is exact for accel order 4
    // and is quite accurate for accel order 6:
    // max(jerk) == 10 / sqrt(3) * accel / accel_time ~=
    //     5.774 * accel / accel_time
    // This leads to the cubic equation
    // (max_v^2 - start_v^2) * (max_v + start_v) / 2 ==
    //     dist^2 * jerk / 3
    // which is solved using Cardano's formula.
    double start_v = ag->start_accel->max_start_v;
    double a = 2./3. * start_v;
    double b = a*a*a;
    double c = dist * dist * ag->max_jerk / 3.;
    double d = sqrt(c * (c + 2. * b));
    double e = pow(b + c + d, 1./3.);
    if (unlikely(e < 0.000000001))
        return start_v;
    double max_v = e + a*a / e - start_v / 3.;
    double max_v2 = max_v * max_v;
    if (unlikely(max_accel_v2 < max_v2))
        max_v2 = max_accel_v2;
    double min_accel_v2 = start_v2 + 2.0 * dist * ag->min_accel;
    if (unlikely(min_accel_v2 > max_v2))
        max_v2 = min_accel_v2;
    return max_v2;
}

inline double
calc_effective_accel(const struct accel_group *ag, double cruise_v)
{
    double effective_accel = sqrt(ag->max_jerk
            * (cruise_v - ag->start_accel->max_start_v) / 6.);
    if (unlikely(effective_accel > ag->max_accel))
        effective_accel = ag->max_accel;
    if (unlikely(effective_accel < ag->min_accel))
        effective_accel = ag->min_accel;
    return effective_accel;
}

inline double
calc_min_accel_time(const struct accel_group *ag, double cruise_v)
{
    double delta_v = cruise_v - ag->start_accel->max_start_v;
    if (fabs(delta_v) < 0.000000001)
        return 0.;
    double min_accel_time = delta_v / ag->max_accel;
    double accel_t = sqrt(6. * delta_v / ag->max_jerk);
    if (likely(accel_t > min_accel_time))
        min_accel_time = accel_t;
    if (likely(ag->min_accel)) {
        double accel_t = delta_v / ag->min_accel;
        if (unlikely(accel_t < min_accel_time))
            min_accel_time = accel_t;
    }
    return min_accel_time;
}

inline double
calc_min_accel_dist(const struct accel_group *ag, double cruise_v)
{
    double start_v = ag->start_accel->max_start_v;
    if (unlikely(cruise_v <= start_v))
        return 0.;
    double accel_t = calc_min_accel_time(ag, cruise_v);
    return (start_v + cruise_v) * 0.5 * accel_t;
}

inline double
calc_max_safe_v2(const struct accel_group *ag)
{
    double dist = ag->combined_d;
    double start_v2 = ag->start_accel->max_start_v2;
    double max_v2 = 2. * ag->max_accel * dist + start_v2;
    // It is possible to accelerate from any velocity to this one over the
    // accumulated distance dist.
    double v2 = pow((9./16.) * dist * dist * ag->max_jerk, (2./3.));
    // Such min v2 is achieved when accelerating from v2 / 9 velocity.
    // But if start_v2 is smaller than v2 / 9, it is sufficient to
    // consider the worst-case acceleration from start_v2 only.
    if (unlikely(start_v2 * 9. < v2))
        v2 = calc_max_v2(ag);
    max_v2 = MIN(max_v2, v2);
    return max_v2;
}

inline double
calc_min_accel_group_time(const struct accel_group *ag, double cruise_v)
{
    if (ag->start_accel->max_start_v >= cruise_v)
        // No acceleration possible - just cruising
        return ag->combined_d / cruise_v;

    double start_v = ag->start_accel->max_start_v;
    double accel_t = calc_min_accel_time(ag, cruise_v);
    double accel_d = (start_v + cruise_v) * 0.5 * accel_t;
    double cruise_t = (ag->combined_d - accel_d) / cruise_v;

    return accel_t + cruise_t;
}
