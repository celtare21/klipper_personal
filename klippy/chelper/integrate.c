// Helpers to integrate the smoothing weight function.
//
// Copyright (C) 2019-2020  Kevin O'Connor <kevin@koconnor.net>
// Copyright (C) 2020  Dmitry Butyugin <dmbutyugin@google.com>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include "integrate.h"
#include "scurve.h"

#include <stdlib.h> // malloc
#include <string.h> // memset

// Calculate the inverse of the norm of the weight function
static double
calc_inv_norm(double hst)
{
    // Inverse norm of the weight function ((t-T)^2-h^2)^2
    return 15. / (16. * hst * hst * hst * hst * hst);
}

// Calculate (t^2-h^2)^2
static inline double
w(const struct smoother *sm, double t)
{
    double t2 = t*t;
    double v = t2 - sm->h2;
    return v*v;
}

static const double w_antideriv_coeffs[][3] = {
    {1./5., -2./3., 1./1.,},
    {1./6., -2./4., 1./2.,},
    {1./7., -2./5., 1./3.,},
    {1./8., -2./6., 1./4.,},
    {1./9., -2./7., 1./5.,},
    {1./10., -2./8., 1./6.,},
    {1./11., -2./9., 1./7.,},
};

// Integrate t^n * (t^2-h^2)^2
static inline double
iwtn(const struct smoother *sm, int n, double t)
{
    const double *coeffs = w_antideriv_coeffs[n];
    double t2 = t*t;
    double v = (coeffs[0] * t2 + coeffs[1] * sm->h2) * t2 + coeffs[2] * sm->h4;
    for (; n >= 0; --n)
        v *= t;
    return v;
}

// Integrate scurve s(t) with smoothing weight function
// ((t-T)^2-h^2)^2 over the range [start; end] with T == -toff
double
integrate_weighted(const struct smoother *sm, double pos, struct scurve *s
                   , double start, double end, double toff)
{
    double toff2 = toff * toff;
    double v = toff2 - sm->h2;
    // Calculate s(t) * w(t) integral as either expansion of s(t) or w(t)
    // over powers of t. w(t) expansion becomes numerically unstable when
    // abs(toff) >> hst, and s(t) - when abs(toff) >> total_accel_t.
    // Note that when abs(toff) >> hst, it means that abs(toff) ~ move_t, so
    // it is not possible that abs(toff) >> total_accel_t at the same time.
    if (toff2 > sm->h2) {
        pos += scurve_eval(s, -toff);
        scurve_offset(s, -toff);

        start += toff; end += toff;
        double res = s->c6 * (iwtn(sm, 6, end) - iwtn(sm, 6, start));
        res += s->c5 * (iwtn(sm, 5, end) - iwtn(sm, 5, start));
        res += s->c4 * (iwtn(sm, 4, end) - iwtn(sm, 4, start));
        res += s->c3 * (iwtn(sm, 3, end) - iwtn(sm, 3, start));
        res += s->c2 * (iwtn(sm, 2, end) - iwtn(sm, 2, start));
        res += s->c1 * (iwtn(sm, 1, end) - iwtn(sm, 1, start));
        res += pos * (iwtn(sm, 0, end) - iwtn(sm, 0, start));
        return res;
    } else {
        double res = (scurve_tn_antiderivative(s, 4, end)
                      - scurve_tn_antiderivative(s, 4, start));
        res += 4. * toff * (scurve_tn_antiderivative(s, 3, end)
                            - scurve_tn_antiderivative(s, 3, start));
        res += 2. * (3. * toff2 - sm->h2) *
            (scurve_tn_antiderivative(s, 2, end)
             - scurve_tn_antiderivative(s, 2, start));
        res += 4. * toff * v * (scurve_tn_antiderivative(s, 1, end)
                                - scurve_tn_antiderivative(s, 1, start));
        res += v * v * (scurve_tn_antiderivative(s, 0, end)
                        - scurve_tn_antiderivative(s, 0, start));
        start += toff; end += toff;
        res += pos * (iwtn(sm, 0, end) - iwtn(sm, 0, start));
        return res;
    }
}

// Integrate velocity jumps near the ends of the range [start; end] with
// smoothing weight function. To get correct results it is required to sum up
// the returned values over the full integration range [T-hst; T+hst].
double
integrate_velocity_jumps(const struct smoother *sm, const struct scurve *s
                         , double start, double end, double toff)
{
    double start_v = scurve_velocity(s, start);
    double end_v = scurve_velocity(s, end);
    // Velocity jumps integration assumes that the weight function vanishes at
    // the integration bounds T-hst and T+hst to ignore velocity jumps there.
    return start_v * w(sm, start + toff) - end_v * w(sm, end + toff);
}

struct smoother *
alloc_smoother(double hst) {
    struct smoother *sm = malloc(sizeof(*sm));
    memset(sm, 0, sizeof(*sm));
    sm->hst = hst;
    sm->inv_norm = calc_inv_norm(hst);
    sm->h2 = hst * hst;
    sm->h4 = sm->h2 * sm->h2;
    return sm;
}
