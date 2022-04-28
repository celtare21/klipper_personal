#ifndef INTEGRATE_H
#define INTEGRATE_H

struct scurve;

struct smoother {
    double inv_norm, hst;
    double h2, h4;
};

struct smoother * alloc_smoother(double hst);
double integrate_weighted(const struct smoother *sm,
                          double pos, struct scurve *s,
                          double start, double end, double toff);
double integrate_velocity_jumps(const struct smoother *sm,
                                const struct scurve *s,
                                double start, double end, double toff);

#endif // integrate.h
