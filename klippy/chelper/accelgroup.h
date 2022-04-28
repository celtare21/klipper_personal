#ifndef ACCELGROUP_H
#define ACCELGROUP_H

struct qmove;

// A group of moves accelerating (or decelerating) together
struct accel_group {
    int accel_order;
    double max_accel, min_accel;
    double max_jerk;
    double min_jerk_limit_time;
    double combined_d, accel_d;
    double accel_t, accel_offset_t, total_accel_t;
    double start_accel_v;
    double effective_accel;
    struct accel_group *start_accel, *next_accel;
    struct qmove *move;
    double max_start_v, max_start_v2, max_end_v2;
};

void fill_accel_group(struct accel_group *ag, struct qmove *m, int accel_order
                      , double accel, double jerk, double min_jerk_limit_time);
void limit_accel(struct accel_group *ag, double accel, double jerk);
void set_max_start_v2(struct accel_group *ag, double start_v2);
double calc_max_v2(const struct accel_group* ag);
double calc_effective_accel(const struct accel_group *ag, double cruise_v);
double calc_min_accel_time(const struct accel_group *ag, double cruise_v);
double calc_min_accel_dist(const struct accel_group *ag, double cruise_v);
double calc_min_accel_group_time(const struct accel_group *ag, double cruise_v);
double calc_max_safe_v2(const struct accel_group *ag);

#endif  // accelgroup.h
