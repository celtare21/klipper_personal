#ifndef ACCELCOMBINE_H
#define ACCELCOMBINE_H

#include "accelgroup.h"
#include "list.h"

struct accel_group;
struct qmove;

struct junction_point {
    struct list_node node;
    // Combined acceleration limits that must be respected
    // from this junction point
    struct accel_group accel;
    struct accel_group *move_ag;
    double min_start_time, min_end_time;
    double max_cruise_end_v2;
};

struct accel_combiner {
    struct list_head junctions;
    struct junction_point *prev_best_jp;
    double junct_start_v2;
};

void init_combiner(struct accel_combiner *ac);
void reset_combiner(struct accel_combiner *ac);
void reset_junctions(struct accel_combiner *ac, double start_v2);
void process_next_accel(struct accel_combiner *ac, struct qmove *move
                        , struct accel_group *ag, double junction_max_v2);
int process_fallback_decel(struct accel_combiner *ac, struct qmove *move
                           , double next_junction_max_v2);

#endif  // accelcombine.h
