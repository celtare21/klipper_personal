#ifndef MOVEQ_H
#define MOVEQ_H

#include "accelcombine.h"
#include "accelgroup.h"
#include "itersolve.h"
#include "list.h"

struct trap_accel_decel;

struct qmove {
    struct list_node node;

    double cruise_v;
    double move_d;

    struct accel_group accel_group, decel_group, fallback_decel, default_accel;
    double smooth_delta_v2, max_smoothed_v2;
    double max_cruise_v2, junction_max_v2;

    struct junction_point jp;

    // Only used to track smootheness, can be deleted
    double start_v, end_v;
};

struct move_accel_decel {
    double accel_t, accel_offset_t, total_accel_t;
    double cruise_t;
    double decel_t, decel_offset_t, total_decel_t;
    double start_accel_v, cruise_v;
    double effective_accel, effective_decel;
    int accel_order;
};

struct moveq {
    double prev_end_v2;
    struct list_head moves;
    struct accel_combiner accel_combiner;
    struct qmove *smoothed_pass_limit;
    double prev_move_end_v;
};

struct move_accel_decel *move_accel_decel_alloc(void);

struct moveq *moveq_alloc(void);
void moveq_reset(struct moveq *mq);

int moveq_add(struct moveq *mq, double move_d
              , double junction_max_v2, double max_cruise_v2
              , int accel_order, double accel, double smoothed_accel
              , double jerk, double min_jerk_limit_time);
int moveq_plan(struct moveq *mq, int lazy);
int moveq_getmove(struct moveq *mq, struct move_accel_decel *accel_decel);

#endif  // moveq.h
