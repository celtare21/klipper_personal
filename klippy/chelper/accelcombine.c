// Moves combiner computing groups of moves that should accelerate together
//
// Copyright (C) 2019  Dmitry Butyugin <dmbutyugin@google.com>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <assert.h> // assert
#include <math.h> // sqrt
#include <stddef.h> // offsetof
#include <string.h> // memset
#include "accelcombine.h"
#include "accelgroup.h"
#include "compiler.h" // likely
#include "moveq.h" // qmove

static const double EPSILON = 0.000000001;

void
reset_junctions(struct accel_combiner *ac, double start_v2)
{
    while (likely(!list_empty(&ac->junctions))) {
        struct junction_point *jp = list_first_entry(
                &ac->junctions, struct junction_point, node);
        list_del(&jp->node);
    }
    ac->junct_start_v2 = start_v2;
    ac->prev_best_jp = NULL;
}

static struct junction_point *
init_junction_point(struct accel_combiner *ac, struct qmove *move
                    , struct accel_group *ag, double junction_max_v2)
{
    struct junction_point *new_jp = &move->jp;
    memset(new_jp, 0, sizeof(*new_jp));
    new_jp->accel = *ag;
    new_jp->accel.start_accel = &new_jp->accel;
    new_jp->move_ag = ag;
    struct junction_point *prev_jp = ac->prev_best_jp;
    double start_v2;
    if (likely(prev_jp)) {
        double prev_end_v2 = MIN(prev_jp->accel.max_end_v2
                , prev_jp->max_cruise_end_v2);
        start_v2 = MIN(junction_max_v2, prev_end_v2);
        new_jp->min_start_time = prev_jp->min_end_time;
    } else
        start_v2 = MIN(junction_max_v2, ac->junct_start_v2);
    set_max_start_v2(&new_jp->accel, start_v2);
    return new_jp;
}

inline static int
check_can_combine(struct accel_combiner *ac
                  , const struct accel_group *next_accel)
{
    if (unlikely(list_empty(&ac->junctions)))
        return 0;
    struct junction_point *prev_jp = list_last_entry(
            &ac->junctions, struct junction_point, node);
    struct accel_group *prev_accel = &prev_jp->accel;
    return prev_accel->accel_order == next_accel->accel_order;
}

static void
drop_decelerating_jps(struct accel_combiner *ac, double accel_limit_v2)
{
    while (likely(!list_empty(&ac->junctions))) {
        struct junction_point *last_jp = list_last_entry(
                &ac->junctions, struct junction_point, node);
        if (unlikely(last_jp->accel.max_start_v2 < accel_limit_v2 + EPSILON))
            // First point from which deceleration is not required
            return;
        // This point must decelerate
        list_del(&last_jp->node);
    }
}

static void
drop_nonaccelerating_jps(struct accel_combiner *ac, double accel_limit_v2)
{
    drop_decelerating_jps(ac, accel_limit_v2 - 2. * EPSILON);
}

static void
limit_accel_jps(struct accel_combiner *ac, struct accel_group *ag
                , double junction_max_v2)
{
    struct junction_point *jp = NULL;
    list_for_each_entry(jp, &ac->junctions, node) {
        // Make sure to not exceed junction_max_v2 during acceleration
        double junction_accel_limit = 0.5 * (junction_max_v2
                - jp->accel.max_start_v2) / jp->accel.combined_d;
        limit_accel(&jp->accel, MIN(junction_accel_limit, ag->max_accel)
                , ag->max_jerk);
    }
}

inline static double
calc_min_accel_end_time(struct junction_point *jp, double cruise_v2)
{
    return jp->min_start_time
        + calc_min_accel_group_time(&jp->accel, sqrt(cruise_v2));
}

static struct junction_point *
calc_best_jp(struct accel_combiner *ac, struct qmove *move
             , struct accel_group *new_ag)
{
    double max_cruise_v2 = move->max_cruise_v2;
    struct junction_point *jp = NULL, *best_jp = NULL;
    list_for_each_entry(jp, &ac->junctions, node) {
        // Choose the best acceleration option
        jp->accel.combined_d += move->move_d;
        jp->accel.max_end_v2 = calc_max_v2(&jp->accel);
        jp->max_cruise_end_v2 = max_cruise_v2;
        jp->min_end_time = calc_min_accel_end_time(
                jp, MIN(jp->accel.max_end_v2, max_cruise_v2));
        if (!best_jp || best_jp->min_end_time > jp->min_end_time + EPSILON)
            best_jp = jp;
    }
    assert(best_jp);
    return best_jp;
}

void
process_next_accel(struct accel_combiner *ac, struct qmove *move
                   , struct accel_group *ag, double junction_max_v2)
{
    struct junction_point *new_jp = init_junction_point(
            ac, move, ag, junction_max_v2);
    double start_v2 = new_jp->accel.max_start_v2;
    if (unlikely(!check_can_combine(ac, ag)))
        reset_junctions(ac, start_v2);

    drop_nonaccelerating_jps(ac, MIN(start_v2, junction_max_v2));
    limit_accel_jps(ac, ag, junction_max_v2);

    // Add the current move to the list (with combined_d == 0)
    list_add_tail(&new_jp->node, &ac->junctions);
    struct junction_point *best_jp = calc_best_jp(ac, move, ag);
    ac->prev_best_jp = best_jp;

    limit_accel(ag, best_jp->accel.max_accel, best_jp->accel.max_jerk);
    set_max_start_v2(ag, start_v2);
    ag->max_end_v2 = best_jp->accel.max_end_v2;
    ag->combined_d = best_jp->accel.combined_d;
    // Point to the real accel_group instance.
    ag->start_accel = best_jp->move_ag;
}

static void
maybe_add_new_fallback_decel_jp(struct accel_combiner *ac, struct qmove *move
                                , double next_junction_max_v2)
{
    double start_v2 = move->decel_group.max_start_v2;
    if (next_junction_max_v2 > start_v2 + EPSILON)
        // This junction point does not reach next_junction_max_v2
        return;
    if (!list_empty(&ac->junctions)) {
        struct junction_point *last_jp = list_last_entry(
                &ac->junctions, struct junction_point, node);
        if (next_junction_max_v2 < last_jp->accel.max_start_v2 + EPSILON)
            // Last junction point already covers this next_junction_max_v2
            return;
    }
    struct junction_point *new_jp = &move->jp;
    memset(new_jp, 0, sizeof(*new_jp));
    new_jp->accel = move->default_accel;
    new_jp->accel.start_accel = &new_jp->accel;
    new_jp->move_ag = &move->decel_group;
    set_max_start_v2(&new_jp->accel, next_junction_max_v2);

    // Add the current move to the list (with combined_d == 0)
    list_add_tail(&new_jp->node, &ac->junctions);
}

static int
find_fallback_decel(struct accel_combiner *ac, struct qmove *move
                    , double max_end_v2)
{
    struct junction_point *jp = NULL;
    struct accel_group *fallback_decel = &move->fallback_decel;
    list_for_each_entry_reversed(jp, &ac->junctions, node) {
        jp->accel.combined_d += move->move_d;
        double safe_end_v2 = calc_max_safe_v2(&jp->accel);
        if (max_end_v2 <= safe_end_v2 + EPSILON) {
            *fallback_decel = jp->accel;
            fallback_decel->max_end_v2 = safe_end_v2;
            set_max_start_v2(
                    fallback_decel, jp->accel.start_accel->max_start_v2);
            // Point to the real accel_group instance.
            fallback_decel->start_accel = jp->move_ag;
            fallback_decel->move = move;
            return 1;
        }
    }
    memset(fallback_decel, 0, sizeof(*fallback_decel));
    return 0;
}

int
process_fallback_decel(struct accel_combiner *ac, struct qmove *move
                       , double next_junction_max_v2)
{
    if (unlikely(!check_can_combine(ac, &move->default_accel)))
        reset_combiner(ac);

    struct accel_group *decel = &move->decel_group;
    double start_v2 = decel->max_start_v2;
    double max_end_v2 = MIN(decel->max_end_v2, move->junction_max_v2);

    if (next_junction_max_v2 > EPSILON)
        drop_decelerating_jps(ac, MIN(start_v2, next_junction_max_v2));
    else
        reset_junctions(ac, 0.);
    limit_accel_jps(ac, &move->default_accel, next_junction_max_v2);
    maybe_add_new_fallback_decel_jp(ac, move, next_junction_max_v2);

    return find_fallback_decel(ac, move, max_end_v2);
}

void
init_combiner(struct accel_combiner *ac)
{
    memset(ac, 0, sizeof(*ac));
    list_init(&ac->junctions);
}

void
reset_combiner(struct accel_combiner *ac)
{
    reset_junctions(ac, 0.);
}
