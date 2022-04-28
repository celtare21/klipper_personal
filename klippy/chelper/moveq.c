// Look-ahead movement planning queue
//
// Copyright (C) 2016-2019  Kevin O'Connor <kevin@koconnor.net>
// Copyright (C) 2019  Dmitry Butyugin <dmbutyugin@google.com>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <assert.h> // assert
#include <math.h> // sqrt
#include <stddef.h> // offsetof
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "accelgroup.h"
#include "compiler.h" // __visible
#include "itersolve.h"
#include "moveq.h"
#include "pyhelper.h" // get_monotonic
#include "scurve.h" // scurve
#include "trapq.h" // trap_accel_decel
#include "trapbuild.h"

#define MOVE_DEBUG 0

static const double EPSILON = 0.000000001;

static struct qmove *
qmove_alloc(void)
{
    struct qmove *m = malloc(sizeof(*m));
    memset(m, 0, sizeof(*m));
    return m;
}

struct moveq * __visible
moveq_alloc(void)
{
    struct moveq *mq = malloc(sizeof(*mq));
    memset(mq, 0, sizeof(*mq));
    list_init(&mq->moves);
    init_combiner(&mq->accel_combiner);
    return mq;
}

// Allocate a new 'move_accel_decel' object
struct move_accel_decel * __visible
move_accel_decel_alloc(void)
{
    struct move_accel_decel *accel_decel = malloc(sizeof(*accel_decel));
    memset(accel_decel, 0, sizeof(*accel_decel));
    return accel_decel;
}


void __visible
moveq_reset(struct moveq *mq)
{
    struct qmove *m = NULL, *nm = NULL;
    list_for_each_entry_safe(m, nm, &mq->moves, node) {
        list_del(&m->node);
        free(m);
    }
    reset_combiner(&mq->accel_combiner);
    memset(mq, 0, sizeof(*mq));
    list_init(&mq->moves);
    init_combiner(&mq->accel_combiner);
}

static struct qmove *
backward_smoothed_pass(struct moveq *mq, int lazy, int *ret)
{
    *ret = 0;
    int update_flush_limit = lazy;
    // Traverse queue from last to first move and determine maximum
    // junction speed assuming the robot comes to a complete stop
    // after the last move.
    struct list_head delayed;
    list_init(&delayed);
    double next_smoothed_v2 = 0., peak_cruise_v2 = 0.;
    reset_junctions(&mq->accel_combiner, 0.);
    struct qmove *move = NULL, *pm = NULL, *flush_limit = NULL;
    list_for_each_entry_reversed_safe(move, pm, &mq->moves, node) {
        // Determine peak cruise velocity
        double reachable_smoothed_v2 = next_smoothed_v2 + move->smooth_delta_v2;
        double smoothed_v2 = MIN(move->max_smoothed_v2, reachable_smoothed_v2);
        if (smoothed_v2 < reachable_smoothed_v2) {
            // It's possible for this move to accelerate
            if (smoothed_v2 + move->smooth_delta_v2 > next_smoothed_v2
                    || !list_empty(&delayed)) {
                // This move can decelerate or this is a full accel
                // move after a full decel move
                if (update_flush_limit && peak_cruise_v2) {
                    flush_limit = move;
                    update_flush_limit = 0;
                }
                peak_cruise_v2 = (smoothed_v2 + reachable_smoothed_v2) * .5;
                peak_cruise_v2 = MIN(move->max_cruise_v2, peak_cruise_v2);
            }
            struct qmove *m = NULL;
            if (!update_flush_limit && move != flush_limit) {
                move->max_cruise_v2 = MIN(move->max_cruise_v2
                        , peak_cruise_v2);
                move->junction_max_v2 = MIN(move->junction_max_v2
                        , peak_cruise_v2);
                list_for_each_entry(m, &delayed, node) {
                    m->max_cruise_v2 = MIN(m->max_cruise_v2, peak_cruise_v2);
                    m->junction_max_v2 = MIN(m->junction_max_v2
                            , peak_cruise_v2);
                }
                m = list_next_entry(move, node);
                if (lazy && list_at_end(m, &mq->moves, node)) {
                    errorf("Logic error: smoothed peak velocity trapezoid at "
                            "the end of the move queue");
                    *ret = ERROR_RET;
                    return NULL;
                }
                if (!list_at_end(m, &mq->moves, node)) {
                    m->junction_max_v2 = MIN(m->junction_max_v2
                            , peak_cruise_v2);
                }
            }
            struct qmove *nm = NULL, *qm = move;
            // Put delayed moves back to their places in mq->moves list.
            list_for_each_entry_safe(m, nm, &delayed, node) {
                list_del(&m->node);
                list_add_after(&m->node, &qm->node);
                qm = m;
            }
        } else {
            // Delay calculating this move until peak_cruise_v2 is known
            list_del(&move->node);
            list_add_head(&move->node, &delayed);
        }
        if (mq->smoothed_pass_limit == move)
            break;
        next_smoothed_v2 = smoothed_v2;
    }
    if (!list_empty(&delayed)) {
        errorf("Non-empty 'delayed' queue after the smoothed pass");
        *ret = ERROR_RET;
        return NULL;
    }
    mq->smoothed_pass_limit = flush_limit;
    if (update_flush_limit) return NULL;
    return flush_limit;
}

static void
backward_pass(struct moveq *mq, struct qmove *end)
{
    // Backward pass for full acceleration
    double junction_max_v2 = 0.;
    reset_junctions(&mq->accel_combiner, junction_max_v2);
    struct qmove *move = (end == NULL
            ? list_last_entry(&mq->moves, struct qmove, node)
            : list_prev_entry(end, node));
    for (; !list_at_end(move, &mq->moves, node);
            move = list_prev_entry(move, node)) {
        // Restore the default accel and decel values if they were modified
        // on previous backward pass
        move->decel_group = move->accel_group = move->default_accel;

        struct accel_group *decel = &move->decel_group;
        process_next_accel(&mq->accel_combiner, move, decel, junction_max_v2);
        junction_max_v2 = move->junction_max_v2;
    }
}

static struct qmove *
compute_safe_flush_limit(struct moveq *mq, int lazy, struct qmove *end)
{
    if (!lazy) return end;
    // Last move has unknown end junction_max_v2
    double junction_max_v2 = 1e100;
    reset_junctions(&mq->accel_combiner, junction_max_v2);

    struct qmove *move = (end == NULL)
        ? list_last_entry(&mq->moves, struct qmove, node)
        : list_prev_entry(end, node);
    struct qmove *start = list_first_entry(&mq->moves, struct qmove, node);
    struct qmove *flush_limit = NULL;
    // Go over all moves from the end of the queue. Find all moves in the queue
    // that have a sufficiently distant junction point on their deceleration
    // paths with junction_max_v2 reached, and keep references to such points in
    // fallback_decel in case they are needed in a forward pass.
    for (; !list_at_end(move, &mq->moves, node) && move != start;
            move = list_prev_entry(move, node)) {
        // Check if the current 'move' has a junction point on its deceleration
        // path where junction_max_v2 of that point is reached. Such junction
        // point must be farther than the minimum safe distance. If such point
        // does not exists, this means that this move can change its max_end_v2
        // if more moves are added to the queue later, so it is not safe to
        // flush the queue until this move.
        if (process_fallback_decel(&mq->accel_combiner, move, junction_max_v2)
                && !flush_limit)
            flush_limit = list_next_entry(move, node);
        junction_max_v2 = move->junction_max_v2;
    }
    return flush_limit ? flush_limit : start;
}

static struct qmove *
forward_pass(struct moveq *mq, struct qmove *end, int lazy)
{
    struct qmove *move = list_first_entry(&mq->moves, struct qmove, node);
    if (move == end)
        return NULL;
    double start_v2 = mq->prev_end_v2;
    double max_end_v2 = move->decel_group.max_end_v2;
    if (max_end_v2 + EPSILON < start_v2) {
        double decel_start_v2 = MIN(move->fallback_decel.max_start_v2
                , move->fallback_decel.start_accel->max_start_v2);
        errorf("Warning: impossible to reach the committed v2 = %.3f"
                ", max velocity = %.3lf, fallback to suboptimal planning"
                ", decelerate to %.3lf", start_v2, max_end_v2, decel_start_v2);
        struct accel_group *decel = &move->decel_group;
        *decel = move->fallback_decel;
        decel->max_end_v2 = start_v2;
        set_max_start_v2(decel->start_accel, MIN(start_v2, decel_start_v2));
    }
    struct vtrap vt;
    vtrap_init(&vt);
    reset_junctions(&mq->accel_combiner, start_v2);
    double prev_cruise_v2 = start_v2;
    struct qmove *last_flushed_move = NULL, *next_move = NULL;
    for (; !list_at_end(move, &mq->moves, node) && move != end;
            move = next_move) {
        // Track next_move early because move will be moved to trapezoid list
        next_move = list_next_entry(move, node);
        struct accel_group *accel = &move->accel_group;
        struct accel_group *decel = &move->decel_group;

        process_next_accel(&mq->accel_combiner, move, accel
                , MIN(move->junction_max_v2, prev_cruise_v2));
        int can_accelerate = decel->max_end_v2 > accel->max_start_v2 + EPSILON;
        if (can_accelerate) {
            // This move can accelerate
            if (vt.decel_head) {
                // Complete the previously combined velocity trapezoid
                if (move->fallback_decel.move || !lazy)
                    // It is safe to flush this velocity trapezoid
                    last_flushed_move = vtrap_flush(&vt, &move->node
                            , &mq->prev_end_v2);
                else
                    // It might be unsafe to flush this velocity trapezoid, so
                    // do not update last_flushed_move, but wait if subsequent
                    // velocity trapezoid will have fallback_decel computed.
                    vtrap_flush(&vt, &move->node, NULL);
            }
            vtrap_add_as_accel(&vt, move);
        }
        int must_decelerate = accel->max_end_v2 + EPSILON > decel->max_start_v2;
        if (must_decelerate || !can_accelerate) {
            // This move must decelerate after acceleration,
            // or this is a full decel move after full accel move.
            for (; move != end; move = next_move
                    , next_move = list_next_entry(move, node)) {
                vtrap_add_as_decel(&vt, move);
                if (move == decel->start_accel->move) break;
            }
            if (move == end) break;
            reset_junctions(&mq->accel_combiner
                    , decel->start_accel->max_start_v2);
        }
        prev_cruise_v2 = move->max_cruise_v2;
    }
    if (!lazy) {
        struct list_node *next_pos = end == NULL ? &mq->moves.root : &end->node;
        if (vt.decel_head || vt.accel_head)
            last_flushed_move = vtrap_flush(&vt, next_pos, &mq->prev_end_v2);
    } else {
        assert(end != NULL);
        // Just put the remaining moves back to mq->moves
        vtrap_clear(&vt, &end->node);
    }
    return last_flushed_move;
}

int __visible
moveq_add(struct moveq *mq, double move_d
          , double junction_max_v2, double max_cruise_v2
          , int accel_order, double accel, double smoothed_accel
          , double jerk, double min_jerk_limit_time)
{
    struct qmove *m = qmove_alloc();
    m->move_d = move_d;
    fill_accel_group(&m->default_accel, m, accel_order, accel, jerk
            , min_jerk_limit_time);
    m->max_cruise_v2 = max_cruise_v2;
    m->junction_max_v2 = junction_max_v2;
    m->smooth_delta_v2 = 2. * smoothed_accel * move_d;

    if (!list_empty(&mq->moves)) {
        struct qmove *prev_move = list_last_entry(
                &mq->moves, struct qmove, node);
        m->max_smoothed_v2 =
            prev_move->max_smoothed_v2 + prev_move->smooth_delta_v2;
        m->max_smoothed_v2 = MIN(
                MIN(m->max_smoothed_v2, junction_max_v2),
                MIN(m->max_cruise_v2, prev_move->max_cruise_v2));
    }
    list_add_tail(&m->node, &mq->moves);
    return 0;
}

int __visible
moveq_getmove(struct moveq *mq, struct move_accel_decel *accel_decel)
{
    memset(accel_decel, 0, sizeof(*accel_decel));
    if (list_empty(&mq->moves)) {
        errorf("Move queue is empty");
        return ERROR_RET;
    }
    struct qmove *move = list_first_entry(&mq->moves, struct qmove, node);
    struct accel_group *accel = &move->accel_group;
    struct accel_group *decel = &move->decel_group;
    accel_decel->accel_order = accel->accel_order;
    // Determine move velocities
    accel_decel->start_accel_v = accel->start_accel_v;
    accel_decel->cruise_v = move->cruise_v;
    // Determine the effective accel and decel
    accel_decel->effective_accel = accel->effective_accel;
    accel_decel->effective_decel = decel->effective_accel;
    // Determine time spent in each portion of move (time is the
    // distance divided by average velocity)
    accel_decel->accel_t = accel->accel_t;
    accel_decel->accel_offset_t = accel->accel_offset_t;
    accel_decel->total_accel_t = accel->total_accel_t;
    accel_decel->decel_t = decel->accel_t;
    accel_decel->decel_offset_t = decel->accel_offset_t;
    accel_decel->total_decel_t = decel->total_accel_t;
    double cruise_d = move->move_d - accel->accel_d - decel->accel_d;
    accel_decel->cruise_t = cruise_d / move->cruise_v;
    // Only used to track smootheness, can be deleted
    struct scurve s_acc, s_dec;
    scurve_fill(&s_acc, accel_decel->accel_order, accel_decel->accel_t,
            accel_decel->accel_offset_t, accel_decel->total_accel_t,
            accel_decel->start_accel_v, accel_decel->effective_accel);
    scurve_fill(&s_dec, accel_decel->accel_order, accel_decel->decel_t,
            accel_decel->decel_offset_t, accel_decel->total_decel_t,
            accel_decel->cruise_v, -accel_decel->effective_decel);
    double start_v, end_v;
    if (accel_decel->accel_t > EPSILON)
        start_v = scurve_velocity(&s_acc, 0.);
    else if (accel_decel->cruise_t > EPSILON)
        start_v = accel_decel->cruise_v;
    else
        start_v = scurve_velocity(&s_dec, 0.);
    if (accel_decel->decel_t > EPSILON)
        end_v = scurve_velocity(&s_dec, accel_decel->decel_t);
    else if (accel_decel->cruise_t > EPSILON)
        end_v = accel_decel->cruise_v;
    else
        end_v = scurve_velocity(&s_acc, accel_decel->accel_t);
    #if MOVE_DEBUG
    static volatile long long move_idx = 0;
    errorf("Move [%lld] ms_v2=%.3f, mc_v2=%.3f, me_v2=%.3f with"
            " move_d=%.6f, max_c_v2=%.3f, jct_v2=%.3f, accel=%.3f, decel=%.3f"
            ", accel_t=%.6f, cruise_t=%.6f, decel_t=%.6f", move_idx
            , start_v*start_v, move->cruise_v*move->cruise_v, end_v*end_v
            , move->move_d, move->max_cruise_v2, move->junction_max_v2
            , accel_decel->effective_accel, accel_decel->effective_decel
            , accel_decel->accel_t, accel_decel->cruise_t
            , accel_decel->decel_t);
    ++move_idx;
    #endif
    if (accel_decel->cruise_t < -EPSILON) {
        errorf("Logic error: impossible move ms_v=%.3f, mc_v=%.3f"
                ", me_v=%.3f, accel_d = %.3f, decel_d = %.3f"
                " with move_d=%.3lf, accel=%.3f, decel=%.3f"
                ", jerk=%.3f", start_v, move->cruise_v, end_v
                    , accel->accel_d, decel->accel_d, move->move_d
                    , accel->max_accel, decel->max_accel
                    , accel->max_jerk);
        return ERROR_RET;
    }
    accel_decel->cruise_t = MAX(0., accel_decel->cruise_t);
    if (fabs(mq->prev_move_end_v - start_v) > 0.0001) {
        errorf("Logic error: velocity jump from %.6f to %.6f"
                , mq->prev_move_end_v, start_v);
        return ERROR_RET;
    }
    // Remove processed move from the queue
    list_del(&move->node);
    free(move);
    mq->prev_move_end_v = end_v;
    return 0;
}

int __visible
moveq_plan(struct moveq *mq, int lazy)
{
    if (list_empty(&mq->moves))
        return 0;
    int ret;
    struct qmove *flush_limit = backward_smoothed_pass(mq, lazy, &ret);
    if (ret) return ret;
    if (lazy && !flush_limit)
        return 0;
    backward_pass(mq, flush_limit);
    struct qmove *safe_flush_limit = compute_safe_flush_limit(
            mq, lazy, flush_limit);
    struct qmove *last_flushed_move = forward_pass(mq, safe_flush_limit, lazy);
    if (!last_flushed_move)
        return 0;
    int flush_count = 0;
    struct qmove *move = NULL;
    list_for_each_entry(move, &mq->moves, node) {
        ++flush_count;
        if (move == last_flushed_move) break;
    }
    return flush_count;
}
