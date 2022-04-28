// Velocity trapezoid builder
//
// Copyright (C) 2019  Dmitry Butyugin <dmbutyugin@google.com>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // sqrt
#include <stddef.h> // offsetof
#include <string.h> // memset
#include "compiler.h" // MIN
#include "list.h"
#include "moveq.h" // qmove
#include "scurve.h" // scurve_get_time
#include "trapbuild.h"

static double
calc_move_peak_v2(struct qmove *m)
{
    struct accel_group *accel = &m->accel_group;
    struct accel_group *decel = &m->decel_group;
    double total_d = accel->combined_d + decel->combined_d - m->move_d;
    double high_v = sqrt(MAX(accel->max_end_v2, decel->max_end_v2));
    double low_v = 0.;
    while (likely(high_v - low_v > 0.000000001)) {
        double guess_v = (high_v + low_v) * 0.5;
        double accel_d = calc_min_accel_dist(accel, guess_v);
        double decel_d = calc_min_accel_dist(decel, guess_v);
        if (accel_d <= accel->combined_d && decel_d <= decel->combined_d &&
                accel_d + decel_d <= total_d)
            low_v = guess_v;
        else
            high_v = guess_v;
    }
    return low_v * low_v;
}

static double
calc_trap_peak_v2(struct qmove *accel_head, struct qmove *decel_head)
{
    if (decel_head != accel_head) {
        if (!decel_head)
            return accel_head->accel_group.max_end_v2;
        double peak_v2 = MIN(decel_head->decel_group.max_end_v2
                , decel_head->junction_max_v2);
        if (accel_head)
            peak_v2 = MIN(peak_v2, accel_head->accel_group.max_end_v2);
        return peak_v2;
    }
    double peak_v2 = calc_move_peak_v2(decel_head);
    return MIN(peak_v2, decel_head->max_cruise_v2);
}

static void
set_accel(struct accel_group* combined, double cruise_v2
          , int time_offset_from_start)
{
    if (combined->start_accel->max_start_v2 > cruise_v2)
        set_max_start_v2(combined->start_accel, cruise_v2);
    double start_accel_v = combined->start_accel->max_start_v;
    double cruise_v = sqrt(cruise_v2);
    double avg_v = (cruise_v + start_accel_v) * 0.5;
    double combined_accel_t = calc_min_accel_time(combined, cruise_v);
    double combined_accel_d = avg_v * combined_accel_t;
    double effective_accel = calc_effective_accel(combined, cruise_v);
    struct scurve s;
    scurve_fill(&s, combined->accel_order,
            combined_accel_t, 0., combined_accel_t,
            start_accel_v, effective_accel);
    double remaining_accel_t = combined_accel_t;
    double remaining_accel_d = combined_accel_d;
    struct accel_group *a = combined->start_accel;
    for (;;) {
        a->move->cruise_v = cruise_v;
        if (likely(remaining_accel_d > 0)) {
            a->effective_accel = effective_accel;
            a->total_accel_t = combined_accel_t;
            a->accel_d = MIN(a->move->move_d, remaining_accel_d);
            a->start_accel_v = start_accel_v;
            double next_pos = a->accel_d + combined_accel_d - remaining_accel_d;
            if (time_offset_from_start) {
                a->accel_offset_t = combined_accel_t - remaining_accel_t;
                a->accel_t = scurve_get_time(&s, next_pos) - a->accel_offset_t;
            } else {
                a->accel_offset_t =
                    combined_accel_t - scurve_get_time(&s, next_pos);
                a->accel_t = remaining_accel_t - a->accel_offset_t;
            }
            remaining_accel_t -= a->accel_t;
            remaining_accel_d -= a->move->move_d;
        }
        if (unlikely(a == combined)) break;
        a = a->next_accel;
    }
}

static void
set_trap_decel(struct qmove *decel_head, struct list_head *trapezoid
               , double cruise_v2, double *end_v2)
{
    struct qmove *m = decel_head;
    while (!list_at_end(m, trapezoid, node)) {
        set_accel(&m->decel_group, cruise_v2, /*time_offset_from_start=*/0);
        m = m->decel_group.start_accel->move;
        cruise_v2 = MIN(cruise_v2, m->decel_group.max_start_v2);
        if (end_v2 != NULL)
            *end_v2 = m->decel_group.max_start_v2;
        m = list_next_entry(m, node);
    }
}

static void
set_trap_accel(struct qmove *accel_head, struct list_head *trapezoid
               , double cruise_v2)
{
    struct qmove *m = accel_head;
    while (!list_at_end(m, trapezoid, node)) {
        set_accel(&m->accel_group, cruise_v2, /*time_offset_from_start=*/1);
        m = m->accel_group.start_accel->move;
        cruise_v2 = MIN(cruise_v2, m->accel_group.max_start_v2);
        m = list_prev_entry(m, node);
    }
}

struct qmove *
vtrap_flush(struct vtrap *vt, struct list_node *next_pos, double *end_v2)
{
    double peak_cruise_v2 = calc_trap_peak_v2(vt->accel_head, vt->decel_head);
    if (vt->decel_head)
        set_trap_decel(vt->decel_head, &vt->trapezoid, peak_cruise_v2, end_v2);
    else if (end_v2 != NULL)
        *end_v2 = peak_cruise_v2;
    if (vt->accel_head)
        set_trap_accel(vt->accel_head, &vt->trapezoid, peak_cruise_v2);
    return vtrap_clear(vt, next_pos);
}

void
vtrap_add_as_accel(struct vtrap *vt, struct qmove *move)
{
    if (vt->accel_head) {
        struct qmove *prev_move =
            list_last_entry(&vt->trapezoid, struct qmove, node);
        prev_move->accel_group.next_accel = &move->accel_group;
    }
    list_del(&move->node);
    list_add_tail(&move->node, &vt->trapezoid);
    vt->accel_head = move;
}

void
vtrap_add_as_decel(struct vtrap *vt, struct qmove *move)
{
    if (vt->decel_head) {
        struct qmove *next_move =
            list_last_entry(&vt->trapezoid, struct qmove, node);
        move->decel_group.next_accel = &next_move->decel_group;
    }
    if (!vt->decel_head)
        vt->decel_head = move;
    if (vt->accel_head != move) {
        list_del(&move->node);
        list_add_tail(&move->node, &vt->trapezoid);
    }
}


void
vtrap_init(struct vtrap *vt)
{
    list_init(&vt->trapezoid);
    vt->accel_head = NULL;
    vt->decel_head = NULL;
}

struct qmove *
vtrap_clear(struct vtrap *vt, struct list_node *next_pos)
{
    struct qmove *move = NULL, *next = NULL, *prev = NULL;
    list_for_each_entry_safe(move, next, &vt->trapezoid, node) {
        move->accel_group.next_accel = move->decel_group.next_accel = NULL;
        list_del(&move->node);
        list_add_before(&move->node, next_pos);
        prev = move;
    }
    vt->accel_head = vt->decel_head = NULL;
    return prev;
}
