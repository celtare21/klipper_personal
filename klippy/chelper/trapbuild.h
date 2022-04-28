#ifndef TRAPBUILD_H
#define TRAPBUILD_H

#include "list.h"

struct qmove;

struct vtrap {
    struct list_head trapezoid;
    struct qmove *accel_head, *decel_head;
};

void vtrap_init(struct vtrap *vt);
struct qmove * vtrap_flush(struct vtrap *vt, struct list_node *next_pos
                           , double *end_v2);
struct qmove * vtrap_clear(struct vtrap *vt, struct list_node *next_pos);
void vtrap_add_as_accel(struct vtrap *vt, struct qmove *move);
void vtrap_add_as_decel(struct vtrap *vt, struct qmove *move);

#endif  // trapbuild.h
