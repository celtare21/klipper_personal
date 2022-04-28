# S-Curve move planning
#
# Copyright (C) 2020  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import logging
import chelper

class error(Exception):
    pass

# Class to track a list of pending move requests and to facilitate
# "look-ahead" across moves to combine acceleration between moves.
class AccelCombiningMoveQueue:
    def __init__(self, scurve, toolhead):
        self.scurve = scurve
        self.toolhead = toolhead
        old_queue = toolhead.get_move_queue()
        self.queue = []
        self._LOOKAHEAD_FLUSH_TIME = old_queue._LOOKAHEAD_FLUSH_TIME
        self.junction_flush = old_queue.junction_flush
        self.last_step_gen_time = self.toolhead.reactor.monotonic()
        ffi_main, ffi_lib = chelper.get_ffi()
        self.cqueue = ffi_main.gc(ffi_lib.moveq_alloc(), ffi_lib.free)
        self.moveq_add = ffi_lib.moveq_add
        self.moveq_plan = ffi_lib.moveq_plan
        self.moveq_getmove = ffi_lib.moveq_getmove
        self.moveq_reset = ffi_lib.moveq_reset
        self.cmove_accel_decel = ffi_main.gc(
                ffi_lib.move_accel_decel_alloc(), ffi_lib.free)
    def reset(self):
        del self.queue[:]
        self.moveq_reset(self.cqueue)
        self.junction_flush = self._LOOKAHEAD_FLUSH_TIME
        self.last_step_gen_time = self.toolhead.reactor.monotonic()
    def set_flush_time(self, flush_time):
        self.junction_flush = flush_time
    def get_last(self):
        if self.queue:
            return self.queue[-1]
        return None
    def is_empty(self):
        return not self.queue
    def flush(self, lazy=False):
        self.junction_flush = self._LOOKAHEAD_FLUSH_TIME
        queue = self.queue
        qsize = len(queue)
        start_moveq_plan = self.toolhead.reactor.monotonic()
        flush_count = self.moveq_plan(self.cqueue, lazy)
        end_moveq_plan = self.toolhead.reactor.monotonic()
        if flush_count < 0:
            raise error('Internal error in moveq_plan')
        elif not flush_count:
            return
        self._get_planned_moves(flush_count)
        logging.info("lazy = %s, qsize = %d, flush_count = %d, plan_time = %.6f"
                , lazy, qsize, flush_count, end_moveq_plan - start_moveq_plan);
        # Generate step times for all moves ready to be flushed
        start_process_moves = self.toolhead.reactor.monotonic()
        self.toolhead._process_moves(queue[:flush_count])
        end_process_moves = self.toolhead.reactor.monotonic()
        # Remove processed moves from the queue
        del queue[:flush_count]
        logging.info('moveq_flush: time between step gen = %.6f,'
                     ' process_moves time = %.6f',
                     end_process_moves - self.last_step_gen_time,
                     end_process_moves - start_process_moves)
        self.last_step_gen_time = self.toolhead.reactor.monotonic()
    def _get_planned_moves(self, flush_count):
        cmove_accel_decel = self.cmove_accel_decel
        for i in range(flush_count):
            move = self.queue[i]
            if self.moveq_getmove(self.cqueue, cmove_accel_decel) < 0:
                raise error('Internal error in moveq_getmove')
            move.total_accel_t = cmove_accel_decel.total_accel_t
            move.accel_t = cmove_accel_decel.accel_t
            move.accel_offset_t = cmove_accel_decel.accel_offset_t
            move.cruise_t = cmove_accel_decel.cruise_t
            move.total_decel_t = cmove_accel_decel.total_decel_t
            move.decel_t = cmove_accel_decel.decel_t
            move.decel_offset_t = cmove_accel_decel.decel_offset_t
            move.start_accel_v = cmove_accel_decel.start_accel_v
            move.cruise_v = cmove_accel_decel.cruise_v
            move.effective_accel = cmove_accel_decel.effective_accel
            move.effective_decel = cmove_accel_decel.effective_decel
    def add_move(self, move):
        scurve = self.scurve
        if self.queue:
            move.calc_junction(self.queue[-1])
        move.accel_order = scurve.accel_order
        self.queue.append(move)
        jerk = scurve.max_jerk if move.is_kinematic_move else 9999999999999999.9
        ret = self.moveq_add(
                self.cqueue, move.move_d,
                move.junction_max_v2, move.max_cruise_v2,
                scurve.accel_order, move.accel, move.accel_to_decel,
                jerk, scurve.min_jerk_limit_time)
        if ret:
            raise error('Internal error in moveq_add')
        self.junction_flush -= move.min_move_t
        if self.junction_flush <= 0.:
            # Enough moves have been queued to reach the target flush time.
            self.flush(lazy=True)

RINGING_REDUCTION_FACTOR = 10.

class SCurve:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.printer.register_event_handler("klippy:connect", self.connect)
        self.toolhead = None
        self.min_jerk_limit_time = config.getfloat(
                'min_jerk_limit_time', 0.02, above=0.)
        self.min_accel = config.getfloat('min_accel', None, above=0.)
        self.max_jerk = config.getfloat('max_jerk', None, above=0.)
        self.accel_order = config.getchoice(
                'acceleration_order', { "2": 2, "4": 4, "6": 6 }, "4")
        # Register gcode commands
        gcode = self.printer.lookup_object('gcode')
        gcode.register_command("SET_SCURVE",
                               self.cmd_SET_SCURVE,
                               desc=self.cmd_SET_SCURVE_help)
    def connect(self):
        self.toolhead = self.printer.lookup_object("toolhead")
        self._set_params(self.min_accel, self.max_jerk)
        # Inject a new move queue
        new_move_queue = AccelCombiningMoveQueue(self, self.toolhead)
        self.toolhead.replace_move_queue(new_move_queue)
    def _set_params(self, min_accel, max_jerk):
        mjlt = self.min_jerk_limit_time
        if not max_jerk:
            max_velocity, max_accel = self.toolhead.get_max_velocity()
            if not min_accel:
                min_accel = max_accel / RINGING_REDUCTION_FACTOR
            max_jerk = min_accel * 6. / mjlt
        else:
            min_accel = max_jerk * mjlt / 6.
        self.min_accel = min_accel
        self.max_jerk = max_jerk
    cmd_SET_SCURVE_help = "Set S-Curve parameters"
    def cmd_SET_SCURVE(self, gcmd):
        accel_order = gcmd.get_int(
            'ACCEL_ORDER', self.accel_order, minval=2, maxval=6)
        if accel_order not in [2, 4, 6]:
            raise gcmd.error(
                    "ACCEL_ORDER = %s is not a valid choice" % (accel_order,))
        self.accel_order = accel_order
        min_accel = gcmd.get_float('MIN_ACCEL', None, above=0.)
        max_jerk = gcmd.get_float('JERK', None, above=0.)
        if min_accel is not None or max_jerk is not None:
            self._set_params(min_accel, max_jerk)
        msg = ("accel_order: %d max_jerk: %.6f min_accel: %.6f"%(
                   accel_order, self.max_jerk, self.min_accel))
        self.printer.set_rollover_info("scurve", "scurve: %s" % (msg,))
        gcmd.respond_info(msg, log=False)

def load_config(config):
    return SCurve(config)
