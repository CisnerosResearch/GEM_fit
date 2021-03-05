#ifndef _TIMER_H
#define _TIMER_H

#define MAX_TIMERS 256

#ifdef __STDC__
#include <time.h>

typedef struct {
    void *parent;
    void *child;
    clock_t sum_utime;
    clock_t start_utime;
} calltimer_t;

/* PROTOTYPES */
void start_timer(void *parent, void *child);
void stop_timer(void *parent, void *child);
calltimer_t *find_timer(void *parent, void *child);

#endif /* __STDC__ */

#endif /* _TIMER_H */
