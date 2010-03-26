#define _GNU_SOURCE

#include <stdio.h>
#include <sched.h>
#include <syscall.h>
#include "pin_to_cpu.h"

void place_me_(int *cpu, int *cpu_rank)
{
  cpu_set_t mask;

    CPU_ZERO(&mask);
    CPU_SET(cpu_rank[*cpu], &mask);

    sched_setaffinity(syscall(SYS_gettid),
                      sizeof(cpu_set_t) / sizeof(__cpu_mask), &mask);

}
