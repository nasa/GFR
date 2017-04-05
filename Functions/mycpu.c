#include <utmpx.h>
#include <unistd.h>
int sched_getcpu();
long int sysconf(int inputvalue);

int findmycpu_ ()
{
    int cpu;
    cpu = sched_getcpu();
    return cpu;
}

int get_num_cpus_per_host_ ()
{
    long int nprocs = -1;
    int num_cpus = -1;
#ifdef _SC_NPROCESSORS_ONLN
    nprocs = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    if (nprocs >= 1)
    {
        num_cpus = nprocs;
    }
    return num_cpus;
}


int get_max_cpus_per_host_ ()
{
    long int nprocs_max = -1;
    int max_cpus = -1;
#ifdef _SC_NPROCESSORS_CONF
    nprocs_max = sysconf(_SC_NPROCESSORS_CONF);
#endif
    if (nprocs_max >= 1)
    {
        max_cpus = nprocs_max;
    }
    return max_cpus;
}
