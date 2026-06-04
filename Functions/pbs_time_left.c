#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define __USE_GNU
#include <dlfcn.h>

static size_t save_time = -1;
static size_t time_left = -1;

static int unset_envar(char *env, char *libv)
{
   char *ep, *es, *ee;

   /* get current value */
   ep = getenv(env);
   if (!ep) return 0;

   /* delete libv if exists */
   es = strstr(ep, libv);
   if (!es) return 0;
   ee = es + strlen(libv);
   while (es != ep && *es != ':' && *es != ' ') es--;
   if (es == ep && *ee) ee++;
   strcpy(es, ee);

   return 1;
}

FILE *mypopen(const char *command, const char *type)
{
   FILE *ret, *(*real_popen)();
   char *es, *ev;
   char *lpd="LD_PRELOAD";

   dlerror();
   real_popen = (FILE *(*)())dlsym(RTLD_NEXT, "popen");
   if (!real_popen){
      fprintf(stderr, "Symbol popen not defined\n");
      exit(1);
   }
   ev = getenv(lpd);
   if (ev) {
      es = (char *)malloc(strlen(ev)+1);
      if (!es) {
         perror("malloc");
         exit(1);
      }
      strcpy(es, ev);
      unsetenv(lpd);
   }

   ret = real_popen(command, type);
   if (ev) {
      setenv(lpd, es, 1);
      free(es);
   }
   return ret;
}

void pbs_time_left_(const char *command, long *time_ret) {

  FILE *fpo;
  char *jobid, *qstat_wait_time;
  char cmd[1024],time_str[1024], short_str[10];
  char line[1024], *env_ret;
  size_t time_now, wait_time;

  int ht,mt,st,hu,mu,su;

  if ( (jobid = getenv("PBS_JOBID")) == NULL ) {
    *time_ret = 0;
    return;
  }

  if ( (qstat_wait_time = getenv("PBS_WAIT_TIME")) == NULL ) {
     wait_time = 600;
  } else {
     wait_time = atoi(qstat_wait_time);
  }

  time_now = time(&time_now);

  if ( save_time != -1 && time_now - save_time < wait_time ) {
    *time_ret = time_left - (time_now - save_time);
    return;
  }

//sprintf(cmd,"/PBS/bin/qstat -f %s\n",jobid);
  sprintf(cmd,"%s -f %s\n",command,jobid);

  if( fpo = (FILE *)mypopen(cmd,"r") ) {
    
    while( fgets(line, 1024, fpo) != NULL ) {
      if ( strcasestr(line,"used.walltime =") != NULL ) {
	strncpy(time_str,strstr(line," = "),1024);
        if ( !strncmp(&time_str[6],":",1) ) {
	  strncpy(short_str,&time_str[3],3); short_str[3] = '\0'; hu = atoi(short_str);
	  strncpy(short_str,&time_str[7],2); short_str[2] = '\0'; mu = atoi(short_str);
	  strncpy(short_str,&time_str[10],2); short_str[2] = '\0'; su = atoi(short_str);
        } else {
	  strncpy(short_str,&time_str[3],2); short_str[2] = '\0'; hu = atoi(short_str);
	  strncpy(short_str,&time_str[6],2); short_str[2] = '\0'; mu = atoi(short_str);
	  strncpy(short_str,&time_str[9],2); short_str[2] = '\0'; su = atoi(short_str);
        }
      }

      if ( strcasestr(line,"List.walltime =") != NULL ) {
	strncpy(time_str,strstr(line," = "),1024);
        if ( !strncmp(&time_str[6],":",1) ) {
	  strncpy(short_str,&time_str[3],3); short_str[3] = '\0'; ht = atoi(short_str);
	  strncpy(short_str,&time_str[7],2); short_str[2] = '\0'; mt = atoi(short_str);
	  strncpy(short_str,&time_str[10],2); short_str[2] = '\0'; st = atoi(short_str);
        } else {
	  strncpy(short_str,&time_str[3],2); short_str[2] = '\0'; ht = atoi(short_str);
	  strncpy(short_str,&time_str[6],2); short_str[2] = '\0'; mt = atoi(short_str);
	  strncpy(short_str,&time_str[9],2); short_str[2] = '\0'; st = atoi(short_str);
        }
      }

    }

    pclose(fpo);

  }
  
  printf("From pbs_time_left.c: \n");
  printf("    Resource_List.walltime = %3d:%02d:%02d \n",ht,mt,st);
  printf("   resources_used.walltime = %3d:%02d:%02d \n",hu,mu,su);
  *time_ret = (long)3600 * (long)(ht-hu) + (long)60 * (long)(mt-mu) + (long)(st-su);
  time_left = *time_ret;
  save_time = time(&save_time);
  return;

}

