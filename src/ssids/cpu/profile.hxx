/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */
#pragma once

//#define PROFILE

#ifdef PROFILE
extern "C" {
#include <GTG.h>
}
#include "time.h"
#include <omp.h>
#endif

namespace spral { namespace ssids { namespace cpu {

class Profile {
public:
   class Task {
   public:
      Task(char const* name, int thread)
      : name(name), thread(thread), t1(Profile::now())
      {}

      void done() {
#ifdef PROFILE
         double t2 = Profile::now();
         setState(t1, "ST_TASK", Profile::get_thread_name(thread), name);
         setState(t2, "ST_TASK", Profile::get_thread_name(thread), "0");
#endif
      }

   private:
      char const* name;
      int thread;
      double t1;
   };

   static
   void init(void) {
#ifdef PROFILE
      // Initialise profiling
      setTraceType(PAJE);
      initTrace("ssids", 0, GTG_FLAG_NONE);
      // Define containers (i.e. nodes/threads)
      addContType("CT_NODE", "0", "Node");
      addContType("CT_THREAD", "CT_NODE", "Thread");
      addContainer(0.0, "C_Node0", "CT_NODE", "0", "Node 0", "0");
      for(int i=0; i<omp_get_max_threads(); ++i)
         addContainer(0.0, get_thread_name(i), "CT_THREAD", "C_Node0", get_thread_name(i), "0");
      // Define states (i.e. task types)
      addStateType("ST_TASK", "CT_THREAD", "Task");
      addEntityValue("TA_SUBTREE", "ST_TASK", "Subtree", GTG_RED);
      addEntityValue("TA_ASSEMBLE", "ST_TASK", "Assemble", GTG_GREEN);
      addEntityValue("TA_FACTOR", "ST_TASK", "Factor", GTG_BLUE);
      addEntityValue("TA_UPDATE", "ST_TASK", "Update", GTG_ORANGE);
      addEntityValue("TA_CHOL_DIAG", "ST_TASK", "CholDiag", GTG_PURPLE);
      addEntityValue("TA_CHOL_TRSM", "ST_TASK", "CholTrsm", GTG_PINK);
      addEntityValue("TA_CHOL_UPD", "ST_TASK", "CholUpd", GTG_SEABLUE);
      addEntityValue("TA_LDLT_DIAG", "ST_TASK", "LDLTDiag", GTG_PURPLE);
      addEntityValue("TA_LDLT_APPLY", "ST_TASK", "LDLTTrsm", GTG_PINK);
      addEntityValue("TA_LDLT_ADJUST", "ST_TASK", "LDLTTrsm", GTG_GRENAT);
      addEntityValue("TA_LDLT_UPD", "ST_TASK", "LDLTUpd", GTG_SEABLUE);
      clock_gettime(CLOCK_REALTIME, &tstart);
#endif
   }

   static
   void end(void) {
#ifdef PROFILE
      endTrace();
#endif
   }

   static
   double now() {
#ifdef PROFILE
      struct timespec t;
      clock_gettime(CLOCK_REALTIME, &t);
      return tdiff(tstart, t);
#else
      return 0;
#endif
   }

private:
#ifdef PROFILE
   static struct timespec tstart;
   static
   char const* get_thread_name(int thread) {
      char const* thread_name[] = {
         "Thread0", "Thread1", "Thread2", "Thread3",
         "Thread4", "Thread5", "Thread6", "Thread7",
         "Thread8", "Thread9", "Thread10","Thread11",
         "Thread12", "Thread13", "Thread14", "Thread15",
         "Thread16", "Thread17", "Thread18", "Thread19",
         "Thread20", "Thread21", "Thread22", "Thread23",
         "Thread24", "Thread25", "Thread26", "Thread27",
         "Thread28", "Thread29", "Thread30", "Thread31"
      };
      return thread_name[thread];
   }
   static
   double tdiff(struct timespec t1, struct timespec t2) {
      return (t2.tv_sec - t1.tv_sec) + 1e-9*(t2.tv_nsec - t1.tv_nsec);
   }
#endif
};


}}} /* namespace spral::ssids::cpu */
