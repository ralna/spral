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
         ::setState(t1, "ST_TASK", Profile::get_thread_name(thread), name);
         ::setState(t2, "ST_TASK", Profile::get_thread_name(thread), "0");
#endif
      }

   private:
      char const* name;
      int thread;
      double t1;
   };

   static
   void setState(char const* name, int thread) {
#ifdef PROFILE
      double t = Profile::now();
      ::setState(t, "ST_TASK", Profile::get_thread_name(thread), name);
#endif
   }

   static
   void setNullState(int thread) {
#ifdef PROFILE
      setState("0", thread);
#endif
   }

   static
   void addEvent(char const* type, int thread, char const*val) {
#ifdef PROFILE
      addEvent(now(), type, get_thread_name(thread), val);
#endif
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
      // GTG_WHITE, GTG_BLACK, GTG_DARKGREY,
      // GTG_LIGHTBROWN, GTG_LIGHTGREY, GTG_DARKBLUE, GTG_DARKPINK
      // GTG_LIGHTPINK
      addStateType("ST_TASK", "CT_THREAD", "Task");
      addEntityValue("TA_SUBTREE", "ST_TASK", "Subtree", GTG_RED);
      addEntityValue("TA_ASSEMBLE", "ST_TASK", "Assemble", GTG_GREEN);
      addEntityValue("TA_CHOL_DIAG", "ST_TASK", "CholDiag", GTG_PURPLE);
      addEntityValue("TA_CHOL_TRSM", "ST_TASK", "CholTrsm", GTG_PINK);
      addEntityValue("TA_CHOL_UPD", "ST_TASK", "CholUpd", GTG_SEABLUE);
      addEntityValue("TA_LDLT_DIAG", "ST_TASK", "LDLTDiag", GTG_PURPLE);
      addEntityValue("TA_LDLT_APPLY", "ST_TASK", "LDLTTrsm", GTG_PINK);
      addEntityValue("TA_LDLT_ADJUST", "ST_TASK", "LDLTTrsm", GTG_GRENAT);
      addEntityValue("TA_LDLT_UPDA", "ST_TASK", "LDLT Upd A", GTG_SEABLUE);
      addEntityValue("TA_LDLT_UPDC", "ST_TASK", "LDLT Upd C", GTG_ORANGE);
      addEntityValue("TA_LDLT_POST", "ST_TASK", "LDLT TPP", GTG_YELLOW);
      addEntityValue("TA_LDLT_TPP", "ST_TASK", "LDLT TPP", GTG_BLUE);
      addEntityValue("TA_ASM_PRE", "ST_TASK", "Assembly Pre", GTG_TEAL);
      addEntityValue("TA_ASM_POST", "ST_TASK", "Assembly Post", GTG_MAUVE);
      addEntityValue("TA_MISC1", "ST_TASK", "Misc 1", GTG_KAKI);
      addEntityValue("TA_MISC2", "ST_TASK", "Misc 2", GTG_REDBLOOD);
      // Define events
      addEventType("EV_AGG_FAIL", "CT_THREAD", "Aggressive pivot fail");
      // Initialise start time
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
