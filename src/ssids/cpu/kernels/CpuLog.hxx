#pragma once

#include <algorithm>
#include <ctime>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spral { namespace ssids { namespace cpu {

class CpuLog {
public:
   class LogTask {
      public:
      unsigned int proc;
      struct timespec tstart;
      struct timespec tend;
      int id[4];
      uint64_t function;
      uint64_t instance;

      inline
      void end(void) {
         clock_gettime(CLOCK_REALTIME, &tend);
      }
   };

   /** Constructor */
   CpuLog(int maxtasks)
   : maxtasks(maxtasks), tasks(maxtasks), nexttask(0)
   {}

   /** Start a task executing */
   inline
   LogTask &tstart(int id0, int id1=0, int id2=0, int id3=0, uint64_t function=0, uint64_t instance=0) {
      unsigned int taskidx;
#pragma omp atomic capture
      taskidx = nexttask++;
      if(taskidx >= maxtasks) taskidx = maxtasks-1; // Out of space
      LogTask &task = tasks[taskidx];
#ifdef _OPENMP
      task.proc = omp_get_thread_num();
#else
      task.proc = 0;
#endif
      task.id[0]=id0; task.id[1]=id1; task.id[2]=id2; task.id[3]=id3;
      task.function = function; task.instance=instance;
      clock_gettime(CLOCK_REALTIME, &task.tstart);
      return task;
   }

   /** Write out an xfig profile */
   void writeFig(std::ostream &out, std::string (*taskname)(int)=getDefaultTaskName, int (*taskcolor)(int, int, int, int)=getDefaultTaskColor);

   /** Write out an aftermath profile */
   void writeAftermath(FILE *fp);

private:
   unsigned int maxtasks; //< maximum number of tasks
   std::vector<LogTask> tasks; //< vector of tasks
   unsigned int nexttask; //< atomically increment task posn counter

   static
   std::string getDefaultTaskName(int id) {
      std::ostringstream ss;
      ss << "Task ID " << id;
      return ss.str();
   }

   static
   int getDefaultTaskColor(int id0, int id1, int id2, int id3) {
      return (id0+1) % 32;
   }

};

}}} /* namespaces spral::ssids::cpu */
