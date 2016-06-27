#include "CpuLog.hxx"

#ifdef BUB_AFTERMATH
extern "C" {
#include <trace_file.h>
#include <ansi_extras.h>
#include <convert.h>
}
#endif

namespace spral { namespace ssids { namespace cpu {

/* Anonymous namespace to hold xfig profile writing routines */
namespace {

   /** Returns true if t1 < t2 */
   bool timeSpecLT(const struct timespec &t1, const struct timespec &t2) {
      if(t1.tv_sec < t2.tv_sec) return true;
      if(t1.tv_sec > t2.tv_sec) return false;
      // Otherwise, t1.tv_sec == t2.tv_sec
      return (t1.tv_nsec < t2.tv_nsec);
   }

   /** Return number of nanoseconds between t2 and t1 */
   long tdiff(const struct timespec &t1, const struct timespec &t2) {
      return 1000000000L*(t2.tv_sec-t1.tv_sec) + t2.tv_nsec-t1.tv_nsec;
   }

   class FigWriter {
   private:
      std::ostream* output;
   public:
      FigWriter(std::ostream &output) :
            output(&output) {
         output << "#FIG 3.2"    << std::endl // Filetype version
                << "Landscape"   << std::endl // Orientation
                << "Center"      << std::endl // Centered
                << "Metric"      << std::endl // Use cm
                << "A4"          << std::endl // Paper type
                << "100.00"      << std::endl // Scaling factor
                << "Single"      << std::endl // One page only
                << "-2"          << std::endl // Mp transparent colour
                << "1200 2"      << std::endl;// Std resolution; upr left origin
      }

      void box(int x1, int y1, int x2, int y2, int pen_colour, int fill_colour, int depth=50) {
         //printf("Box (%d, %d), (%d, %d) pen %d fill %d depth %d\n", x1, y1, x2, y2, pen_colour, fill_colour, depth);
         if(x1==x2 || y1==y2) return; // Skip 0-size box
         *output << "2 2 0 1 " << pen_colour << " " << fill_colour << " " << depth << "-1 20 0.000 0 0 -1 0 0 5" << std::endl;
         *output << "   " <<
            x1 << " " << y1 << " " <<
            x2 << " " << y1 << " " <<
            x2 << " " << y2 << " " <<
            x1 << " " << y2 << " " <<
            x1 << " " << y1 << std::endl;
      }

      void text(int x, int y, std::string text) {
         *output << "4 0 0 50 -1 0 12 0.0000 4 150 345 " <<
            x << " " << y << " " << text << "\\001" << std::endl;
      }
   };

   bool TaskStartComparator (const typename CpuLog::LogTask &i, const typename CpuLog::LogTask &j) {
      if(i.tstart.tv_sec < j.tstart.tv_sec) return true;
      if(i.tstart.tv_sec > j.tstart.tv_sec) return false;
      // Otherwise, i.tv_sec == j.tv_sec
      return (i.tstart.tv_nsec < j.tstart.tv_nsec);
   }

   struct ProcProfile {
      typedef typename CpuLog::LogTask LogTaskSpec;
      typedef std::vector<LogTaskSpec> TaskVector;
      TaskVector tasks;

      struct ProfileRow {
         std::vector<LogTaskSpec> tasks;
         bool addTask(LogTaskSpec t) {
            // Check we don't overlap with last task
            if(tasks.size() != 0) {
               struct timespec &freefrom = (tasks.end()-1)->tend;
               if( timeSpecLT(t.tstart,freefrom) ) return false;
            }
            // If not, add task to row
            tasks.push_back(t);
            return true;
         }
      };
      typedef std::vector<struct ProfileRow> RowVector;
      RowVector rows;

      void addTask(LogTaskSpec t) { tasks.push_back(t); }

      void process() {
         // Sort tasks by start time
         std::sort(tasks.begin(), tasks.end(), TaskStartComparator);
         // Create rows
         for(typename TaskVector::iterator t=tasks.begin(); t!=tasks.end(); ++t) {
            for(typename RowVector::iterator row=rows.begin(); row!=rows.end(); ++row) {
               // Try adding to row, if succeed go to next loop i
               // FIXME  goto done_outer
               if (row->addTask(*t)) goto done_outer;
            }
            // Failed to add to any row, create a new one
            rows.push_back(ProfileRow());
            (rows.end()-1)->addTask(*t);
            done_outer: continue; // Noop
         }
      }
   };
};

/** Write out an xfig profile */
void CpuLog::writeFig(std::ostream &out, std::string (*taskname)(int), int (*taskcolor)(int, int, int, int)) {
   /* Create per-CTA profiles */
   typedef std::vector< ProcProfile > ProfVector;
#ifdef _OPENMP
   const int nproc = omp_get_max_threads();
#else
   const int nproc = 1;
#endif
   ProfVector profs(nproc);
   struct timespec clkzero;
   clock_gettime(CLOCK_REALTIME, &clkzero);
   struct timespec clkmax;
   clkmax.tv_sec = 0;
   clkmax.tv_nsec = 0;
   for(unsigned int tid=0; tid<nexttask; ++tid) {
      LogTask &t = tasks[tid];
      if( timeSpecLT(clkmax, t.tend) ) clkmax = t.tend;
      if( timeSpecLT(t.tstart, clkzero) ) clkzero = t.tstart;
      //printf("Got task PROC %u (%u) %d.%ld:%d.%ld\n", t.proc, t.id[0], t.tstart.tv_sec, t.tstart.tv_nsec, t.tend.tv_sec, t.tend.tv_nsec);
      profs[t.proc].addTask(t);
   }
   /* Output Fig file */
   const int left=400, right=13200; // Edges of visible xfig area
   const int top=200; // Edges of visible xfig area
   FigWriter fig(out);
   int rtop = top; // Current top of row
   int rheight = 200, rgap=50;
   float xscale = (0.0+right-left) / tdiff(clkzero, clkmax);
   int sm=0;
   int maxid=0;
   bool cused[32];
   for(int i=0; i<32; i++) cused[i] = false;
   for(typename ProfVector::iterator prof=profs.begin(); prof!=profs.end(); ++prof) {
      if(prof->tasks.size() == 0) continue;
      std::ostringstream name;
      name << "SM" << sm++;
      fig.text(0, rtop+rheight, name.str());
      prof->process();
      for(auto row=prof->rows.begin(); row!=prof->rows.end(); ++row) {
         for(auto t=row->tasks.begin(); t!=row->tasks.end(); ++t) {
            fig.box(left+tdiff(clkzero,t->tstart)*xscale, rtop,
                    left+tdiff(clkzero, t->tend)*xscale, rtop+rheight,
                    0, taskcolor(t->id[0], t->id[1], t->id[2], t->id[3]));
            if(t->id[0] > maxid) maxid = t->id[0];
            cused[t->id[0]] = true;
         }
         rtop += rheight + rgap;
      }
   }
   /* Key */
   rtop += 2*(rheight + rgap);
   int scale_top = rtop;
   for(int i=0; i<=maxid; i++) {
      if(i<32 && !cused[i]) continue;
      fig.box(left, rtop, left+rheight, rtop+rheight, 0, taskcolor(i, 0, 0, 0));
      fig.text(left+rheight+rgap, rtop+rheight, taskname(i));
      rtop += rheight + rgap;
   }
   /* Scale */
   int clk;
   for(clk=1000; clk*xscale<(right-left)/10; clk*=2); // ensure length is good
   int scale_right = right - 1000;
   fig.box(scale_right-clk*xscale, scale_top, scale_right, scale_top+rgap, 0, 0);
   fig.box(scale_right-50, scale_top-rgap, scale_right, scale_top+2*rgap, 0, 0);
   fig.box(scale_right-clk*xscale, scale_top-rgap, scale_right-clk*xscale+50, scale_top+2*rgap, 0, 0);
   std::ostringstream ss;
   ss << clk;
   fig.text(scale_right-clk*xscale/2-500, scale_top+2*rheight, ss.str());
}

/* Anonymous namespace for aftermath support functions */
#ifdef BUB_AFTERMATH
namespace {

class TraceHeader {
public:
   TraceHeader(struct tm *now) {
      header.magic = TRACE_MAGIC;
      header.version = TRACE_VERSION;
      header.day = now->tm_mday;
      header.month = now->tm_mon+1;
      header.year = now->tm_year+1900;
      header.hour = now->tm_hour;
      header.minute = now->tm_min;
   }

   void write(FILE *fp) {
      write_struct_convert(fp, &header, sizeof(header), trace_header_conversion_table, 0);
   }
private:
   struct trace_header header;
};

class TraceSingleEvent {
public:
   TraceSingleEvent(
         /// SINGLE_TYPE_TEXEC_START or SINGLE_TYPE_TEXEC_END
         enum single_event_type type,
         /// Time of event
         uint64_t time,
         /// CPU of event
         uint64_t cpu,
         /// Worker of event
         uint64_t worker=0,
         /// Type of task, e.g. function address
         uint64_t active_task=0,
         /// Instance of task, e.g. pointer to data structure
         uint64_t what=0)
   {
      event.header.type = EVENT_TYPE_SINGLE;
      event.header.time = time;
      event.header.cpu = cpu;
      event.header.worker = worker;
      event.header.active_task = active_task;
      event.header.active_frame = what;
      event.what = what;
      event.type = type;
   }

   void write(FILE *fp) {
      write_struct_convert(fp, &event, sizeof(event), trace_single_event_conversion_table, 0);
   }
private:
   struct trace_single_event event;
};

class TraceStateEvent {
public:
   TraceStateEvent(
         /// WORKER_STATE_SEEKING or WORKER_STATE_TASKEXEC
         enum worker_state state,
         /// Start time of state
         uint64_t start_time,
         /// End time of state
         uint64_t end_time,
         /// CPU of event
         uint64_t cpu,
         /// Worker of event
         uint64_t worker=0,
         /// Type of task, e.g. function address
         uint64_t active_task=0,
         /// Instance of task, e.g. pointer to data structure
         uint64_t what=0)
   {
      event.header.type = EVENT_TYPE_STATE;
      event.header.time = start_time;
      event.header.cpu = cpu;
      event.header.worker = worker;
      event.header.active_task = active_task;
      event.header.active_frame = what;
      event.end_time = end_time;
      event.state = state;
   }

   void write(FILE *fp) {
      write_struct_convert(fp, &event, sizeof(event), trace_state_event_conversion_table, 0);
   }
private:
   struct trace_state_event event;
};

class AfterMath {
public:
   AfterMath(FILE *fp)
   : fp(fp)
   {
      time_t now = time(NULL);
      TraceHeader(localtime(&now)).write(fp);
   }

   void writeTask(uint64_t start_time, uint64_t end_time, uint64_t cpu,
         uint64_t task_type=0, uint64_t task_instance=0) {
      TraceSingleEvent(SINGLE_TYPE_TEXEC_START, start_time, cpu, 0, task_type, task_instance).write(fp);
      TraceSingleEvent(SINGLE_TYPE_TEXEC_END, end_time, cpu, 0, task_type, task_instance).write(fp);
      TraceStateEvent(WORKER_STATE_TASKEXEC, start_time, end_time, cpu, 0, task_type, task_instance);
   }

private:
   FILE *fp;
};

} /* anon */
#endif /* ifdef BUB_AFTERMATH */

/** Write out an aftermath profile */
void CpuLog::writeAftermath(FILE *fp) {
#ifndef BUB_AFTERMATH
   throw std::runtime_error("Aftermath support not selected at compile time\n");
#else
   /* Determine clock zero and configuration space */
   struct timespec clkzero;
   clock_gettime(CLOCK_REALTIME, &clkzero);
   struct timespec clkmax;
   int idmax[3] = {0, 0, 0};
   clkmax.tv_sec = 0;
   clkmax.tv_nsec = 0;
   for(unsigned int tid=0; tid<nexttask; ++tid) {
      LogTask &t = tasks[tid];
      if( timeSpecLT(clkmax, t.tend) ) clkmax = t.tend;
      if( timeSpecLT(t.tstart, clkzero) ) clkzero = t.tstart;
      for(int i=0; i<3; i++)
         if(t.id[i+1] > idmax[i]) idmax[i] = t.id[i+1];
   }
   
   /* Initialize file */
   AfterMath trace(fp);

   /* Loop over tasks, storing them to disk */
   for(unsigned int tid=0; tid<nexttask; ++tid) {
      LogTask &t = tasks[tid];
      if( timeSpecLT(clkmax, t.tend) ) clkmax = t.tend;
      if( timeSpecLT(t.tstart, clkzero) ) clkzero = t.tstart;
      if(!t.function) t.function = t.id[0];
      if(!t.instance) t.instance = idmax[1]*t.id[3]+idmax[0]*t.id[2]+t.id[1];
      trace.writeTask( tdiff(clkzero, t.tstart), tdiff(clkzero, t.tend), t.proc,
            t.function, t.instance);
   }
#endif
}

}}} /* namespaces spral::ssids::cpu */
