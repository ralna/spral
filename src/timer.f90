! Module spral_timespec is seperate so it can be USEd in interface block
module spral_timespec
   use, intrinsic :: iso_c_binding

   private
   public :: timespec

   type, bind(C) :: timespec
      integer(C_INT) :: tv_sec
      integer(C_INT64_T) :: tv_nsec
   end type timespec
end module spral_timespec

module spral_timer
!$ use omp_lib
   use, intrinsic :: iso_c_binding
   use spral_timespec
   implicit none

   private
   public :: log_type, timespec ! Datatypes
   public :: clock_gettime,   & ! Interface to C routine
             log_start,       & ! Start logging a region
             log_stop,        & ! Stop logging a region
             tdiff              ! Calculate time difference between two timespec

   type log_type
      character(len=2) :: id
      integer :: thread
      integer :: narg
      integer :: arg(3)
      type(timespec) :: start
   end type log_type

   interface
      integer(C_INT) function clock_gettime(clk_id, tp) bind(C)
         use, intrinsic :: iso_c_binding
         use spral_timespec
         integer(C_INT), value, intent(in) :: clk_id
         type(timespec), intent(out) :: tp
      end function clock_gettime
   end interface

contains

real function tdiff(tp1, tp2)
   type(timespec) :: tp1, tp2

   tdiff = tp2%tv_sec - tp1%tv_sec
   tdiff = tdiff + 1e-9 * real(tp2%tv_nsec-tp1%tv_nsec)
end function tdiff

subroutine log_start(task, id, opt1, opt2, opt3)
   type(log_type), intent(out) :: task
   character(len=2), intent(in) :: id
   integer, optional, intent(in) :: opt1
   integer, optional, intent(in) :: opt2
   integer, optional, intent(in) :: opt3

   integer :: dummy

   task%thread = 0
!$ task%thread = omp_get_thread_num()

   task%id = id
   task%narg = 0
   if(present(opt1)) then
      task%narg = task%narg + 1
      task%arg(task%narg) = opt1
   endif
   if(present(opt2)) then
      task%narg = task%narg + 1
      task%arg(task%narg) = opt2
   endif
   if(present(opt3)) then
      task%narg = task%narg + 1
      task%arg(task%narg) = opt3
   endif

   dummy = clock_gettime(0, task%start)
end subroutine log_start

subroutine log_stop(task, unit_log)
   type(log_type), intent(in) :: task
   integer, intent(in) :: unit_log

   integer :: dummy
   type(timespec) :: finish

   dummy = clock_gettime(0, finish)

   select case(task%narg)
   case(0)
      write(unit_log, "(i4, 2(2i12), 1x, a2)") &
         task%thread, task%start%tv_sec, task%start%tv_nsec, finish%tv_sec, &
         finish%tv_nsec, task%id
   case(1)
      write(unit_log, "(i4, 2(2i12), 1x, a2, i12)") &
         task%thread, task%start%tv_sec, task%start%tv_nsec, finish%tv_sec, &
         finish%tv_nsec, task%id, task%arg(1:1)
   case(2)
      write(unit_log, "(i4, 2(2i12), 1x, a2, 2i12)") &
         task%thread, task%start%tv_sec, task%start%tv_nsec, finish%tv_sec, &
         finish%tv_nsec, task%id, task%arg(1:2)
   case(3)
      write(unit_log, "(i4, 2(2i12), 1x, a2, 3i12)") &
         task%thread, task%start%tv_sec, task%start%tv_nsec, finish%tv_sec, &
         finish%tv_nsec, task%id, task%arg(1:3)
   end select
end subroutine log_stop
end module spral_timer
