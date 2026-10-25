module timer_module
    use kinds_module
    use globals_module
    use settings_module
    use teletype_module
    implicit none

    real(4) :: delta_t
    real(4), dimension(2) :: times
    ! wall clock, started by start_klocka()
    integer(8) :: wall_start=-1, wall_rate=1

    public start_klocka, stockas_klocka, swatch, wall_seconds, time_string
    private
contains
    ! Start the wall clock. Call once, at program start.
    subroutine start_klocka()
        call system_clock(wall_start, wall_rate)
    end subroutine

    ! Wall-clock seconds since start_klocka()
    function wall_seconds() result(t)
        real(DP) :: t
        integer(8) :: now

        if (wall_start < 0) then
            call start_klocka()
            t=0.d0
            return
        end if
        call system_clock(now)
        t=real(now-wall_start, DP)/real(wall_rate, DP)
    end function

    ! '<label><seconds>sec ( <hours> h )', plus ' ( <days> days )' beyond two days
    function time_string(label, seconds) result(s)
        character(*), intent(in) :: label
        real(DP), intent(in) :: seconds
        character(BUFLEN) :: s

        character(16) :: days_str
        real(DP) :: hours

        hours=seconds/3600.d0
        write(s, '(a,f11.2,a,a,f6.1,a)') label, seconds, 'sec', ' (', hours, ' h )'
        if (hours > 48.d0) then
            write(days_str, '(a,f6.1,a)') ' (', hours/24.d0, ' days )'
            s=trim(s)//days_str
        end if
    end function

    subroutine stockas_klocka()
        real(4) :: etime
        character(BUFLEN) :: fdate

        ! etime: CPU time of the process, user+sys, summed over all threads
        delta_t=etime(times)
        call msg_out(repeat('-', 70))
        call msg_out(time_string('   wall time:', wall_seconds()))
        call msg_out(time_string('    cpu time:', real(delta_t, DP)))
        call msg_out(time_string('        user:', real(times(1), DP)))
        call msg_out(time_string('         sys:', real(times(2), DP)))
        call msg_out(repeat('-', 70))
        call msg_out(fdate())
    end subroutine

    ! CPU time since the previous call (dtime)
    subroutine swatch(foo)
        integer(I4), intent(in) :: foo

        real(4) :: dtime

        delta_t=dtime(times)
        if (foo > 0) then
            call msg_out(repeat('-', 70))
            call msg_out(time_string('    cpu time:', real(delta_t, DP)))
            call msg_out(time_string('        user:', real(times(1), DP)))
            call msg_out(time_string('         sys:', real(times(2), DP)))
            call msg_out(repeat('-', 70))
        end if
    end subroutine
end module


! vim:et:sw=4:ts=4
