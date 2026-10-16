      subroutine string_left(c,time,dt)
*
*
*     remove left blank in 20 length string
*      
      character*20 c,form
      integer i,ip
      real*8 dt,dtmp,time

*     Counts for decimal number precision: the smallest ip (at most 9)
*     for which dt*10**ip is a whole number to the rounding of dt.
      ip = 0
      dtmp = dt
 2    if(ip.lt.9.and.
     &   abs(dtmp-anint(dtmp)).gt.1.D-12*max(1.D0,abs(dtmp))) then
         ip = ip + 1
         dtmp = dtmp*10.D0
         go to 2
      end if

      if(ip.ge.10) then
         write(form,'(A5,I2,A1)') '(F20.',ip,')'
      else
         write(form,'(A5,I1,A1)') '(F20.',ip,')'
      end if

      write(c,form) time
      
*     Counts for shift
      i = 1
 1    if(c(i:i).eq.' ') then
         i = i + 1
         go to 1
      end if
      if(i.ge.2) then
         c(1:21-i) = c(i:20)
         if(ip.eq.0) then
            c(21-i:20)=' '
         else
            c(22-i:20)=' '
         end if
      end if

      return

      end
