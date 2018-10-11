module hours_above

  implicit none

contains

function hrsabove(tmin, tmax)
	!calculates cumulative hours above 17 degrees for the day
    real                    :: hrsabove !cumulative hours above 17degC
    real, intent(in)        :: tmin, tmax !daily tmin and tmax, degC
    real, dimension(24)     :: td !hourly temperature, degC
    integer                 :: i, couni
    integer, dimension(24)  :: h !hour

    data h/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,   &
			17, 18, 19, 20, 21, 22, 23/

    couni = 0

	!first generate hourly temperature from input tmin and tmax
	!Sinusoidal temperature equation from 
	!Reicosky et al. 1989 Agric. For. Meteorology 46:193-209
    do i = 1, 24

      if (((h(i) .ge. 0) .and. (h(i) .lt. 5)) .or.                     &
			((h(i) .gt. 14) .and. (h(i) .le. 24))) then

        if (h(i) .lt. 5) then
          td(i) = (tmin + tmax)/2 + ((tmax - tmin)/2)*                 &
				(cos((pi*(float(h(i))+10.0))/15.0))
        else if (h(i) .gt. 14) then
          td(i) = (tmin + tmax)/2 - ((tmax - tmin)/2)*                 &
				(cos(pi*(float(h(i)) + 14.0)/15.0))
        endif

      else if ((h(i) .ge. 5) .and. (h(i) .le. 14)) then
        
        td(i) = (tmin + tmax)/2 - ((tmax - tmin)/2)*                   &
				(cos(pi*(float(h(i)-5)/float(9))))
      
      end if

    end do

	!now calculate cumulative hours above 17
	do i = 1, 24
     
     if (td(i) .gt. 17.0) then
		couni = couni + 1
     endif
   
   end do

   hrsabove = couni

end function hrsabove
      
      
end module hours_above
