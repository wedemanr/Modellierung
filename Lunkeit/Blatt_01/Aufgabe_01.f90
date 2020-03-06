program Evolutionsgleichung
!Berechnung der analytischen Loesung der Evolutionsgleichung

implicit none
integer        :: i_F = 1
integer        :: i_T0 = 0
integer        :: i_dt = 1
integer        :: j1, j2
real					 :: z_alpha = 0.05
real           :: z_T

call fileopen(10,"Evolutionsgleichung.txt", "replace", "formatted")


!Outputdatei
!write(10,*) "time ; Temp"

!Berechnung

do j1 = 0,200*i_dt,i_dt
	z_T = (i_F/z_alpha)+(-i_F/z_alpha)*exp(-z_alpha*j1)
	write(10,*) j1,";", z_T
end do
close(10)
call system("gnuplot -p data_plot.plt")


end program Evolutionsgleichung




!######################################################################

subroutine fileopen(k_UN,y_filename,y_status,y_form)


implicit none

!Parameter Deklaration 
integer   , intent(in)  :: k_UN
character (len=*) , intent(in)  :: y_filename
character (len=*) , intent(in)  :: y_status
character (len=*) , intent(in) :: y_form

!Lokale Variable
integer                 :: i_error

open(unit=k_UN, file=y_filename ,iostat=i_error, status=y_status, form=y_form)
if (i_error /= 0) then
   write(*,*)"*** FEHLER ***"
   write(*,*)"beim oeffnen der Datei ", y_filename," . Errorcode:", i_error 
   stop
endif

return
end
