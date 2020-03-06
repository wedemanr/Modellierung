program Euler
!Berechnung der analytischen Loesung der Evolutionsgleichung

implicit none
integer        			   :: i_T0 = 0, i_c = 0
integer        			   :: j, i_nt
real					 			   :: z_a = 0.05, z_F = 1.
real           			   :: z_Tan, z_Tnum=0, z_RMSE=0
real                   :: z_dt = 1, z_error=0, z_max=0, z_diff=0, z_t=0

call fileopen(10,"Evolutionsgleichung.txt", "replace", "formatted")
open(20, file="Output.txt")


!Berechnung
!Explicit (Euler)	
	
do while (z_t <= 200)
	z_t = i_c * z_dt
	
	z_Tan = (z_F/z_a)+(-z_F/z_a)*exp(-z_a*z_t)
	write(10,*) z_t,";",z_Tan 
		
	z_diff = abs(z_Tan - z_Tnum)
  if (z_max <= z_diff) then
		z_max = z_diff	
	end if 
	
	z_error = z_error + (z_Tan - z_Tnum)**2
	i_c = i_c + 1

	write(20,*) z_t,";",z_Tnum
	z_Tnum = z_Tnum + z_dt * (-z_a*z_Tnum + z_F)

end do

z_RMSE = sqrt(z_error/i_c)
write(*,*) "RMSE:",z_RMSE
write(*,*) "Maxinumsnorm:", z_max


end program Euler
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
