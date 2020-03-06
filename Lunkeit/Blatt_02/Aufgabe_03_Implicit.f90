program Implicit_Evo

implicit none
integer        			   :: i_T0 = 0, i_c = 0
integer        			   :: j, i_nt
real					 			   :: z_a = 0.05, z_F = 1.
real           			   :: z_Tan, z_Tnum=0, z_RMSE=0
real                   :: z_dt = 1, z_error=0, z_max=0, z_diff=0, z_t=0

open(10,file = "Evolutionsgleichung.txt")
open(20,file = "Output.txt")

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
	z_Tnum = (z_Tnum + (z_dt * z_F))/(1+(z_dt * z_a))

end do

z_RMSE = sqrt(z_error/i_c)
write(*,*) "RMSE:",z_RMSE
write(*,*) "Maxinumsnorm:", z_max



end program Implicit_Evo
