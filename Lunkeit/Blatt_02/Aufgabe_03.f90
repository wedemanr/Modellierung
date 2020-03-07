program Implicit_Evo

implicit none
character(2)                       :: y_scheme
integer                            :: i_T0 = 0, i_c = 0
integer                            :: j, i_nt
real,parameter                     :: z_a = 0.05, z_F = 1.
real                               :: z_Tan, z_Tnum=0, z_RMSE=0
real                               :: z_dt = 1, z_error=0, z_max=0, z_diff=0, z_t=0
real                               :: z_k1, z_k2, z_k3, z_k4
real                               :: Runge_Kutta

write(*,*) 'Choose your scheme -- IM:Implicit, CN: Cranḱ-Nicolsen, RK: Runge-Kutta, BF: Bashford'
read(*,*) y_scheme

if (y_scheme == 'IM') then
        open(11,file = "Analytisch.txt")
        open(12,file = "Numerisch_Implicit.txt")

        do while (z_t <= 200)
                z_t = i_c * z_dt

                z_Tan = (z_F/z_a)+(-z_F/z_a)*exp(-z_a*z_t)
                write(11,*) z_t,";",z_Tan 

                z_diff = abs(z_Tan - z_Tnum)
                if (z_max <= z_diff) then
                        z_max = z_diff
                end if 
  
                z_error = z_error + (z_Tan - z_Tnum)**2
                i_c = i_c + 1

                write(12,*) z_t,";",z_Tnum
                z_Tnum = (z_Tnum + (z_dt * z_F))/(1+(z_dt * z_a))

        end do

        close(11)
        close(12)

else if (y_scheme == 'CN') then
        write(*,*) 'No Calculation Code'

else if (y_scheme == 'RK') then
        open(31, file = "Analytisch.txt")
        open(32, file = "Numerisch_RungeKutta.txt")
        do while (z_t <= 200)
                z_t = i_c * z_dt
                
                z_Tan = (z_F/z_a)+(-z_F/z_a)*exp(-z_a*z_t)
                write(31,*) z_t,";",z_Tan 
                
                z_diff = abs(z_Tan - z_Tnum)
                if (z_max <= z_diff) then
                        z_max = z_diff
                end if

                z_error = z_error + (z_Tan -z_Tnum)**2
                i_c = i_c + 1

                
                !Runge Kutta
                write(32,*) z_t,";",z_Tnum
                z_k1 = z_dt * Runge_Kutta(z_Tnum)
                z_k2 = z_dt * Runge_Kutta(z_Tnum + z_k1/2)
                z_k3 = z_dt * Runge_Kutta(z_Tnum + z_k2/2)
                z_k4 = z_dt * Runge_Kutta(z_Tnum + z_k3)
                z_Tnum = z_Tnum + (z_k1 + 2*z_k2 + 2*z_k3 + z_k4)/6
                write(*,*) z_k1, z_k2, z_k3, z_k4, z_Tnum
        end do

else if (y_scheme == 'BF') then
        write(*,*) 'No Calculation Code'

else
       write(*,*) 'ERROR 404: No scheme found'
       stop
endif

z_RMSE = sqrt(z_error/i_c)
write(*,*) "Schema:", y_scheme
write(*,*) "RMSE:",z_RMSE
write(*,*) "Maximumsnorm:", z_max

end program Implicit_Evo




!Funktionen
real function Runge_Kutta(p_T)

implicit none

real, intent(in) :: p_T
real, parameter :: Z_Af = 0.05, Z_Ff = 1.

Runge_Kutta = -Z_Af*p_T+Z_Ff
return
end function Runge_Kutta
