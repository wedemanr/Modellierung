program Implicit_Evo

implicit none
character(2)                       :: y_scheme
integer                            :: i_T0 = 0, i_c = 0
integer                            :: j, i_nt
real,parameter                     :: z_a = 0.05, z_F = 1., z_g=0.5
real                               :: z_Tan, z_Tnum=0,z_Tnum2=0, z_Tnum3, z_RMSE=0
real                               :: z_dt = 1, z_error=0, z_max=0, z_diff=0, z_t=0
real                               :: z_k1, z_k2, z_k3, z_k4
real                               :: Runge_Kutta

write(*,*) 'Choose your scheme -- IM:Implicit, CN: Crank-Nicolson, RK: Runge-Kutta, BF: Bashforth'
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
        open(21, file = "Analytisch.txt")
        open(22, file = "Numerisch_CrankNicolson.txt")

         do while (z_t <= 200)
                z_t = i_c * z_dt

                z_Tan = (z_F/z_a)+(-z_F/z_a)*exp(-z_a*z_t)
                write(21,*) z_t,";",z_Tan 

                z_diff = abs(z_Tan - z_Tnum)
                if (z_max <= z_diff) then
                        z_max = z_diff
                end if 
  
                z_error = z_error + (z_Tan - z_Tnum)**2
                i_c = i_c + 1
                
                
                write(22,*) z_t,";",z_Tnum
                z_Tnum = (z_Tnum*(1-z_a*z_g*z_dt)+z_F*z_dt)/(1+z_a*z_dt*(1-z_g))

        end do

        close(21)
        close(22)

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
                
        end do

else if (y_scheme == 'BF') then
        open(41, file = "Analytisch.txt")
        open(42, file = "Numerisch_Bashforth.txt")
        

        do while (z_t <= 200)
                z_t = i_c * z_dt

                z_Tan = (z_F/z_a)+(-z_F/z_a)*exp(-z_a*z_t)
                write(41,*) z_t,";",z_Tan 

                z_diff = abs(z_Tan - z_Tnum)
                if (z_max <= z_diff) then
                        z_max = z_diff
                end if 
  
                z_error = z_error + (z_Tan - z_Tnum)**2
                i_c = i_c + 1

                write(42,*) z_t,";",z_Tnum,";",z_Tnum2
                if (i_c == 1) then      !Berechne ersten Zeitschritt mit EF, da BF-Schema 2 Zeitschritte benötigt
                        z_Tnum2 = z_Tnum
                        z_Tnum = z_Tnum+z_dt*(-z_a*z_Tnum+z_F)
                        cycle
                endif
                z_Tnum3 = z_Tnum
                z_Tnum = z_Tnum*(1-(3*z_a*z_dt)/2)+z_dt*(z_F+(z_a*z_Tnum2)/2)
                z_Tnum2 = z_Tnum3
        end do
        
        close(41)
        close(42)


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
