!u301027 Ralf Wedemann 04.18


!******************************************************************************************************!
!**********************                          Oromod                          **********************!
!******************************************************************************************************!


module oromod
implicit none

real, parameter :: ZPI = 3.141592653589793_8 ! Pi
real, parameter :: Z2PI = ZPI + ZPI
real, parameter :: ZG = 9.80665       ! Gravity
real, parameter :: ZK = 273.15        ! 0° in Kelvin

end module oromod

!******************************************************************************************************!
!**********************                         Subroutine                       **********************!
!******************************************************************************************************!

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

!******************************************************************************************************!

subroutine get_index(k_lon, k_lat,p_lon_koordinate, p_lat_koordinate,k_anzahl, k_lon_out, k_lat_out)
use oromod
implicit none

!Eingelesene Variablen
integer, intent(in)       :: k_lon, k_lat, k_anzahl
real , intent(in)					:: p_lon_koordinate(k_anzahl), p_lat_koordinate(k_anzahl)
integer, intent(out)  		:: k_lon_out(k_anzahl), k_lat_out(k_anzahl)

!Lokale Variablen
real                      :: p_lon_neu(k_anzahl), p_lat_neu(k_anzahl)
integer                   :: j1

p_lon_neu(:) = 0
p_lat_neu(:) = 0

do j1 = 1,k_anzahl
	p_lat_neu(j1) = 90. - p_lat_koordinate(j1)
	if (p_lon_koordinate(j1) < 0) then
		p_lon_neu(j1) = 360 + p_lon_koordinate(j1)
	else
		p_lon_neu(j1) = p_lon_koordinate(j1)
	endif
	k_lat_out(j1) = AINT(p_lat_neu(j1) / (180. / (k_lat-1)) +1)
  k_lon_out(j1) = AINT(p_lon_neu(j1) / (360. / k_lat) +1)
end do


return
end

!******************************************************************************************************!

subroutine get_value(klon,klat,plambda,pphi,pfield,pval)
use oromod
implicit none

!Eingelesen Variablen
integer, intent(in)     :: klon, klat
real, intent(in)        :: plambda, pphi,pfield(klon,klat)
real, intent(out)       :: pval

!Lokale Variablen
real                    :: plambda_neu, pphi_neu, plambda_neu1, pphi_neu1
real                    :: pdl, pdp, pl1, pl2
integer                 :: klambda_neu, kphi_neu
integer                 :: j1

plambda_neu = 0
pphi_neu = 0

pphi_neu = 90. - pphi
if (plambda < 0) then
	plambda_neu = 360 + plambda
else
	plambda_neu = plambda
endif

pphi_neu1 = pphi_neu / (180. / (klat-1)) +1
plambda_neu1 = plambda_neu / (360. / klon) +1

kphi_neu = AINT(pphi_neu / (180. / (klat-1)) +1)
klambda_neu = AINT(plambda_neu / (360. / klon) +1)

pdp = pphi_neu1 - kphi_neu
pdl = plambda_neu1 - klambda_neu


!write(*,*) "plambda", plambda
!write(*,*) "pphi", pphi
!write(*,*) "pphi_neu1",pphi_neu1
!write(*,*) "kphi_neu", kphi_neu
!write(*,*) "pdp", pdp
!write(*,*) "plambda_neu1",plambda_neu1
!write(*,*) "klambda_neu",klambda_neu
!write(*,*) "pdl", pdl

!Berechnung für bestimmte Gitterpunkte
pl1 = (1-pdl)*pfield(klambda_neu,kphi_neu) + pdl*pfield(klambda_neu+1,kphi_neu)
pl2 = (1-pdl)*pfield(klambda_neu,kphi_neu+1) + pdl*pfield(klambda_neu+1,kphi_neu+1)
pval = (1-pdp)*pl1 + pdp*pl2

!write(*,*)"pfield(klambda_neu,kphi_neu)",pfield(klambda_neu,kphi_neu)
!write(*,*)"pfield(klambda_neu+1,kphi_neu)",pfield(klambda_neu+1,kphi_neu)
!write(*,*)"pfield(klambda_neu,kphi_neu+1)",pfield(klambda_neu,kphi_neu+1)
!write(*,*)"pfield(klambda_neu+1,kphi_neu+1)",pfield(klambda_neu+1,kphi_neu+1)


!write(*,*) "pl1",pl1 
!write(*,*) "pl2",pl2
!write(*,*) "pval",pval
!write(*,*) "(1-pdp)",(1-pdp)
!write(*,*) 

return
end

!******************************************************************************************************!
!**********************                          Function                        **********************!
!******************************************************************************************************!


real function windstaerke(p_u,p_v)

implicit none


!Eingelesene Variablen
real, intent(in)     :: p_u       ! u-Windkomponente
real, intent(in)     :: p_v       ! v-Windkomponente

                         
windstaerke = sqrt(p_u**2 + p_v**2)

return
end function windstaerke 

!******************************************************************************************************!


real function windrichtung(p_u,p_v)

implicit none

!Eingelesene Variablen
real, intent(in)    :: p_u
real, intent(in)    :: p_v

!Lokale Variablen
real, parameter :: ZPI = 3.141592653589793_8 ! Pi


windrichtung = atan2(-p_u,-p_v)
windrichtung = windrichtung*(180/ZPI)
windrichtung = windrichtung+360
windrichtung = mod(windrichtung,real(360))


return
end function windrichtung

!******************************************************************************************************!

real function zonalwind(p_s,p_r)

implicit none

!Eingelesene Variablen
real, intent(in)   :: p_s
real, intent(in)   :: p_r

!Lokale Variablen
real, parameter :: ZPI = 3.141592653589793_8 ! Pi

zonalwind = -p_s*sin(p_r*(ZPI/180))

return
end function zonalwind

!******************************************************************************************************!

real function meridionalwind(p_s,p_r)

implicit none

!Eingelesene Variablen
real, intent(in)   :: p_s
real, intent(in)   :: p_r

!Lokale Variablen
real, parameter :: ZPI = 3.141592653589793_8 ! Pi

meridionalwind = -cos(p_r*(ZPI/180))*p_s

return
end function meridionalwind

!******************************************************************************************************!

real function land_max(p_feld, p_mask, k_dim1, k_dim2)

implicit none 

!Eingelesen Variablen
real, intent(in)      				:: p_feld(k_dim1,k_dim2)
real, intent(in)      				:: p_mask(k_dim1,k_dim2)
integer, intent(in)   				:: k_dim1, k_dim2

!Lokale Variablen
integer               				:: j1,j2
real,dimension(k_dim1,k_dim2) :: z_help




do j1=1,k_dim1
	do j2=1,k_dim2
		z_help(j1,j2)=p_feld(j1,j2)*p_mask(j1,j2)  !Temperatur über Land
	end do	
end do
land_max=maxval(z_help(:,:))

return
end function land_max

!******************************************************************************************************!


real function see_max(p_feld, p_mask, k_dim1, k_dim2)

implicit none 

!Eingelesen Variablen
real, intent(in)      					:: p_feld(k_dim1, k_dim2)
real, intent(in)     					  :: p_mask(k_dim1, k_dim2)
integer, intent(in)   					:: k_dim1, k_dim2

!Lokale Variablen
integer               					:: j1,j2
real,dimension(k_dim1,k_dim2) 	:: z_help



do j1=1,k_dim1
	do j2=1,k_dim2
		z_help(j1,j2)=p_feld(j1,j2)*(-1*(p_mask(j1,j2)-1))  !Temperatur über Wasser
	end do
end do
see_max=maxval(z_help(:,:))

return
end function see_max

