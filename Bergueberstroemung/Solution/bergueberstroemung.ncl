;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

   x = ispan(-80000,80000,1000)
   y = ispan(0,10000,100)
   ; hill dimensions
   Lzb = 300. ;m
   Lxb = 3000. ;m

   ; atmospheric conditions
   U = 10 ;m/s
   theta0 = 290 ;K
   dthetadz = .005 ;K/m

   ; misc
   g = 9.81 ;m/s^2

   ; definition bell shaped hill
   zs = Lzb * ( Lxb^2 ) / ( Lxb^2 + x^2 )

   ; Brunt-Vaisala frequency
   N = sqrt( g/theta0 * dthetadz )

   ; displacement of stream lines
   delta = zs*cos(N/U*(y-zs))-x/Lxb*zs*sin(N/U*(y-zs))

   ; temperature
   theta = theta0 * ( 1 + ( N^2 * (y-delta)/g ) )

   ; horizontal wind
   u = U * ( 1 - dddz )

   ; dddz
   dddz = -zs * sin( N/U * (y-zs) ) * N/U - x/Lxb * zs * cos( N/U * (y-zs) ) * N/U

   ; vertical wind
   w = U * dddx

   dddx = - ( 2 * Lxb^3 * Lzb^2 * N * x^2 * cos(N/U*(y-zs(x))) ) / ( U * (Lxb^2+x^2)^3 ) \
          - ( 2 * Lxb^2 * Lzb       * x   * cos(N/U*(y-zs(x))) ) /       (Lxb^2+x^2)^2   \
          - ( 2 * Lxb^4 * Lzb^2 * N * x   * sin(N/U*(y-zs(x))) ) / ( U * (Lxb^2+x^2)^3 ) \
          + ( 2 * Lxb   * Lzb       * x^2 * sin(N/U*(y-zs(x))) ) /       (Lxb^2+x^2)^2   \
          - (     Lxb   * Lzb             * sin(N/U*(y-zs(x))) ) /       (Lxb^2+x^2)

   ; plot stuff

      x@long_name = "x [km]"
      y@long_name = "z [m]"
   wks = gsn_open_wks("pdf", "bergueberstroemung_ncl")
      res = True
      ;res@xyLineColor = "indianred4"
      ;res@xyLineThicknessF = 3
      ;res@tmXBMode = "Explicit"
      ;res@tmXBValues = (/-50,-40,-30/)
      ;res@tmXBLabel = (/"A","B","C"/)
   plot = gsn_csm_xy(wks,x,delta,res)
   ;splot delta
   ;splot theta
   ;splot u
   ;splot w

end
