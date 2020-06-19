       program Lorenz
        
       implicit none
       real, parameter                 :: Z_dt=0.01, Z_TS=100
       real                            :: z_t=0, z_x=1, z_y=1, z_z=1
       real                            :: z_i1,z_i2,z_i3,z_i4
       real                            :: z_k1,z_k2,z_k3,z_k4
       real                            :: z_j1,z_j2,z_j3,z_j4
       real                            :: RK_X,RK_Y,RK_Z
       integer                         :: i_c=0
       
       
       open(10,file="Lorenz.txt")

       do while (z_t <= Z_TS)
               write(*,*) "In progress:",(z_t/Z_TS)*100,"%"
               z_t = i_c*Z_dt

               write(10,*) z_x,";",z_y,";",z_z

               z_i1 = Z_dt * RK_X(z_x,z_y)
               z_j1 = Z_dt * RK_Y(z_x,z_y,z_z)
               z_k1 = Z_dt * RK_Z(z_x,z_y,z_z)
               
               z_i2 = Z_dt * RK_X(z_x+z_i1/2,z_y+z_j1/2)
               z_j2 = Z_dt * RK_Y(z_x+z_i1/2,z_y+z_j1/2,z_z+z_k1/2) 
               z_k2 = Z_dt * RK_Y(z_x+z_i1/2,z_y+z_j1/2,z_z+z_k1/2) 
                  
               z_i3 = Z_dt * RK_X(z_x+z_i2/2,z_y+z_j2/2)
               z_j3 = Z_dt * RK_Y(z_x+z_i2/2,z_y+z_j2/2,z_z+z_k2/2) 
               z_k3 = Z_dt * RK_Y(z_x+z_i2/2,z_y+z_j2/2,z_z+z_k2/2)

               z_i4 = Z_dt * RK_X(z_x+z_i3,z_y+z_j3)
               z_j4 = Z_dt * RK_Y(z_x+z_i3,z_y+z_j3,z_z+z_k3)
               z_k4 = Z_dt * RK_Y(z_x+z_i3,z_y+z_j3,z_z+z_k3)
 
               z_x = z_x + (z_i1 + 2*z_i2 + 2*z_i3 + z_i4)/6
               z_y = z_y + (z_j1 + 2*z_j2 + 2*z_j3 + z_j4)/6
               z_z = z_z + (z_k1 + 2*z_k2 + 2*z_k3 + z_k4)/6
               
               i_c = i_c + 1 
       end do
       close(10)

       end program Lorenz



       !Functions
       real function RK_X(p_x,p_y)
               implicit none
               real, intent(in) :: p_x,p_y
               real, parameter  :: Z_sigma=10
                
               RK_X = -Z_sigma*p_x + Z_sigma*p_y
               return
       end function RK_X

       real function RK_Y(p_x,p_y,p_z)
               implicit none
               real, intent(in) :: p_x,p_y,p_z
               real, parameter  :: Z_r=28

               RK_Y = -p_x*p_z + Z_r*p_x - p_y
               return
       end function RK_Y

       real function RK_Z(p_x,p_y,p_z)
               implicit none
               real, intent(in) :: p_x,p_y,p_z
               real, parameter  :: Z_b=8/3

               RK_Z = p_x*p_y - Z_b*p_z
               return
       end function RK_Z


