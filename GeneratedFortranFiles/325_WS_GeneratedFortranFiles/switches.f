c      
c     SUBROUTINE switches
c      
      subroutine switches(j)
          include 'common_geo.inc'
          include 'common.inc'
          integer j
          real*8 x_pos
          real*8 y_pos
          real*8 z_pos
          x_pos = x(j)
            if (0.lt.0.61D0-z_pos) then
              loc_z = 0.1D1
              else
              loc_z = 0.D0
            endif
      end
