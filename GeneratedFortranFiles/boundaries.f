c      
c     SUBROUTINE boundaries
c      
      subroutine boundaries()
        include 'common_geo.inc'
        include 'common.inc'
        j = 1
          spb(1,1) = 0
          spb(2,1) = 2
          spb(3,1) = 800
          spb(4,1) = 250
          spb(5,1) = 60
          spb(6,1) = 0
          spb(7,1) = 2
          spb(8,1) = 0
          spb(9,1) = 2
          spb(10,1) = 250
          spb(11,1) = 0
          spb(12,1) = 5
          spb(13,1) = 0
          spb(14,1) = 2
          spb(15,1) = 1500
          spb(16,1) = 0
          spb(17,1) = 2
          spb(18,1) = 5
          spb(19,1) = 20
          spb(20,1) = 0
          spb(21,1) = 2
          spb(22,1) = 0
          spb(23,1) = 2
          ibc(1,1) = 2
          ibc(2,1) = 0
          ibc(3,1) = 0
          ibc(4,1) = 0
          ibc(5,1) = 0
          ibc(6,1) = 2
          ibc(7,1) = 0
          ibc(8,1) = 2
          ibc(9,1) = 0
          ibc(10,1) = 0
          ibc(11,1) = 2
          ibc(12,1) = 0
          ibc(13,1) = 2
          ibc(14,1) = 0
          ibc(15,1) = 0
          ibc(16,1) = 2
          ibc(17,1) = 0
          ibc(18,1) = 0
          ibc(19,1) = 0
          ibc(20,1) = 2
          ibc(21,1) = 0
          ibc(22,1) = 2
          ibc(23,1) = 0
        j = nx
          spb(1,2) = 0
          spb(2,2) = 0
          spb(3,2) = 0
          spb(4,2) = 0
          spb(5,2) = 0
          spb(6,2) = 0
          spb(7,2) = 0
          spb(8,2) = 0
          spb(9,2) = 0
          spb(10,2) = 0
          spb(11,2) = 0
          spb(12,2) = 0
          spb(13,2) = 0
          spb(14,2) = 0
          spb(15,2) = 0
          spb(16,2) = 0
          spb(17,2) = 0
          spb(18,2) = 0
          spb(19,2) = 0
          spb(20,2) = 0
          spb(21,2) = 0
          spb(22,2) = 0
          spb(23,2) = 0
          ibc(1,2) = 1
          ibc(2,2) = 1
          ibc(3,2) = 1
          ibc(4,2) = 1
          ibc(5,2) = 1
          ibc(6,2) = 1
          ibc(7,2) = 1
          ibc(8,2) = 1
          ibc(9,2) = 1
          ibc(10,2) = 1
          ibc(11,2) = 1
          ibc(12,2) = 1
          ibc(13,2) = 1
          ibc(14,2) = 1
          ibc(15,2) = 1
          ibc(16,2) = 1
          ibc(17,2) = 1
          ibc(18,2) = 1
          ibc(19,2) = 1
          ibc(20,2) = 1
          ibc(21,2) = 1
          ibc(22,2) = 1
          ibc(23,2) = 1
      end
