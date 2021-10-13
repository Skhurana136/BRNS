c      
c     SUBROUTINE rates
c      
      subroutine rates(j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            r(1,j) = loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))
            r(2,j) = loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))
            r(3,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxm
     +in-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)
            r(4,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxm
     +in-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)
            r(5,j) = loc_z*(kl*sp(1,j)*vel**0.58E0+kdet*sp(1,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1))
            r(6,j) = loc_z*(kl*sp(6,j)*vel**0.58E0+kdet*sp(6,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1))
            r(7,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2,
     +j)
            r(8,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,
     +j)
            r(9,j) = loc_z*kdeac*sp(1,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))
            r(10,j) = loc_z*kdeac*sp(2,j)*(1-1/(exp((doxmin-kmax1*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))
            r(11,j) = loc_z*kreac*sp(6,j)/(exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
            r(12,j) = loc_z*kreac*sp(7,j)/(exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
            r(13,j) = loc_z*km*sp(1,j)
            r(14,j) = loc_z*km*sp(2,j)
            r(15,j) = 0.1E-1*loc_z*km*sp(6,j)
            r(16,j) = 0.1E-1*loc_z*km*sp(7,j)
            r(17,j) = loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*k
     +max2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j))
            r(18,j) = loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*k
     +max2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j))
            r(19,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2
     +*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
            r(20,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2
     +*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
            r(21,j) = loc_z*(kl*sp(8,j)*vel**0.58E0+kdet*sp(8,j)/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1))
            r(22,j) = loc_z*(kl*sp(11,j)*vel**0.58E0+kdet*sp(11,j)/(exp(
     +(Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)
     +-sp(22,j))/st/Bfmax)+1))
            r(23,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(9
     +,j)
            r(24,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(1
     +2,j)
            r(25,j) = loc_z*kdeac*sp(8,j)*(1-1/(exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+1))
            r(26,j) = loc_z*kdeac*sp(9,j)*(1-1/(exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+1))
            r(27,j) = loc_z*kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/
     +st/no3min)+1)
            r(28,j) = loc_z*kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/
     +st/no3min)+1)
            r(29,j) = loc_z*km*sp(8,j)
            r(30,j) = loc_z*km*sp(9,j)
            r(31,j) = 0.1E-1*loc_z*km*sp(11,j)
            r(32,j) = 0.1E-1*loc_z*km*sp(12,j)
            r(33,j) = loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))
            r(34,j) = loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))
            r(35,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3
     +*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)
            r(36,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3
     +*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)
            r(37,j) = loc_z*(kl*sp(13,j)*vel**0.58E0+kdet*sp(13,j)/(exp(
     +(Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)
     +-sp(22,j))/st/Bfmax)+1))
            r(38,j) = loc_z*(kl*sp(16,j)*vel**0.58E0+kdet*sp(16,j)/(exp(
     +(Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)
     +-sp(22,j))/st/Bfmax)+1))
            r(39,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(1
     +4,j)
            r(40,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(1
     +7,j)
            r(41,j) = loc_z*kdeac*sp(13,j)*(1-1/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))
            r(42,j) = loc_z*kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))
            r(43,j) = loc_z*kreac*sp(16,j)/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)
            r(44,j) = loc_z*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)
            r(45,j) = loc_z*km*sp(13,j)
            r(46,j) = loc_z*km*sp(14,j)
            r(47,j) = 0.1E-1*loc_z*km*sp(16,j)
            r(48,j) = 0.1E-1*loc_z*km*sp(17,j)
            r(49,j) = loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,
     +j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/ammin)+1)
            r(50,j) = loc_z*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,
     +j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/ammin)+1)
            r(51,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20
     +,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     ++1)
            r(52,j) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21
     +,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     ++1)
            r(53,j) = loc_z*(kl*sp(20,j)*vel**0.58E0+kdet*sp(20,j)/(exp(
     +(Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)
     +-sp(22,j))/st/Bfmax)+1))
            r(54,j) = loc_z*(kl*sp(22,j)*vel**0.58E0+kdet*sp(22,j)/(exp(
     +(Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)
     +-sp(22,j))/st/Bfmax)+1))
            r(55,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2
     +1,j)
            r(56,j) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2
     +3,j)
            r(57,j) = loc_z*kdeac*sp(20,j)*(1-1/(exp((ammin-kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
            r(58,j) = loc_z*kdeac*sp(21,j)*(1-1/(exp((ammin-kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
            r(59,j) = loc_z*kreac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
            r(60,j) = loc_z*kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
            r(61,j) = loc_z*km*sp(20,j)
            r(62,j) = loc_z*km*sp(21,j)
            r(63,j) = 0.1E-1*loc_z*km*sp(22,j)
            r(64,j) = 0.1E-1*loc_z*km*sp(23,j)
            r(65,j) = loc_z*kpd*sp(18,j)
            r(66,j) = loc_z*kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1E2*(0.
     +1E0*amming-sp(5,j))/st/amming)+1)
            r(67,j) = loc_z*(250-sp(4,j))*(1-1/(exp(0.25E1*(0.4E0-ws)/st
     +)+1))
      end
