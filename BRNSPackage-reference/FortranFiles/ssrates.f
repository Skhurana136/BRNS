c      
c     SUBROUTINE ssrates
c      
      subroutine ssrates(rat,drdc,isp,j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            if (isp.eq.1) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(
     +1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-
     +kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)-loc_z*(kl*sp(1,j)*vel**0.58D0+kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1))+loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2,j)-lo
     +c_z*kdeac*sp(1,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))+loc_z*kreac*sp(6,j)/(exp((
     +doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)-loc_z*km*sp(1,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_
     +z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(1,j)/(
     +exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(1
     +6,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)
     +-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-
     +loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)*sp(2,j)-loc_z*kdeac*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))-loc_z*km
            endif
            if (isp.eq.2) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(
     +2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-
     +kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)+loc_z*(kl*sp(1,j)*vel**0.58D0+kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1))-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2,j)-lo
     +c_z*kdeac*sp(2,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))+loc_z*kreac*sp(7,j)/(exp((
     +doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)-loc_z*km*sp(2,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_
     +z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*kdeac*(1-1/(exp((d
     +oxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/do
     +xmin)+1))-loc_z*km
            endif
            if (isp.eq.3) then
            rat = -fdoco*(loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+loc_z/(exp((doxmin-kmax1*sp(
     +3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*
     +sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))-fdocn*(lo
     +c_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j))+loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j)))-fdocs*(loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     ++loc_z*kpd*sp(18,j)
            drdc = -fdoco*(-loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(1,j)*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+loc_z/(exp((doxmin-kmax1*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1
     +*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))-loc_z/(exp((doxmi
     +n-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin
     +)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,
     +j))-loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j))*(-kmax1/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/doxmin)+loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))-loc_z/(exp((doxmin-kmax1*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)
     +*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))-fdocn*(-loc_z
     +/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*sp(8,j)*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j))*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksn
     +doc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*km
     +ax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksn
     +doc+sp(3,j))-loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax
     +2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j))**2-loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +)+1)**2*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3m
     +in*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+loc_z/(exp((no3min-kmax2
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))-loc_z/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)-fdocs*(-loc_z/(e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*
     +*2*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)+loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/
     +st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-loc_z/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kma
     +x3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-loc_z/(exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*
     +sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc
     +_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-loc_z/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,
     +j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
            endif
            if (isp.eq.4) then
            rat = -foo*(loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+loc_z/(exp((doxmin-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp
     +(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))-foa*(loc_z*
     +sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/a
     +mmin)+1)+loc_z*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(kso
     +x+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/ammin)+1))+loc_z*(250-sp(4,j))*(1-1/(exp(0.25D1*(0.4
     +D0-ws)/st)+1))
            drdc = -foo*(-loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(1,j)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+loc_z/(exp((doxmin-kmax1*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*s
     +p(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))-loc_z/(exp((doxmin-
     +kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+
     +1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**
     +2-loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))*(-kmax1*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4
     +,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/
     +doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j)))/st/doxmin)+loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ks
     +odoc+sp(3,j))/(ksox+sp(4,j))-loc_z/(exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)-foa*(loc_z*sp(2
     +0,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))/(exp((ammin-kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc
     +_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**
     +2/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/ammin)+1)-loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(
     +ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))
     +**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +ox+sp(4,j)))/st/ammin)+loc_z*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j)
     +)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*sp(21,j)*kmax4*sp(5,j)/(ksamm+
     +sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*sp(21,j)*kma
     +x4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*
     +(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin))-loc_z*(1-
     +1/(exp(0.25D1*(0.4D0-ws)/st)+1))
            endif
            if (isp.eq.5) then
            rat = fcn*loc_z*kpd*sp(18,j)-fcn*(loc_z/(exp((amming-sp(5,j)
     +)/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)+loc_z/(exp((amming-sp(5,j))/st/ammi
     +ng)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j)))/st/doxmin)+1)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn
     +*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
     ++loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(9,j)*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z/(exp((amming-s
     +p(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14
     +,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)+1)+loc_z/(exp((amming-sp(5,j)
     +)/st/amming)+1)*Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/ammin)+1)+loc_z/(exp((amming-sp(5,j))/st/amming)
     ++1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j
     +))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)
     +))/st/ammin)+1))-faa*(loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))
     +*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+loc_z*sp(21,j)*kmax4*sp(5,j)/
     +(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(
     +ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
            drdc = -fcn*(loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Yo
     +*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(ex
     +p((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+loc_z/(exp(
     +(amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/st/amming*exp((am
     +ming-sp(5,j))/st/amming)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)
     +**2*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+loc_z/(exp((ammin
     +g-sp(5,j))/st/amming)+1)**2*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)+1)/st/amming*exp((amming-sp(5,j))/st/a
     +mming)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ys*kmax3*sp(13
     +,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/st/amming*exp((amming-sp(5
     +,j))/st/amming)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ys*km
     +ax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/st/amming*exp((am
     +ming-sp(5,j))/st/amming)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)
     +**2*Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j
     +))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)
     +))/st/ammin)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+loc_z/(e
     +xp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax4/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/
     +st/amming)+1)*Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)
     +/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming
     +)+1)*Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,
     +j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))
     ++kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin
     +*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/ammin)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ya*sp(21,j)
     +*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
     +/st/amming*exp((amming-sp(5,j))/st/amming)+loc_z/(exp((amming-sp(5
     +,j))/st/amming)+1)*Ya*sp(21,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox
     ++sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +a*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j))
     +/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j
     +)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1
     +)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/
     +(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kma
     +x4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin))-faa*
     +(loc_z*sp(20,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/am
     +min)+1)-loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(k
     +sox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +ox+sp(4,j)))/st/ammin)+1)-loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/
     +(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin)+loc_z*sp(21,j)*kmax4/(ksamm+sp(5
     +,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*sp(21,j)*kmax4*sp(5
     +,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*s
     +p(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/am
     +min)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp
     +(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     +)
            endif
            if (isp.eq.6) then
            rat = -loc_z*(kl*sp(6,j)*vel**0.58D0+kdet*sp(6,j)/(exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)+1))+loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*
     +sp(7,j)+loc_z*kdeac*sp(1,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))-loc_z*kreac*sp(6
     +,j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/doxmin)+1)-0.1D-1*loc_z*km*sp(6,j)
            drdc = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)-loc_z*kreac/(exp((doxmin-kmax
     +1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-0
     +.1D-1*loc_z*km
            endif
            if (isp.eq.7) then
            rat = loc_z*(kl*sp(6,j)*vel**0.58D0+kdet*sp(6,j)/(exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)+1))-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*s
     +p(7,j)+loc_z*kdeac*sp(2,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))-loc_z*kreac*sp(7,
     +j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)-0.1D-1*loc_z*km*sp(7,j)
            drdc = -loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*
     +kreac/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp
     +(4,j)))/st/doxmin)+1)-0.1D-1*loc_z*km
            endif
            if (isp.eq.8) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(
     +8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-loc_z*(kl
     +*sp(8,j)*vel**0.58D0+kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+lo
     +c_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(9,j)-loc_z*kdeac*sp
     +(8,j)*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))+loc_z*kreac
     +*sp(11,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-loc_z*km*sp(8
     +,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-loc_z*(kl*vel**0
     +.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(8,j)/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-loc_z*katt/(
     +exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(1
     +6,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)
     +-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*s
     +p(9,j)-loc_z*kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
     +)-loc_z*km
            endif
            if (isp.eq.9) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(
     +9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z*(kl
     +*sp(8,j)*vel**0.58D0+kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-lo
     +c_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(9,j)-loc_z*kdeac*sp
     +(9,j)*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))+loc_z*kreac
     +*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-loc_z*km*sp(9
     +,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-loc_z*katt*(1-1/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*kdeac*(1-1/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1))-loc_z*km
            endif
            if (isp.eq.10) then
            rat = -fnitra*(loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
     ++1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j))+loc_z/(exp((no3min-kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/n
     +o3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j)))
            drdc = -fnitra*(-loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)**2*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+loc_z/(exp((no3min-kmax2
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))-loc_z/(exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j))-loc_z/(exp((no3min
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*sp(9,j)*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*k
     +indox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+k
     +max2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/
     +(ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
     ++loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*ki
     +ndox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))-lo
     +c_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp
     +(3,j)))
            endif
            if (isp.eq.11) then
            rat = -loc_z*(kl*sp(11,j)*vel**0.58D0+kdet*sp(11,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1))+loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)
     +)*sp(12,j)+loc_z*kdeac*sp(8,j)*(1-1/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1))-loc_z*kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)-0.1D-1*loc_z*km*sp(11,j)
            drdc = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)+kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-loc_z*kreac/(exp((no3min-km
     +ax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksn
     +doc+sp(3,j)))/st/no3min)+1)-0.1D-1*loc_z*km
            endif
            if (isp.eq.12) then
            rat = loc_z*(kl*sp(11,j)*vel**0.58D0+kdet*sp(11,j)/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1))-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
     +*sp(12,j)+loc_z*kdeac*sp(9,j)*(1-1/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1))-loc_z*kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min)+1)-0.1D-1*loc_z*km*sp(12,j)
            drdc = -loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*
     +kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-0.1D-1*loc_z*km
            endif
            if (isp.eq.13) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(
     +13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z*(kl*sp(13,j)*vel**
     +0.58D0+kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z*katt*(1-1
     +/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+1))*sp(14,j)-loc_z*kdeac*sp(13,j)*(1-1/
     +(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1
     +))+loc_z*kreac*sp(16,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)-loc_z*km*sp(13,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z*(kl*vel**0.58D0+kdet/(exp(
     +(Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)
     +-sp(22,j))/st/Bfmax)+1)+kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)
     +**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)-loc_z*k
     +deac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+1))-loc_z*km
            endif
            if (isp.eq.14) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(
     +14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)+loc_z*(kl*sp(13,j)*vel**
     +0.58D0+kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*katt*(1-1
     +/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+1))*sp(14,j)-loc_z*kdeac*sp(14,j)*(1-1/
     +(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1
     +))+loc_z*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)-loc_z*km*sp(14,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z*katt*(1-1/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1))-loc_z*kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1))-loc_z*km
            endif
            if (isp.eq.15) then
            rat = -fsulph*(loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j))+loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(
     +3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))
            drdc = -fsulph*(-loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(13,j)*sp(15,j)/(ksso4+s
     +p(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(
     +kinno3+sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(
     +ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+loc_z/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j))-loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))
     +**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j))-loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(
     +3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j)))/st/so4min)+1)**2*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15
     +,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)+loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j))-loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     +/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)
            endif
            if (isp.eq.16) then
            rat = -loc_z*(kl*sp(16,j)*vel**0.58D0+kdet*sp(16,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1))+loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)
     +)*sp(17,j)+loc_z*kdeac*sp(13,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1))-loc_z*kreac*sp(16,j)/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-0.1
     +D-1*loc_z*km*sp(16,j)
            drdc = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)+kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)-loc_z*kreac/(exp((so4min-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-0.1D-1*loc_z*
     +km
            endif
            if (isp.eq.17) then
            rat = loc_z*(kl*sp(16,j)*vel**0.58D0+kdet*sp(16,j)/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1))-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
     +*sp(17,j)+loc_z*kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1))-loc_z*kreac*sp(17,j)/(exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-0.1D
     +-1*loc_z*km*sp(17,j)
            drdc = -loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*
     +kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc
     ++sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4
     +min)+1)-0.1D-1*loc_z*km
            endif
            if (isp.eq.18) then
            rat = loc_z*km*sp(1,j)+loc_z*km*sp(2,j)+0.1D-1*loc_z*km*sp(6
     +,j)+0.1D-1*loc_z*km*sp(7,j)+loc_z*km*sp(8,j)+loc_z*km*sp(9,j)+0.1D
     +-1*loc_z*km*sp(11,j)+loc_z*km*sp(13,j)+loc_z*km*sp(14,j)+0.1D-1*lo
     +c_z*km*sp(16,j)+0.1D-1*loc_z*km*sp(17,j)-loc_z*kpd*sp(18,j)+loc_z*
     +kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st
     +/amming)+1)
            drdc = -loc_z*kpd
            endif
            if (isp.eq.19) then
            rat = 0
            drdc = 0
            endif
            if (isp.eq.20) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-
     +loc_z*(kl*sp(20,j)*vel**0.58D0+kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1))+loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)-lo
     +c_z*kdeac*sp(20,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))+loc_z*kreac*sp(22,j)/(exp((a
     +mmin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammi
     +n)+1)-loc_z*km*sp(20,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*(k
     +l*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(20,j)/(exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-loc
     +_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j
     +)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st
     +/Bfmax)*sp(21,j)-loc_z*kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm
     ++sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))-loc_z*km
            endif
            if (isp.eq.21) then
            rat = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+
     +loc_z*(kl*sp(20,j)*vel**0.58D0+kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1))-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)-lo
     +c_z*kdeac*sp(21,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))+loc_z*kreac*sp(23,j)/(exp((a
     +mmin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammi
     +n)+1)-loc_z*km*sp(21,j)
            drdc = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*ka
     +tt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*kdeac*(1-1/(exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1
     +))-loc_z*km
            endif
            if (isp.eq.22) then
            rat = -loc_z*(kl*sp(22,j)*vel**0.58D0+kdet*sp(22,j)/(exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax)+1))+loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)
     +)*sp(23,j)+loc_z*kdeac*sp(20,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))-loc_z*kreac*sp(
     +22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/ammin)+1)-0.1D-1*loc_z*km*sp(22,j)
            drdc = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)+kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bf
     +max-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp
     +(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/B
     +fmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)-loc_z*kreac/(exp((ammin-kma
     +x4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.
     +1D-1*loc_z*km
            endif
            if (isp.eq.23) then
            rat = loc_z*(kl*sp(22,j)*vel**0.58D0+kdet*sp(22,j)/(exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)+1))-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
     +*sp(23,j)+loc_z*kdeac*sp(21,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))-loc_z*kreac*sp(2
     +3,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/ammin)+1)-0.1D-1*loc_z*km*sp(23,j)
            drdc = -loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*
     +kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4
     +,j)))/st/ammin)+1)-0.1D-1*loc_z*km
            endif
      end
