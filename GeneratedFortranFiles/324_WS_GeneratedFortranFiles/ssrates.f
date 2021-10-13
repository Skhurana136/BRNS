c      
c     SUBROUTINE ssrates
c      
      subroutine ssrates(rat,drdc,isp,j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            if (isp.eq.1) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax
     +1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-k
     +l*sp(1,j)*vel**0.58D0-kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+ka
     +tt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2,j)-kdeac*sp(1,j)*(1-1/(
     +exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/doxmin)+1))+kreac*sp(6,j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-km*sp(1,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-kl*vel**
     +0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kdet*sp(1,j)/(exp((Bfmax-
     +sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,
     +j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)
     +-kdeac*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1))-km
            endif
            if (isp.eq.2) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax
     +1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+k
     +l*sp(1,j)*vel**0.58D0+kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-ka
     +tt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2,j)-kdeac*sp(2,j)*(1-1/(
     +exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/doxmin)+1))+kreac*sp(7,j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-km*sp(2,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(3,j
     +)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-katt*(1-
     +1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)+1))-kdeac*(1-1/(exp((doxmin-kmax1*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))-km
            endif
            if (isp.eq.3) then
            rat = -fdoco*(1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+1/(exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))-fdocn*(1/(exp((no
     +3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,
     +j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+1/(exp((
     +no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))-fdocs
     +*(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4mi
     +n)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+1/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*s
     +p(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))+kpd*sp(18,j)
            drdc = -fdoco*(-1/((exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*sp(1,j)*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1
     +*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))-1/((ex
     +p((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1)**2)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))*(-kmax1/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp
     +((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st
     +/doxmin)+1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j))-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(
     +3,j))**2*sp(4,j)/(ksox+sp(4,j)))-fdocn*(-1/((exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1)**2)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))-1/(exp((no3mi
     +n-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2-1/((exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2)*kmax2*sp(9,j)*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(
     +-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+s
     +p(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,
     +j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,
     +j))-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j))**2)-fdocs*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+1)**2)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+s
     +p(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/s
     +t/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-1
     +/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
     ++1)**2)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/
     +st/so4min)+1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/
     +(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
     +/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-1/(exp((so
     +4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*
     +sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
            endif
            if (isp.eq.4) then
            rat = -foo*(1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j))+1/(exp((doxmin-kmax1*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp
     +(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))-foa*(sp(20,j)*kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*
     +sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+sp(21
     +,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     ++1))+(250-sp(4,j))*(1-1/(exp(0.25D1*(0.4D0-ws)/st)+1))
            drdc = -foo*(-1/((exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2)*kmax1*sp(1,j)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*sp(3,j)/(ksodoc+sp(
     +3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1/(exp((doxmin-kmax1*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)
     +*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))-1/(exp((doxmin-kmax1*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*s
     +p(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2-1/((exp(
     +(doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +doxmin)+1)**2)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))*(-kmax1*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((
     +doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))/(k
     +sox+sp(4,j))-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j))**2)-foa*(sp(20,j)*kmax4*sp(5,j)/(ksamm+
     +sp(5,j))/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp(20,j)*kmax4*sp(5,j)/(ksamm
     ++sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp(20,j)*kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kma
     +x4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(
     +5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)
     +/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+sp(21,j)*kmax4*
     +sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(
     +ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp(21,j)*kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((ammin-kma
     +x4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp
     +(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((
     +ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/amm
     +in)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammi
     +n-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin))
     +-1+1/(exp(0.25D1*(0.4D0-ws)/st)+1)
            endif
            if (isp.eq.5) then
            rat = fcn*kpd*sp(18,j)-fcn*(1/(exp((amming-sp(5,j))/st/ammin
     +g)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4
     +,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/doxmin)+1)+1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax
     +1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((do
     +xmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/dox
     +min)+1)+1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+1/(exp((amming-sp
     +(5,j))/st/amming)+1)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)+1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*k
     +max3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)+1/(exp((amming-s
     +p(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)+1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax4*
     +sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+1/(exp
     +((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+s
     +p(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))-faa*(sp(20,j)*kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+sp(21,j
     +)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1
     +))
            drdc = -fcn*(1/((exp((amming-sp(5,j))/st/amming)+1)**2)*Yo*k
     +max1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp(
     +(doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +doxmin)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+1/((exp((ammi
     +ng-sp(5,j))/st/amming)+1)**2)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/st/amming*exp((amming
     +-sp(5,j))/st/amming)+1/((exp((amming-sp(5,j))/st/amming)+1)**2)*Yn
     +*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
     +/st/amming*exp((amming-sp(5,j))/st/amming)+1/((exp((amming-sp(5,j)
     +)/st/amming)+1)**2)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+1
     +/((exp((amming-sp(5,j))/st/amming)+1)**2)*Ys*kmax3*sp(13,j)*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)/st/amming*exp((amming-sp(5,j))/st/am
     +ming)+1/((exp((amming-sp(5,j))/st/amming)+1)**2)*Ys*kmax3*sp(14,j)
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)/st/amming*exp((amming-sp(5,j)
     +)/st/amming)+1/((exp((amming-sp(5,j))/st/amming)+1)**2)*Ya*sp(20,j
     +)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1
     +)/st/amming*exp((amming-sp(5,j))/st/amming)+1/(exp((amming-sp(5,j)
     +)/st/amming)+1)*Ya*sp(20,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp
     +(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/ammin)+1)-1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20
     +,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((
     +ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/amm
     +in)+1)-1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5
     +,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax
     +4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5
     +,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/
     +(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1/((exp((amming-
     +sp(5,j))/st/amming)+1)**2)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j
     +))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)/st/amming*exp((amming-sp(5,
     +j))/st/amming)+1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*k
     +max4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5
     +,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-1/(exp((a
     +mming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5
     +,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-1/(exp((amming-sp(5,j)
     +)/st/amming)+1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ks
     +ox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)
     +))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/ammin))-faa*(sp(20,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/ammin)+1)-sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))
     +**2*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp(20,j)*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm
     ++sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+sp(21,j)*kmax4/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(
     +5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp(21,j)*kmax4*sp(5,j)/
     +(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-sp(21,j)*km
     +ax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kma
     +x4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2
     +*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin))
            endif
            if (isp.eq.6) then
            rat = -kl*sp(6,j)*vel**0.58D0-kdet*sp(6,j)/(exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,j)+kdeac*
     +sp(1,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/doxmin)+1))-kreac*sp(6,j)/(exp((doxmin-kmax1*sp(
     +3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-0.1D-1
     +*km*sp(6,j)
            drdc = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kde
     +t*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)*sp(7,j)-kreac/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-0.1D-1*km
            endif
            if (isp.eq.7) then
            rat = kl*sp(6,j)*vel**0.58D0+kdet*sp(6,j)/(exp((Bfmax-sp(1,j
     +)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st
     +/Bfmax)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,j)+kdeac*s
     +p(2,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1))-kreac*sp(7,j)/(exp((doxmin-kmax1*sp(3
     +,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-0.1D-1*
     +km*sp(7,j)
            drdc = -katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-kreac/(exp((
     +doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin)+1)-0.1D-1*km
            endif
            if (isp.eq.8) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-kl*sp(8,j)*ve
     +l**0.58D0-kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+katt*(1-1/(exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+1))*sp(9,j)-kdeac*sp(8,j)*(1-1/(exp((no3min-
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(k
     +sndoc+sp(3,j)))/st/no3min)+1))+kreac*sp(11,j)/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)-km*sp(8,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-kl*vel**0.58D0-kdet/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)-kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8
     +,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax
     +)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-kdeac*(1-1/
     +(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))-km
            endif
            if (isp.eq.9) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(9,j)
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+kl*sp(8,j)*ve
     +l**0.58D0+kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-katt*(1-1/(exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+1))*sp(9,j)-kdeac*sp(9,j)*(1-1/(exp((no3min-
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(k
     +sndoc+sp(3,j)))/st/no3min)+1))+kreac*sp(12,j)/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)-km*sp(9,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-katt*(1-1/(exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)+1))-kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1))-km
            endif
            if (isp.eq.10) then
            rat = -fnitra*(1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*
     +kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))+1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1
     +)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))
            drdc = -fnitra*(-1/((exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)**2)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10
     +,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4,j))/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min
     +*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j
     +))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp((no3min-kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))-1/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j))-1/((exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)+1)**2)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindo
     +x+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/
     +(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3
     +,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/(exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(
     +4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))-1/(exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))
            endif
            if (isp.eq.11) then
            rat = -kl*sp(11,j)*vel**0.58D0-kdet*sp(11,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(12,j)+kde
     +ac*sp(8,j)*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))-kreac*
     +sp(11,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-0.1D-1*km*sp(1
     +1,j)
            drdc = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kde
     +t*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(12,j)-kreac/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1)-0.1D-1*km
            endif
            if (isp.eq.12) then
            rat = kl*sp(11,j)*vel**0.58D0+kdet*sp(11,j)/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(12,j)+kdea
     +c*sp(9,j)*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))-kreac*s
     +p(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-0.1D-1*km*sp(12
     +,j)
            drdc = -katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-kreac/(exp((
     +no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-0.1D-1*km
            endif
            if (isp.eq.13) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j
     +)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)-kl*sp(13,j)*vel**0.58D0-kdet
     +*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax)+1))*sp(14,j)-kdeac*sp(13,j)*(1-1/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))+kreac*sp(16,j)/(exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-km*s
     +p(13,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)-kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1)-kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfma
     +x*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)-kdeac*(1-1/(exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))-km
            endif
            if (isp.eq.14) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j
     +)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)+kl*sp(13,j)*vel**0.58D0+kdet
     +*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-katt*(1-1/(exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax)+1))*sp(14,j)-kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))+kreac*sp(17,j)/(exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-km*s
     +p(14,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1))-kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1))-km
            endif
            if (isp.eq.15) then
            rat = -fsulph*(1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
            drdc = -fsulph*(-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)+1)**2)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ks
     +so4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)/(ksso4+sp(15,
     +j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno
     +3+sp(10,j))-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +)-1/((exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+
     +sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4m
     +in)+1)**2)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))*(-kmax
     +3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)/st/so4min)+1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)*kmax3*sp(14,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+s
     +p(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-1/(exp((s
     +o4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3
     +*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))
            endif
            if (isp.eq.16) then
            rat = -kl*sp(16,j)*vel**0.58D0-kdet*sp(16,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(17,j)+kde
     +ac*sp(13,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)+1))-kreac*sp(16,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)-0.1D-1*km*sp(16,j)
            drdc = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kde
     +t*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(17,j)-kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)-0.1D-1*km
            endif
            if (isp.eq.17) then
            rat = kl*sp(16,j)*vel**0.58D0+kdet*sp(16,j)/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(17,j)+kdea
     +c*sp(14,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1))-kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-0.1D-1*km*sp(17,j)
            drdc = -katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-kreac/(exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-0.1D
     +-1*km
            endif
            if (isp.eq.18) then
            rat = km*sp(1,j)+km*sp(2,j)+0.1D-1*km*sp(6,j)+0.1D-1*km*sp(7
     +,j)+km*sp(8,j)+km*sp(9,j)+0.1D-1*km*sp(11,j)+km*sp(13,j)+km*sp(14,
     +j)+0.1D-1*km*sp(16,j)+0.1D-1*km*sp(17,j)-kpd*sp(18,j)+kmax5*sp(5,j
     +)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)
            drdc = -kpd
            endif
            if (isp.eq.19) then
            rat = 0
            drdc = 0
            endif
            if (isp.eq.20) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-kl*s
     +p(20,j)*vel**0.58D0-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kat
     +t*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)-kdeac*sp(20,j)*(1-1/
     +(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/ammin)+1))+kreac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-km*sp(20,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)
     +/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-kl*vel**0.58
     +D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kdet*sp(20,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-katt/(exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)-k
     +deac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/ammin)+1))-km
            endif
            if (isp.eq.21) then
            rat = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+kl*s
     +p(20,j)*vel**0.58D0+kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kat
     +t*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)-kdeac*sp(21,j)*(1-1/
     +(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/ammin)+1))+kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-km*sp(21,j)
            drdc = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)
     +/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-katt*(1-1/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1))-kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(
     +ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))-km
            endif
            if (isp.eq.22) then
            rat = -kl*sp(22,j)*vel**0.58D0-kdet*sp(22,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(23,j)+kde
     +ac*sp(20,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/ammin)+1))-kreac*sp(22,j)/(exp((ammin-kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.1D-1*
     +km*sp(22,j)
            drdc = -kl*vel**0.58D0-kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-kde
     +t*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)-katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(23,j)-kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.1D-1*km
            endif
            if (isp.eq.23) then
            rat = kl*sp(22,j)*vel**0.58D0+kdet*sp(22,j)/(exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(23,j)+kdea
     +c*sp(21,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/ammin)+1))-kreac*sp(23,j)/(exp((ammin-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.1D-1*k
     +m*sp(23,j)
            drdc = -katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20
     +,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-kreac/(exp((
     +ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/amm
     +in)+1)-0.1D-1*km
            endif
      end
