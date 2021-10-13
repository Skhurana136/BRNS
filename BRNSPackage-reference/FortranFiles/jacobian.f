c      
c     SUBROUTINE jacobian
c      
      subroutine jacobian(pd,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension pd(ncomp,ncomp)
        call switches(j)
         pd(20,23) = -loc_z*kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.1D-1*loc_z*km-1/delt
         pd(17,8) = 0
         pd(4,19) = 0
         pd(7,5) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1
     +*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((dox
     +min-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)/st/amming*exp((amming-sp(5,j))/st/amming)
         pd(8,15) = 0
         pd(13,7) = 0
         pd(16,15) = -loc_z*kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z
     +*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,
     +j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(20,4) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j
     +)*kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z/
     +(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(k
     +samm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z/(exp((amm
     +ing-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j
     +))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox
     ++sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/ammin)+loc_z*kdeac*sp(21,j)/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)*
     +*2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(k
     +samm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+loc_z*k
     +reac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(kso
     +x+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2
     +)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/ammin)
         pd(1,6) = loc_z*kreac/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+0.1D-1*loc_z*km+loc_z*(kl
     +*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(6,j)/(exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))+loc_z
     +*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)*sp(7,j)+1/delt
         pd(4,10) = -foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/
     +st/no3min)+1)**2*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(
     +4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))
     +/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+foo/fdoco*fdocn*
     +loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kin
     +dox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))-foo
     +/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax
     +2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(
     +3,j)/(ksndoc+sp(3,j))-foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1)**2*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kin
     +dox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+foo/fdoco
     +*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9
     +,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j))-foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +**2*sp(3,j)/(ksndoc+sp(3,j))+foo/fdoco/fnitra*fdocn/delt
         pd(7,14) = 0
         pd(9,18) = 1/fcn*faa*foo/fdoco/foa*loc_z*kpd-1/delt
         pd(15,17) = 0
         pd(20,17) = 0
         pd(1,3) = -loc_z*kreac*sp(6,j)/(exp((doxmin-kmax1*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +*2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-loc_z*kdeac*sp(1,
     +j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)**2*(-kmax1/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/doxmin)
         pd(1,19) = 0
         pd(17,14) = loc_z*kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1))
         pd(21,21) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*kmax4*
     +sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*
     +katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-loc_z*kdeac*(1-1/(exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     ++1))-loc_z*km-1/delt
         pd(4,13) = 0
         pd(10,7) = 0
         pd(16,21) = 0
         pd(4,4) = -foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1)**2*kmax2**2*sp(8,j)*kindox**2/(kindox+sp(4,j))**3*sp(
     +10,j)**2/(ksno3+sp(10,j))**2*sp(3,j)**2/(ksndoc+sp(3,j))**2/st/no3
     +min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-foo/fdoco*fdocn*loc_z/(
     +exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kindox/(ki
     +ndox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)
     +)-foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))
     +*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
     +**2*kmax2**2*sp(9,j)*kindox**2/(kindox+sp(4,j))**3*sp(10,j)**2/(ks
     +no3+sp(10,j))**2*sp(3,j)**2/(ksndoc+sp(3,j))**2/st/no3min*exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)-foo/fdoco*fdocn*loc_z/(exp((no3min-
     +kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(k
     +sndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j)
     +)**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
         pd(20,11) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(1,13) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(14,11) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(21,15) = 0
         pd(17,20) = 0
         pd(19,21) = 0
         pd(22,7) = 0
         pd(2,4) = -loc_z*kdeac*sp(2,j)/(exp((doxmin-kmax1*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-loc_z*kreac*sp(7,
     +j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,
     +j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/d
     +oxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/doxmin)
         pd(12,7) = 1/delt
         pd(2,12) = 0
         pd(10,14) = fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j))
         pd(20,5) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ya*sp(2
     +0,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((am
     +min-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin
     +)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+loc_z/(exp((amming-
     +sp(5,j))/st/amming)+1)*Ya*sp(20,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming)+
     +1)*Ya*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4
     +,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(
     +20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((a
     +mmin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammi
     +n)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5
     +,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+l
     +oc_z*kdeac*sp(21,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp
     +(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/ammin)+loc_z*kreac*sp(23,j)/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-k
     +max4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+s
     +p(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,
     +j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
         pd(1,7) = -loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)
     +-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
         pd(11,2) = -1/delt
         pd(17,2) = 0
         pd(21,9) = 0
         pd(11,20) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-loc_z*km-1
     +/delt
         pd(21,3) = 0
         pd(23,2) = 0
         pd(11,7) = -1/delt
         pd(18,9) = -1/delt
         pd(5,22) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)*
     +*2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-0.1D-1
     +*fcn*loc_z*km-fcn/delt
         pd(6,19) = 0
         pd(9,6) = -1/delt
         pd(13,2) = 0
         pd(2,6) = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(5,5) = -fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Yn*
     +kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/
     +st/amming*exp((amming-sp(5,j))/st/amming)-fcn*loc_z/(exp((amming-s
     +p(5,j))/st/amming)+1)**2*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/st/amming*exp((amming-sp(5
     +,j))/st/amming)-fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Y
     +o*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(e
     +xp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/doxmin)+1)/st/amming*exp((amming-sp(5,j))/st/amming)-1/delt
         pd(12,12) = 0
         pd(15,22) = 0
         pd(18,22) = -0.1D-1*loc_z*km-1/delt
         pd(2,1) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(3,17) = 0
         pd(22,9) = 0
         pd(14,16) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(18,5) = loc_z*kmax5/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming
     +-sp(5,j))/st/amming)+1)-loc_z*kmax5*sp(5,j)/(ksamm+sp(5,j))**2/(ex
     +p(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)+0.1D2*loc_z*kmax5*sp(
     +5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+
     +1)**2/st/amming*exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)-1/fcn/
     +delt
         pd(10,2) = 0
         pd(12,16) = -0.1D-1*loc_z*km
         pd(13,16) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)-0.1D-1*loc
     +_z*km-1/delt
         pd(10,10) = -fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(13,j)*sp
     +(15,j)**2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindo
     +x**2/(kindox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/st/so4min*
     +exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j
     +))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-fn
     +itra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))**2-fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(14,j)*sp(15,j
     +)**2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/
     +(kindox+sp(4,j))**2*kinno3**2/(kinno3+sp(10,j))**3/st/so4min*exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-fnitra/
     +fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(
     +3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))**2
         pd(3,23) = 0
         pd(6,16) = -loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*kdet*sp(6,j)/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp
     +(2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(7,j)
         pd(9,1) = -1/delt
         pd(10,19) = 0
         pd(6,7) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
         pd(18,21) = -loc_z*km-1/delt
         pd(8,11) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(22,2) = 0
         pd(2,8) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(6,10) = 0
         pd(22,21) = -loc_z*kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+
     +sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
         pd(23,18) = 0
         pd(15,7) = 0
         pd(18,15) = 1/fcn*faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(3,16) = 0
         pd(6,2) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))
         pd(16,2) = 0
         pd(19,20) = 0
         pd(22,8) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(11,12) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z
     +*kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
         pd(14,15) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*
     +sp(14,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((amming-sp(5,j))/st
     +/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z*kdeac
     +*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ks
     +so4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)-loc_z*kreac*sp(17,j)/(exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax
     +3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)/st/so4min)
         pd(1,4) = -loc_z*kreac*sp(6,j)/(exp((doxmin-kmax1*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-loc_z*kdeac*sp(1,
     +j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,
     +j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/d
     +oxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/doxmin)
         pd(6,20) = -loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*kdet*sp(6,j)/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp
     +(2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(7,j)
         pd(9,5) = loc_z*kmax5/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-
     +sp(5,j))/st/amming)+1)-loc_z*kmax5*sp(5,j)/(ksamm+sp(5,j))**2/(exp
     +(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)+0.1D2*loc_z*kmax5*sp(5
     +,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1
     +)**2/st/amming*exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)-1/fcn/d
     +elt
         pd(3,6) = 0
         pd(12,11) = -0.1D-1*loc_z*km
         pd(23,7) = 0
         pd(22,14) = 0
         pd(11,6) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-1/delt
         pd(14,21) = 0
         pd(9,11) = loc_z*kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)+0.1D-1*loc_z*km
         pd(23,1) = 0
         pd(16,14) = -loc_z*kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min)+1))
         pd(19,8) = 0
         pd(22,20) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(23,19) = 0
         pd(15,8) = 0
         pd(18,16) = -1/delt
         pd(13,11) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(2,13) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(6,15) = 0
         pd(9,2) = -1/delt
         pd(23,13) = 0
         pd(3,12) = 0
         pd(6,6) = -0.1D-1*loc_z*km-loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,
     +j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*(kl*vel**0.
     +58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-loc_z*katt/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp
     +(2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(7,j)-1/delt
         pd(15,2) = 0
         pd(8,6) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(10,3) = -fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(13,j)*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+fnitra/fdocn*fdocs
     +*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-fnitra/fdocn*fdo
     +cs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/
     +so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-fni
     +tra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)**2*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j)
     +)*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+
     +sp(10,j))*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+fnitra/fdocn*fdocs*loc_z/(exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(
     +14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))-fnitra/fdocn*fdocs*loc_z/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*s
     +p(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
         pd(13,17) = -loc_z*kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)-0.1D-1*loc_z*km-1/delt
         pd(5,13) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)*
     +*2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*lo
     +c_z*km-fcn/delt
         pd(7,1) = loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)+kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax))+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+loc_z*katt/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*
     +sp(7,j)
         pd(8,19) = 0
         pd(10,6) = 0
         pd(7,9) = 0
         pd(9,14) = -1/delt
         pd(10,23) = 0
         pd(2,10) = 0
         pd(21,22) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(4,18) = 0
         pd(7,6) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*(kl*vel**0.58D0+kdet/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax))+loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)+loc_z*katt/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*
     +sp(7,j)
         pd(16,16) = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1)+kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(4,9) = foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)*kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j))
         pd(7,15) = 0
         pd(9,19) = 0
         pd(19,15) = 0
         pd(20,22) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)-0.1D-1*loc
     +_z*km-1/delt
         pd(3,2) = -fdoco*loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j))
         pd(8,22) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(17,9) = 0
         pd(12,2) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+1/
     +delt
         pd(11,3) = -loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-k
     +max2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j)))/st/no3min)-loc_z*kreac*sp(12,j)/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*e
     +xp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-1/fcn*faa*foo/fdoco/foa/delt
         pd(14,6) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(21,10) = 0
         pd(19,22) = 0
         pd(5,9) = -fcn*loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-fc
     +n*loc_z*kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))-fcn
     +*loc_z*km-fcn/delt
         pd(12,8) = -loc_z*km
         pd(16,9) = 0
         pd(19,3) = 0
         pd(20,10) = 0
         pd(1,12) = 0
         pd(14,12) = 0
         pd(6,17) = 0
         pd(19,4) = 0
         pd(23,15) = 0
         pd(15,4) = -fsulph/fdocs/foo*fdoco*foa*loc_z*sp(20,j)*kmax4*sp(
     +5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*
     +fdoco*foa*loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +ox+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*loc_z*sp(20
     +,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     ++1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j
     +)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-fsu
     +lph/fdocs/foo*fdoco*foa*loc_z*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j
     +))/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*loc_z*sp
     +(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2/(ex
     +p((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +ammin)+1)+fsulph/fdocs/foo*fdoco*foa*loc_z*sp(21,j)*kmax4*sp(5,j)/
     +(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(
     +ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(
     +5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-fsulph/fdocs/foo*fdoc
     +o*loc_z*(1-1/(exp(0.25D1*(0.4D0-ws)/st)+1))-fsulph/fdocs/foo*fdoco
     +/delt
         pd(18,20) = -loc_z*km-1/delt
         pd(10,1) = 0
         pd(12,17) = -0.1D-1*loc_z*km
         pd(5,10) = -fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*km
     +ax2*sp(8,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+fcn*loc_z/(
     +exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j))/(e
     +xp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+fcn*loc_z/(exp((amming-sp
     +(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*ex
     +p((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*loc_z*kdeac*sp(9,j)/(exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp
     +(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*loc_z*kreac
     +*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*ki
     +ndox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+km
     +ax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(
     +ksndoc+sp(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +faa*foo/fdoco/foa/fnitra*fdocn/delt
         pd(6,11) = -loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*kdet*sp(6,j)/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp
     +(2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(7,j)
         pd(18,6) = -1/delt
         pd(19,10) = 0
         pd(7,23) = 0
         pd(12,19) = 0
         pd(14,5) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ys*kmax
     +3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/st/amming*exp((ammi
     +ng-sp(5,j))/st/amming)
         pd(18,14) = -1/delt
         pd(8,10) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*
     +sp(9,j)*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z/(exp((amm
     +ing-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z/(exp((amming-sp(5,j))/st/am
     +ming)+1)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-km
     +ax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksn
     +doc+sp(3,j)))/st/no3min)+loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)+loc_z*kreac*sp(12,j)/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j)
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
         pd(13,13) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z*(kl*vel**0.58D0+kdet/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)-loc_z
     +*km-1/delt
         pd(22,1) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(5,2) = -fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kma
     +x1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-km
     +ax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
         pd(5,15) = faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(11,23) = -0.1D-1*loc_z*km-1/delt
         pd(12,1) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+1/
     +delt
         pd(8,21) = 0
         pd(11,8) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-1/delt
         pd(4,2) = 0
         pd(9,9) = -loc_z*kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1))-1/delt
         pd(10,11) = 0
         pd(11,19) = 0
         pd(21,2) = 0
         pd(2,23) = 0
         pd(3,18) = loc_z*kpd
         pd(16,4) = -loc_z*kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/
     +(kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)-loc_z*kreac*sp(17,j)/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(22,10) = 0
         pd(2,3) = -loc_z*kdeac*sp(2,j)/(exp((doxmin-kmax1*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +*2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-loc_z*kreac*sp(7,
     +j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin)+1)**2*(-kmax1/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,
     +j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/d
     +oxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/doxmin)
         pd(14,17) = loc_z*kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+1)
         pd(6,18) = 0
         pd(8,20) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(13,1) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(5,4) = fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax
     +2*sp(8,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(
     +3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+f
     +cn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2**2*sp(8,j)*k
     +indox**2/(kindox+sp(4,j))**3*sp(10,j)**2/(ksno3+sp(10,j))**2*sp(3,
     +j)**2/(ksndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)**2/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*loc_z*k
     +deac*sp(9,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*k
     +indox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndo
     +c+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*lo
     +c_z*kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*k
     +max2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/
     +(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-
     +fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(
     +3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(k
     +sodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(ex
     +p((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(exp(
     +(amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j)/(
     +ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fcn*loc_z/(exp((ammin
     +g-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(exp((amming-sp(5,j))/st
     +/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(kso
     +x+sp(4,j))**2/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(exp((amming-sp(5,j))/st/a
     +mming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+
     +sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j)/(ksodoc+sp(3,j))/(ksox
     ++sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2
     +)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j)))/st/doxmin)+faa/foa*loc_z*(1-1/(exp(0.25D1*(0.4D0-ws)/
     +st)+1))+faa/foa/delt
         pd(12,13) = -loc_z*km
         pd(23,5) = loc_z*kdeac*sp(20,j)/(exp((ammin-kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+
     +sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*s
     +p(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+loc_z*kdeac*sp(21,j)/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/a
     +mmin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*s
     +p(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((am
     +min-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin
     +)+loc_z*kreac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(
     +4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox
     ++sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/ammin)+loc_z*kreac*sp(23,j)/(exp((ammin-kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*
     +(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksam
     +m+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
         pd(8,5) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Yn*kmax
     +2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)/st/a
     +mming*exp((amming-sp(5,j))/st/amming)
         pd(17,21) = 0
         pd(12,6) = 1/delt
         pd(16,7) = 0
         pd(20,12) = 0
         pd(1,14) = 0
         pd(14,10) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3
     +*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((ammin
     +g-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+s
     +p(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,
     +j))**2*kinno3**2/(kinno3+sp(10,j))**3/(exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/st/so4min*exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z*kdeac*sp(
     +14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc
     ++sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4
     +min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*ex
     +p((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_
     +z*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/s
     +t/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)
         pd(10,13) = fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j))
         pd(16,13) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(20,6) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(1,8) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(11,1) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-1/delt
         pd(4,20) = 0
         pd(8,16) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(20,3) = 0
         pd(15,16) = 0
         pd(19,13) = 0
         pd(22,15) = 0
         pd(14,22) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(17,7) = 0
         pd(4,14) = 0
         pd(10,8) = 0
         pd(3,5) = 0
         pd(7,11) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(
     +2,j)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(7,j)
         pd(19,19) = 1
         pd(1,2) = 0
         pd(11,13) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-1/delt
         pd(17,13) = loc_z*kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4
     ++sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min)+1))
         pd(21,20) = loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8
     +,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax
     +)+1)+kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax))+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(10,12) = 0
         pd(20,7) = 0
         pd(1,9) = 0
         pd(21,11) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(18,23) = -0.1D-1*loc_z*km-1/delt
         pd(3,22) = 0
         pd(10,18) = 0
         pd(16,8) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(17,22) = 0
         pd(19,23) = 0
         pd(17,12) = 0
         pd(4,15) = 0
         pd(5,21) = -fcn*loc_z*km-fcn/delt
         pd(7,10) = 0
         pd(15,21) = -fsulph/fdocs/foo*fdoco*foa*loc_z*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
         pd(19,14) = 0
         pd(20,21) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z
     +*kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(kso
     +x+sp(4,j)))/st/ammin)+1))
         pd(1,23) = 0
         pd(18,10) = 1/fcn*faa*foo/fdoco/foa/fnitra*fdocn/delt
         pd(16,17) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z
     +*kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)+1)+0.1D-1*loc_z*km+1/delt
         pd(20,2) = 0
         pd(2,15) = 0
         pd(4,8) = foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)*kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j)
     +)*sp(3,j)/(ksndoc+sp(3,j))
         pd(7,16) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(
     +2,j)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(7,j)
         pd(15,15) = -1/delt
         pd(9,15) = 1/fcn*faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(13,12) = 0
         pd(5,14) = -fcn*loc_z*km-fcn/delt
         pd(14,23) = 0
         pd(6,12) = 0
         pd(19,9) = 0
         pd(22,19) = 0
         pd(23,20) = loc_z*kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
         pd(2,14) = 0
         pd(7,22) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(
     +2,j)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(7,j)
         pd(12,18) = loc_z*kpd+1/delt
         pd(15,9) = 0
         pd(17,3) = loc_z*kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**
     +2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*
     +kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(16,j)/(exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*
     +(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+loc_z*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4mi
     +n*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(1,5) = 0
         pd(2,17) = 0
         pd(9,20) = -loc_z*km-1/delt
         pd(13,18) = 0
         pd(22,6) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(9,3) = -loc_z*kdeac*sp(8,j)/(exp((no3min-kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-km
     +ax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksn
     +doc+sp(3,j)))/st/no3min)-loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp
     +(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)-loc_z*kreac*sp(11,j)/(exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/s
     +t/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-loc_z*kreac*sp(12,
     +j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-1/fcn*
     +faa*foo/fdoco/foa/delt
         pd(23,14) = 0
         pd(6,5) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Yo*kmax1
     +*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((dox
     +min-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in)+1)/st/amming*exp((amming-sp(5,j))/st/amming)
         pd(12,23) = 0
         pd(15,3) = fsulph/fdocs/delt
         pd(18,19) = 0
         pd(21,16) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(2,16) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(9,4) = -loc_z*kdeac*sp(8,j)/(exp((no3min-kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindox
     +/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j
     +)))/st/no3min)-loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j)))/st/no3min)-loc_z*kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/
     +st/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)-loc_z*kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(
     +kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +)/st/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1/fcn*faa/foa*loc_z*(1-1/(exp(0.25D1*(0.4D0-ws
     +)/st)+1))+1/fcn*faa/foa/delt
         pd(16,22) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(4,3) = -foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1)**2*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(kindox+sp(4
     +,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)
     +/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+foo/fdoco*fdocn*
     +loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp(8,j)*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))-fo
     +o/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kma
     +x2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,
     +j)/(ksndoc+sp(3,j))**2-foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)+1)**2*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*(-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+foo/fdo
     +co*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2*sp
     +(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp
     +(3,j))-foo/fdoco*fdocn*loc_z/(exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j))**2
         pd(5,3) = -fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kma
     +x2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksnd
     +oc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+fcn*loc_z/(
     +exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox
     ++sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2/(e
     +xp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+fcn*loc_z/(exp((amming-sp
     +(5,j))/st/amming)+1)*Yn*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*
     +sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*e
     +xp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))
     +*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*loc_z*kdeac*sp(9,j)/(exp
     +((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(k
     +indox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))*
     +*2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-fcn*loc_z*kre
     +ac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))
     ++kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(
     +4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3mi
     +n)-fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/
     +(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksod
     +oc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(e
     +xp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**
     +2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fcn*loc_z/(exp((am
     +ming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4
     +,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(exp((amming-sp(5,j))
     +/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j
     +)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/doxmin)+1)+fcn*loc_z/(exp((amming-sp(5,j))/s
     +t/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ks
     +ox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksodoc+sp(3,j))*sp(4,j)/(k
     +sox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,
     +j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/doxmin)-faa*foo/fdoco/foa/delt
         pd(23,6) = 0
         pd(6,1) = -loc_z*km+loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +o*kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxm
     +in-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)-loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet
     +*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,
     +j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/s
     +t/Bfmax))-loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2,j)-loc_z*katt/(exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
     +-1/delt
         pd(16,3) = -loc_z*kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+s
     +p(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +*2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z
     +*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(20,16) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(1,18) = 0
         pd(6,23) = 0
         pd(9,10) = -loc_z*kdeac*sp(8,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-km
     +ax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksn
     +doc+sp(3,j)))/st/no3min)-loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*ki
     +ndox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp(
     +(no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp
     +(3,j)/(ksndoc+sp(3,j)))/st/no3min)-loc_z*kreac*sp(11,j)/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j)
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/
     +no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+s
     +p(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)-loc_z*kreac*sp(12,j)
     +/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kin
     +dox+sp(4,j))/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp
     +(3,j)))/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/fcn*faa
     +*foo/fdoco/foa/fnitra*fdocn/delt
         pd(13,6) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(21,1) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(3,13) = 0
         pd(7,8) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(2
     +,j)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,
     +j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)*sp(7,j)
         pd(16,18) = 0
         pd(3,3) = fdoco*loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(1,j)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fdoco*loc_z/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+fdoco*loc_z/(e
     +xp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(k
     +sox+sp(4,j))+fdoco*loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(2,j)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(k
     +sox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fdoco*loc_z/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+fdoco*loc_z/(e
     +xp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(k
     +sox+sp(4,j))-1/delt
         pd(7,17) = 0
         pd(12,14) = -loc_z*km
         pd(14,18) = 0
         pd(17,11) = 0
         pd(2,7) = loc_z*kreac/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+0.1D-1*loc_z*km+loc_z*kat
     +t*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+1/delt
         pd(10,4) = -fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(13,j)*sp(
     +15,j)**2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox
     +**2/(kindox+sp(4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/st/so4min*e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-fni
     +tra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*s
     +p(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+
     +sp(10,j))-fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3**2*sp(14,j)*sp(15,j)
     +**2/(ksso4+sp(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(
     +kindox+sp(4,j))**3*kinno3**2/(kinno3+sp(10,j))**2/st/so4min*exp((s
     +o4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-fnitra/f
     +docn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10
     +,j))
         pd(3,9) = 0
         pd(17,17) = -loc_z*kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)-0.1D-1*loc_z*km-1/delt
         pd(22,4) = -loc_z*kdeac*sp(21,j)/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j
     +)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-loc_z*kreac*sp(23,j)/(ex
     +p((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*
     +sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((a
     +mmin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammi
     +n)
         pd(5,11) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)*
     +*2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(18,2) = -1/delt
         pd(2,9) = 0
         pd(15,5) = -fsulph/fdocs/foo*fdoco*foa*loc_z*sp(20,j)*kmax4/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*
     +fdoco*foa*loc_z*sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/
     +(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*loc_z*sp(20
     +,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((amm
     +in-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
     ++1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j
     +)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-fsu
     +lph/fdocs/foo*fdoco*foa*loc_z*sp(21,j)*kmax4/(ksamm+sp(5,j))*sp(4,
     +j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/ammin)+1)+fsulph/fdocs/foo*fdoco*foa*loc_z*sp
     +(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j))/(ex
     +p((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +ammin)+1)+fsulph/fdocs/foo*fdoco*foa*loc_z*sp(21,j)*kmax4*sp(5,j)/
     +(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(
     +ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*
     +*2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
         pd(13,14) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z
     +*kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(
     +kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/
     +st/so4min)+1))
         pd(17,23) = 0
         pd(5,16) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)*
     +*2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-0.1D-1
     +*fcn*loc_z*km-fcn/delt
         pd(7,4) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(2,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(exp((doxmin-kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_
     +z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_z/(exp
     +((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j)/
     +(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(10,17) = 0
         pd(3,1) = -fdoco*loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j))
         pd(11,5) = loc_z*kmax5/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming
     +-sp(5,j))/st/amming)+1)-loc_z*kmax5*sp(5,j)/(ksamm+sp(5,j))**2/(ex
     +p(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)+0.1D2*loc_z*kmax5*sp(
     +5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+
     +1)**2/st/amming*exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)-1/fcn/
     +delt
         pd(14,4) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*
     +sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((amming
     +-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(14,j)*sp(15,j)**2/(ksso4+sp
     +(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j
     +))**3*kinno3**2/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/st/so4min*exp((so4min-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z*kdeac*sp(1
     +4,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+
     +sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4m
     +in)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z
     +*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st
     +/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)
         pd(17,5) = 0
         pd(19,12) = 0
         pd(22,16) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(23,23) = -loc_z*kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.1D-1*loc_z*km-1/delt
         pd(9,13) = -1/delt
         pd(10,9) = 0
         pd(13,23) = 0
         pd(21,6) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(3,14) = 0
         pd(19,18) = 0
         pd(11,14) = -1/delt
         pd(14,13) = loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8
     +,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax
     +)+1)+kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((B
     +fmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-s
     +p(22,j))/st/Bfmax))+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/
     +Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j
     +)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(5,19) = 0
         pd(6,22) = -loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*kdet*sp(6,j)/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp
     +(2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(7,j)
         pd(13,5) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ys*kmax
     +3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp
     +(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)/st/amming*exp((ammi
     +ng-sp(5,j))/st/amming)
         pd(5,8) = -fcn*loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(k
     +sno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+fcn*loc_z*(
     +kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(8,j)/(exp
     +((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j
     +)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))+fcn
     +*loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)*sp(9,j)
         pd(7,12) = 0
         pd(12,9) = -loc_z*km
         pd(15,19) = 0
         pd(23,9) = 0
         pd(3,20) = 0
         pd(6,13) = -loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*kdet*sp(6,j)/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp
     +(2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax)*sp(7,j)
         pd(16,6) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(23,11) = 0
         pd(6,4) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(exp((doxmin-kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_
     +z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ks
     +odoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_z/(exp
     +((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j)/
     +(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(8,4) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp
     +(9,j)*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z
     +/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2**2*sp(9,j)*kindox**2
     +/(kindox+sp(4,j))**3*sp(10,j)**2/(ksno3+sp(10,j))**2*sp(3,j)**2/(k
     +sndoc+sp(3,j))**2/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2/st
     +/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+loc_z*kdeac*sp(9,j)
     +/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,
     +j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*kindox/(kindo
     +x+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/s
     +t/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+loc_z*kreac*sp(12,
     +j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*kmax2*kindox/(kin
     +dox+sp(4,j))**2*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))
     +/st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksn
     +o3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
         pd(13,19) = 0
         pd(11,18) = 1/fcn*faa*foo/fdoco/foa*loc_z*kpd-1/delt
         pd(2,5) = 0
         pd(16,12) = 0
         pd(19,6) = 0
         pd(12,22) = 0
         pd(15,10) = -fsulph/fnitra*fdocn/fdocs/delt
         pd(4,5) = 0
         pd(4,16) = 0
         pd(5,20) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)*
     +*2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*lo
     +c_z*km-fcn/delt
         pd(9,8) = -loc_z*kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1))-1/delt
         pd(13,4) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*
     +sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((amming
     +-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(13,j)*sp(15,j)**2/(ksso4+sp
     +(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,j
     +))**3*kinno3**2/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/st/so4min*exp((so4min-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kin
     +dox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kdeac*sp(1
     +4,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+
     +sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4m
     +in)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp
     +((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z
     +*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st
     +/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)
         pd(4,7) = 0
         pd(5,7) = 0
         pd(15,20) = -fsulph/fdocs/foo*fdoco*foa*loc_z*kmax4*sp(5,j)/(ks
     +amm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
         pd(19,17) = 0
         pd(20,20) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*kmax4*
     +sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z*
     +(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+kdet*sp(20,j)/(e
     +xp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16
     +,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax))-l
     +oc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1
     +,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/
     +st/Bfmax)*sp(21,j)-loc_z*km-1/delt
         pd(22,11) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(1,22) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(11,15) = 1/fcn*faa/fsulph*fdocs*foo/fdoco/foa/delt
         pd(15,1) = 0
         pd(18,3) = -1/fcn*faa*foo/fdoco/foa/delt
         pd(13,10) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3
     +*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((ammin
     +g-sp(5,j))/st/amming)+1)*Ys*kmax3**2*sp(13,j)*sp(15,j)**2/(ksso4+s
     +p(15,j))**2*sp(3,j)**2/(kssdoc+sp(3,j))**2*kindox**2/(kindox+sp(4,
     +j))**2*kinno3**2/(kinno3+sp(10,j))**3/(exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2/st/so4min*exp((so4min-k
     +max3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kdeac*sp(
     +14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc
     ++sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4
     +min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*ex
     +p((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_
     +z*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/s
     +t/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)
         pd(20,1) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(2,19) = 0
         pd(15,14) = 0
         pd(21,5) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Ya*sp(2
     +1,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((am
     +min-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin
     +)+1)/st/amming*exp((amming-sp(5,j))/st/amming)+loc_z/(exp((amming-
     +sp(5,j))/st/amming)+1)*Ya*sp(21,j)*kmax4/(ksamm+sp(5,j))*sp(4,j)/(
     +ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming)+
     +1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4
     +,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,
     +j)))/st/ammin)+1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(
     +21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((a
     +mmin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammi
     +n)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5
     +,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin
     +-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-l
     +oc_z*kdeac*sp(21,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j
     +)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp
     +(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/
     +(ksox+sp(4,j)))/st/ammin)-loc_z*kreac*sp(23,j)/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-k
     +max4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+s
     +p(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,
     +j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
         pd(23,4) = loc_z*kdeac*sp(20,j)/(exp((ammin-kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)
     +/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4
     +,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+loc_z*kdeac*sp(21,j)/(exp
     +((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/a
     +mmin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((am
     +min-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin
     +)+loc_z*kreac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*s
     +p(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5
     +,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp
     +(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/ammin)+loc_z*kreac*sp(23,j)/(exp((ammin-kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*
     +(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksam
     +m+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp
     +(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)
         pd(3,10) = 1/fnitra*fdocn/delt
         pd(22,17) = 0
         pd(23,22) = -loc_z*kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,
     +j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-0.1D-1*loc_z*km-1/delt
         pd(11,9) = -loc_z*kdeac*(1-1/(exp((no3min-kmax2*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no
     +3min)+1))-1/delt
         pd(8,9) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z*katt*(1-1
     +/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z*kdeac*(1-1/(exp((no3min-kmax2
     +*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc
     ++sp(3,j)))/st/no3min)+1))+loc_z*km+1/delt
         pd(11,22) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-0.1D-1*loc
     +_z*km-1/delt
         pd(6,8) = -loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(
     +2,j)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(7,j)
         pd(16,11) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(19,5) = 0
         pd(20,8) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(22,23) = loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j
     +)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))+loc_z
     +*kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(
     +4,j)))/st/ammin)+1)+0.1D-1*loc_z*km+1/delt
         pd(23,16) = 0
         pd(1,10) = 0
         pd(18,13) = -1/delt
         pd(2,18) = 0
         pd(8,3) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*s
     +p(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+s
     +p(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3
     ++sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z/(exp((amm
     +ing-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2/(exp((no3m
     +in-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j)))/st/no3min)+1)+loc_z/(exp((amming-sp(5,j))/st/am
     +ming)+1)*Yn*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*exp((no3min-k
     +max2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j)))/st/no3min)+loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*k
     +indox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+s
     +p(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)
     +/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/st/no3min*ex
     +p((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+loc_z*kreac*sp(12,j)/(exp((no
     +3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,
     +j)/(ksndoc+sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))/(ksndoc+sp(3,j))+kmax2*kindox/(kindo
     +x+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))**2)/
     +st/no3min*exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno
     +3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
         pd(3,21) = 0
         pd(23,10) = 0
         pd(21,12) = 0
         pd(17,16) = -loc_z*kreac/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)-0.1D-1*loc_z*km-1/delt
         pd(21,19) = 0
         pd(16,23) = 0
         pd(18,7) = -1/delt
         pd(4,21) = 0
         pd(2,21) = 0
         pd(4,12) = 0
         pd(2,11) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(4,1) = 0
         pd(10,16) = 0
         pd(7,21) = 0
         pd(14,3) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*s
     +p(14,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((amming-sp(5,j))/st
     +/amming)+1)*Ys*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(14,j)*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,
     +j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)-loc_z*kdea
     +c*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/s
     +t/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)-loc_z*kreac*sp(17,j)/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)
         pd(17,4) = loc_z*kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))**2*kinno3/(
     +kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15
     +,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinn
     +o3+sp(10,j)))/st/so4min)+loc_z*kdeac*sp(14,j)/(exp((so4min-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +*2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(16,j)/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(17,j
     +)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(
     +3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
     ++1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))**2*kinno3/(kinno3+sp(10,j))/st/so4min*exp((s
     +o4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(21,7) = 0
         pd(9,16) = -1/delt
         pd(7,3) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(2,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_
     +z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ks
     +odoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_z/(exp
     +((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(8,17) = 0
         pd(10,22) = 0
         pd(20,13) = -loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(1,15) = 0
         pd(14,9) = 0
         pd(18,18) = 1/fcn*faa*foo/fdoco/foa*loc_z*kpd-1/delt
         pd(2,20) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(9,21) = -loc_z*km-1/delt
         pd(17,18) = 0
         pd(22,5) = -loc_z*kdeac*sp(21,j)/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4/(ksamm
     ++sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))**2*
     +sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-loc_z*kreac*sp(23,j)/(ex
     +p((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/
     +ammin)+1)**2*(-kmax4/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))+kmax4*
     +sp(5,j)/(ksamm+sp(5,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/ammin*exp((a
     +mmin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammi
     +n)
         pd(8,12) = -loc_z*kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,
     +j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)
     ++1)
         pd(12,5) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)**2*Yo*kma
     +x1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((d
     +oxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/do
     +xmin)+1)/st/amming*exp((amming-sp(5,j))/st/amming)-loc_z/(exp((amm
     +ing-sp(5,j))/st/amming)+1)**2*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)/st/amming*exp((amming
     +-sp(5,j))/st/amming)-loc_z*kmax5/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0
     +*amming-sp(5,j))/st/amming)+1)+loc_z*kmax5*sp(5,j)/(ksamm+sp(5,j))
     +**2/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)-0.1D2*loc_z*km
     +ax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/a
     +mming)+1)**2/st/amming*exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)
         pd(15,23) = 0
         pd(18,1) = -1/delt
         pd(13,20) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(3,19) = 0
         pd(18,8) = -1/delt
         pd(19,1) = 0
         pd(23,12) = 0
         pd(18,17) = -1/delt
         pd(21,14) = 0
         pd(8,7) = 0
         pd(19,7) = 0
         pd(2,22) = -loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(7,20) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(
     +2,j)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(7,j)
         pd(15,11) = 0
         pd(18,11) = -1/delt
         pd(21,8) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(9,12) = loc_z*kreac/(exp((no3min-kmax2*kindox/(kindox+sp(4,j
     +))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+
     +1)
         pd(11,17) = -1/delt
         pd(13,8) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(4,11) = 0
         pd(13,22) = -loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(11,11) = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1)+kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-1/delt
         pd(5,18) = (fcn*fdoco*foa+faa*foo)/fdoco/foa*loc_z*kpd
         pd(8,1) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(15,18) = -fsulph/fdocs*loc_z*kpd
         pd(23,8) = 0
         pd(16,1) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(20,18) = 0
         pd(22,13) = -loc_z*kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(1,20) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(5,1) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)-fcn*loc
     +_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
         pd(10,20) = 0
         pd(20,15) = 0
         pd(1,17) = 0
         pd(14,7) = 0
         pd(3,4) = fdoco*loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(1,j)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fdoco*loc_z/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+fdoco*loc_z/(e
     +xp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))**2+fdoco*loc_z/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*kmax1*sp(2,j)*sp(3,j)/
     +(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))*(-kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-fdoco*loc_z/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*km
     +ax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+fdoco*loc_z/(e
     +xp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/
     +st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j))**2
         pd(22,3) = 0
         pd(8,14) = 0
         pd(12,3) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+lo
     +c_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(
     +ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+loc_z/(e
     +xp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**
     +2*sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-loc_z/(exp((amming
     +-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)/(ksodoc+sp(3,j))*sp(4,j)/
     +(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)+loc_z/(exp((amming-sp(5,j))/st/ammi
     +ng)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+
     +sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)
     +*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/
     +(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)**2*(-kmax1/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     ++kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j)))/st/doxm
     +in*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)
         pd(9,23) = -0.1D-1*loc_z*km-1/delt
         pd(20,9) = 0
         pd(22,22) = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1)+kdet*sp(22,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((
     +Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-
     +sp(22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(23,j)
         pd(23,17) = 0
         pd(1,11) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(17,6) = 0
         pd(14,19) = 0
         pd(21,23) = loc_z*kreac/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j
     +))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
         pd(4,23) = 0
         pd(5,12) = fcn*loc_z*kreac/(exp((no3min-kmax2*kindox/(kindox+sp
     +(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3m
     +in)+1)
         pd(16,19) = 0
         pd(7,18) = 0
         pd(12,15) = 0
         pd(15,13) = 0
         pd(23,3) = 0
         pd(12,21) = 0
         pd(18,4) = 1/fcn*faa/foa*loc_z*(1-1/(exp(0.25D1*(0.4D0-ws)/st)+
     +1))+1/fcn*faa/foa/delt
         pd(21,17) = 0
         pd(4,17) = 0
         pd(9,7) = -1/delt
         pd(3,8) = 0
         pd(8,18) = 0
         pd(9,22) = -0.1D-1*loc_z*km-1/delt
         pd(10,15) = -fnitra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,
     +j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j
     +))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(13,j)*sp(15
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax
     +3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(ki
     +ndox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+fnitra/fdocn*fdocs
     +*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdo
     +c+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so
     +4min)+1)*kmax3*sp(13,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*
     +kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-fnitra/fdocn*fdoc
     +s*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssd
     +oc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/s
     +o4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kss
     +doc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))-fnit
     +ra/fdocn*fdocs*loc_z/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+1)**2*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j))*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp
     +(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3
     +/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j)))/st/so4min)+fnitra/fdocn*fdocs*loc_z/(exp((so4min-
     +kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14
     +,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,
     +j))*kinno3/(kinno3+sp(10,j))-fnitra/fdocn*fdocs*loc_z/(exp((so4min
     +-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(
     +kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(1
     +4,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+fnitra/fdocn/fsulph*fdoc
     +s/delt
         pd(11,16) = -loc_z*kdet*sp(11,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(12,j)-1/delt
         pd(14,2) = 0
         pd(17,19) = 0
         pd(8,13) = -loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(12,4) = -loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*
     +sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))/(exp((doxmin-kmax1
     +*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+lo
     +c_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(
     +ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((doxmin-kmax1*sp(3,
     +j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+loc_z/(e
     +xp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc
     ++sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1*sp(3,j
     +)/(ksodoc+sp(3,j))/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))*s
     +p(4,j)/(ksox+sp(4,j))**2)/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)-loc_z/(exp((amming
     +-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))/
     +(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j)))/st/doxmin)+1)+loc_z/(exp((amming-sp(5,j))/st/ammi
     +ng)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(
     +4,j))**2/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox
     ++sp(4,j)))/st/doxmin)+1)+loc_z/(exp((amming-sp(5,j))/st/amming)+1)
     +*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/
     +(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)**2*(-kmax1*sp(3,j)/(ksodoc+sp(3,j))/(ksox+sp(4,j))
     ++kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/doxm
     +in*exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)
         pd(6,14) = 0
         pd(8,23) = 0
         pd(10,21) = 0
         pd(16,5) = 0
         pd(20,14) = 0
         pd(1,16) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp
     +(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfm
     +ax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp
     +(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)
         pd(3,11) = 0
         pd(14,8) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(17,1) = 0
         pd(6,21) = 0
         pd(3,7) = 0
         pd(7,13) = loc_z*kdet*sp(1,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/s
     +t/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11
     +,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*kdet*sp(6,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(ex
     +p((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,
     +j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-s
     +p(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(
     +2,j)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6
     +,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfma
     +x-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(2
     +2,j))/st/Bfmax)*sp(7,j)
         pd(9,17) = -1/delt
         pd(12,10) = 0
         pd(3,15) = 1/fsulph*fdocs/delt
         pd(12,20) = 0
         pd(14,14) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z*katt*(1-1/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1))-loc_z*kdeac*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j)))/st/so4min)+1))-loc_z*km-1/delt
         pd(17,15) = loc_z*kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp(
     +(so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*
     +kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j
     +)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)
     +))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,
     +j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j
     +)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(16,j)/(exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(
     +-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp
     +(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2
     +*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+s
     +p(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp
     +(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(1
     +0,j)))/st/so4min)+loc_z*kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j
     +))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3
     ++sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(21,18) = 0
         pd(4,22) = 0
         pd(16,20) = -loc_z*kdet*sp(16,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)-loc_z*katt/(exp((Bfmax-sp(1,j)
     +-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/
     +Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)
     +-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(17,j)
         pd(7,19) = 0
         pd(14,1) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/
     +st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(1
     +1,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-s
     +p(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bf
     +max)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-s
     +p(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(6,3) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp
     +(1,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*s
     +p(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_
     +z/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ks
     +odoc+sp(3,j))**2*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-loc_z/(exp
     +((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)**2*(-kmax1/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kmax1*sp(3,j)/(ksodoc+sp(3,j))**2*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin*exp((doxmin-kmax1*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)
         pd(19,11) = 0
         pd(14,20) = loc_z*kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(14,j)
         pd(6,9) = 0
         pd(16,10) = -loc_z*kdeac*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp
     +(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(k
     +inno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)-loc_z*kreac*sp(17,j)/(exp((so4min-kmax3
     +*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox
     ++sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(11,4) = -loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+sp
     +(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kindo
     +x/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,
     +j)))/st/no3min)-loc_z*kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(ki
     +ndox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/
     +st/no3min)+1)**2*kmax2*kindox/(kindox+sp(4,j))**2*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))/st/no3min*exp((no3min-kmax2*kin
     +dox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(
     +3,j)))/st/no3min)+1/fcn*faa/foa*loc_z*(1-1/(exp(0.25D1*(0.4D0-ws)/
     +st)+1))+1/fcn*faa/foa/delt
         pd(15,12) = 0
         pd(18,12) = -0.1D-1*loc_z*km-1/delt
         pd(8,8) = -loc_z*(kl*vel**0.58D0+kdet/(exp((Bfmax-sp(1,j)-sp(8,
     +j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)
     ++1)+kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bfmax*exp((Bfm
     +ax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(
     +22,j))/st/Bfmax))-loc_z*katt/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-
     +sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st/Bf
     +max*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-
     +sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(13,15) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*
     +sp(13,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((amming-sp(5,j))/st
     +/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4m
     +in-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox
     +/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kdeac
     +*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(ks
     +sdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st
     +/so4min)+1)**2*(-kmax3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ks
     +so4+sp(15,j))**2*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(17,j)/(exp((so4mi
     +n-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/
     +(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax
     +3/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))**2*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)
     +/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +)/st/so4min)
         pd(5,17) = -0.1D-1*fcn*loc_z*km-fcn/delt
         pd(11,21) = -loc_z*km-1/delt
         pd(19,2) = 0
         pd(1,1) = loc_z*kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,
     +j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2/st
     +/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(
     +6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(7,j)-loc_z*kdeac*(1-
     +1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/doxmin)+1))
         pd(2,2) = -loc_z*kdeac*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+
     +sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))
         pd(15,6) = 0
         pd(21,13) = loc_z*kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(
     +13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**2
     +/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(
     +11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+loc_z*katt/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-
     +sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(21,j)
         pd(8,2) = 0
         pd(19,16) = 0
         pd(20,19) = 0
         pd(22,12) = 0
         pd(1,21) = 0
         pd(7,2) = -loc_z*km+loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Y
     +o*kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxm
     +in-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmi
     +n)+1)-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-1/delt
         pd(10,5) = 0
         pd(13,3) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*s
     +p(13,j)*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(kss
     +o4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinn
     +o3/(kinno3+sp(10,j)))/st/so4min)+1)-loc_z/(exp((amming-sp(5,j))/st
     +/amming)+1)*Ys*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))
     +/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3
     +,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+
     +1)-loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(13,j)*sp(
     +15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(
     +4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,
     +j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kdea
     +c*sp(14,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/s
     +t/so4min)+1)**2*(-kmax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))
     +*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(17,j)/(exp((so4
     +min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindo
     +x/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*(-km
     +ax3*sp(15,j)/(ksso4+sp(15,j))/(kssdoc+sp(3,j))*kindox/(kindox+sp(4
     +,j))*kinno3/(kinno3+sp(10,j))+kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))**2*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(
     +10,j)))/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3
     +,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,
     +j)))/st/so4min)
         pd(4,6) = 0
         pd(5,6) = fcn*loc_z*kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp
     +(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)**
     +2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp
     +(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+fcn*loc_z*katt/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)**2/st/Bfmax*exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)*sp(9,j)
         pd(22,18) = 0
         pd(23,21) = loc_z*kdeac*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
         pd(11,10) = -loc_z*kdeac*sp(9,j)/(exp((no3min-kmax2*kindox/(kin
     +dox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/s
     +t/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+sp(10,j))*s
     +p(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ks
     +no3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*exp((no3min-k
     +max2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j)))/st/no3min)-loc_z*kreac*sp(12,j)/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)**2*(-kmax2*kindox/(kindox+sp(4,j))/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+kmax2*kindox/(kindox+sp(4,j))*s
     +p(10,j)/(ksno3+sp(10,j))**2*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min*ex
     +p((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*
     +sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1/fcn*faa*foo/fdoco/foa/fnitr
     +a*fdocn/delt
         pd(17,10) = loc_z*kdeac*sp(13,j)/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/(ksso4+sp(1
     +5,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kin
     +no3+sp(10,j)))/st/so4min)+loc_z*kdeac*sp(14,j)/(exp((so4min-kmax3*
     +sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+
     +sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-kmax3*sp(15,j)/
     +(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*
     +kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(16,j)/(exp((so
     +4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kind
     +ox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)**2*kma
     +x3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kind
     +ox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((so4min-kmax
     +3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindo
     +x+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+loc_z*kreac*sp(17,
     +j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp
     +(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min
     +)+1)**2*kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*k
     +indox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j))**2/st/so4min*exp((
     +so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*ki
     +ndox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)
         pd(5,23) = -0.1D-1*fcn*loc_z*km-fcn/delt
         pd(7,7) = -0.1D-1*loc_z*km-loc_z*katt*(1-1/(exp((Bfmax-sp(1,j)-
     +sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/B
     +fmax)+1))-1/delt
         pd(13,9) = 0
         pd(13,21) = 0
         pd(21,4) = loc_z/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j
     +)*kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))/(exp((ammin-kmax4*s
     +p(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z/
     +(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2/(exp((ammin-kmax4*sp(5,j)/(k
     +samm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-loc_z/(exp((amm
     +ing-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j
     +))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+s
     +p(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox
     ++sp(4,j))**2)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/ammin)-loc_z*kdeac*sp(21,j)/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)*
     +*2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(ksox+sp(4,j))+kmax4*sp(5,j)/(k
     +samm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2)/st/ammin*exp((ammin-kmax4
     +*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)-loc_z*k
     +reac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(k
     +sox+sp(4,j)))/st/ammin)+1)**2*(-kmax4*sp(5,j)/(ksamm+sp(5,j))/(kso
     +x+sp(4,j))+kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))**2
     +)/st/ammin*exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/ammin)
      end
