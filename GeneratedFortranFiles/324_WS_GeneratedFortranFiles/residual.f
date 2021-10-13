c      
c     SUBROUTINE residual
c      
      subroutine residual(funcs,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension funcs(ncomp)
        call switches(j)
          funcs(1) = kreac*sp(6,j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp
     +(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)+0.1D-1*km*sp(6,j)+kl*
     +sp(6,j)*vel**0.58D0+kdet*sp(6,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13
     +,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-katt
     +*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,
     +j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,j)-kdeac*sp(1,j)*(1-1/(ex
     +p((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/s
     +t/doxmin)+1))+sp(6,j)/delt-spold(6,j)/delt
          funcs(2) = -kdeac*sp(2,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))+kreac*sp(7,j)/
     +(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))
     +)/st/doxmin)+1)+0.1D-1*km*sp(7,j)-kl*sp(6,j)*vel**0.58D0-kdet*sp(6
     +,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)
     +-sp(16,j)-sp(22,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(
     +8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfma
     +x)+1))*sp(7,j)+sp(7,j)/delt-spold(7,j)/delt
          funcs(3) = -fdoco/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*
     +sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodo
     +c+sp(3,j))*sp(4,j)/(ksox+sp(4,j))-fdoco/(exp((doxmin-kmax1*sp(3,j)
     +/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2
     +,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))+kpd*sp(18,j)-(
     +sp(3,j)-1/fnitra*fdocn*sp(10,j)-1/fsulph*fdocs*sp(15,j))/delt+(spo
     +ld(3,j)-1/fnitra*fdocn*spold(10,j)-1/fsulph*fdocs*spold(15,j))/del
     +t
          funcs(4) = foo/fdoco*fdocn/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1)*kmax2*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(1
     +0,j))*sp(3,j)/(ksndoc+sp(3,j))+foo/fdoco*fdocn/(exp((no3min-kmax2*
     +kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+
     +sp(3,j)))/st/no3min)+1)*kmax2*sp(9,j)*kindox/(kindox+sp(4,j))*sp(1
     +0,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j))+foo/fdoco/fnitra*fd
     +ocn*sp(10,j)/delt-foo/fdoco/fnitra*fdocn*spold(10,j)/delt
          funcs(5) = -fcn/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*s
     +p(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(
     +ksndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j
     +)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+fcn*(kl
     +*sp(8,j)*vel**0.58D0+kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(1
     +3,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))-fc
     +n*katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-
     +sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(9,j)-fcn*kdeac*sp(9,j
     +)*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+
     +sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))+fcn*kreac*sp(12
     +,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(
     +10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-fcn/(exp((amming-sp
     +(5,j))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(
     +4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp
     +(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-fcn*km*sp(9,j)-fcn/(exp((ammin
     +g-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))
     +*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j)
     +)*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)-fcn*km*sp(13,j)-fcn*km*sp(
     +14,j)-0.1D-1*fcn*km*sp(16,j)-0.1D-1*fcn*km*sp(17,j)-fcn*km*sp(20,j
     +)-fcn*km*sp(21,j)-0.1D-1*fcn*km*sp(22,j)-0.1D-1*fcn*km*sp(23,j)+(f
     +cn*fdoco*foa+faa*foo)/fdoco/foa*kpd*sp(18,j)-faa/foa*(250-sp(4,j))
     +*(1-1/(exp(0.25D1*(0.4D0-ws)/st)+1))-(faa*foo/fdoco/foa*sp(3,j)-fa
     +a/foa*sp(4,j)+sp(5,j)+fcn*sp(9,j)-faa*foo/fdoco/foa/fnitra*fdocn*s
     +p(10,j)+fcn*sp(13,j)+fcn*sp(14,j)-faa/fsulph*fdocs*foo/fdoco/foa*s
     +p(15,j)+fcn*sp(16,j)+fcn*sp(17,j)+fcn*sp(20,j)+fcn*sp(21,j)+fcn*sp
     +(22,j)+fcn*sp(23,j))/delt+(faa*foo/fdoco/foa*spold(3,j)-faa/foa*sp
     +old(4,j)+spold(5,j)+fcn*spold(9,j)-faa*foo/fdoco/foa/fnitra*fdocn*
     +spold(10,j)+fcn*spold(13,j)+fcn*spold(14,j)-faa/fsulph*fdocs*foo/f
     +doco/foa*spold(15,j)+fcn*spold(16,j)+fcn*spold(17,j)+fcn*spold(20,
     +j)+fcn*spold(21,j)+fcn*spold(22,j)+fcn*spold(23,j))/delt
          funcs(6) = -km*sp(1,j)-0.1D-1*km*sp(6,j)+1/(exp((amming-sp(5,j
     +))/st/amming)+1)*Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/doxmin)+1)-kl*sp(1,j)*vel**0.58D0-kdet*sp(1,j
     +)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)+1)-kl*sp(6,j)*vel**0.58D0-kdet*sp(6,j)
     +/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1))*sp(2,j)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,j)-(sp(1,
     +j)+sp(6,j))/delt+(spold(1,j)+spold(6,j))/delt
          funcs(7) = -km*sp(2,j)-0.1D-1*km*sp(7,j)+1/(exp((amming-sp(5,j
     +))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/doxmin)+1)+kl*sp(1,j)*vel**0.58D0+kdet*sp(1,j
     +)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-s
     +p(16,j)-sp(22,j))/st/Bfmax)+1)+kl*sp(6,j)*vel**0.58D0+kdet*sp(6,j)
     +/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp
     +(16,j)-sp(22,j))/st/Bfmax)+1)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j
     +)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+
     +1))*sp(2,j)-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j
     +)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,j)-(sp(2,
     +j)+sp(7,j))/delt+(spold(2,j)+spold(7,j))/delt
          funcs(8) = -1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(
     +9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-kl*sp(8,j
     +)*vel**0.58D0-kdet*sp(8,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp
     +(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+katt*(1-1/
     +(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(
     +16,j)-sp(22,j))/st/Bfmax)+1))*sp(9,j)+kdeac*sp(9,j)*(1-1/(exp((no3
     +min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j
     +)/(ksndoc+sp(3,j)))/st/no3min)+1))-kreac*sp(12,j)/(exp((no3min-kma
     +x2*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksnd
     +oc+sp(3,j)))/st/no3min)+1)+km*sp(9,j)+sp(9,j)/delt-spold(9,j)/delt
          funcs(9) = -kdeac*sp(8,j)*(1-1/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1))-kdeac*sp(9,j)*(1-1/(exp((no3min-kmax2*kindox/(kindox+
     +sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no
     +3min)+1))+kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j)
     +)*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1
     +)+kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,
     +j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)+0.1D-1
     +*km*sp(11,j)-km*sp(20,j)-km*sp(21,j)-0.1D-1*km*sp(22,j)-0.1D-1*km*
     +sp(23,j)+1/fcn*faa*foo/fdoco/foa*kpd*sp(18,j)+kmax5*sp(5,j)/(ksamm
     ++sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)-1/fcn*fa
     +a/foa*(250-sp(4,j))*(1-1/(exp(0.25D1*(0.4D0-ws)/st)+1))-(sp(1,j)+s
     +p(2,j)+1/fcn*faa*foo/fdoco/foa*sp(3,j)-1/fcn*faa/foa*sp(4,j)+1/fcn
     +*sp(5,j)+sp(6,j)+sp(7,j)+sp(8,j)+sp(9,j)-1/fcn*faa*foo/fdoco/foa/f
     +nitra*fdocn*sp(10,j)+sp(13,j)+sp(14,j)-1/fcn*faa/fsulph*fdocs*foo/
     +fdoco/foa*sp(15,j)+sp(16,j)+sp(17,j)+sp(18,j)+sp(20,j)+sp(21,j)+sp
     +(22,j)+sp(23,j))/delt+(spold(1,j)+spold(2,j)+1/fcn*faa*foo/fdoco/f
     +oa*spold(3,j)-1/fcn*faa/foa*spold(4,j)+1/fcn*spold(5,j)+spold(6,j)
     ++spold(7,j)+spold(8,j)+spold(9,j)-1/fcn*faa*foo/fdoco/foa/fnitra*f
     +docn*spold(10,j)+spold(13,j)+spold(14,j)-1/fcn*faa/fsulph*fdocs*fo
     +o/fdoco/foa*spold(15,j)+spold(16,j)+spold(17,j)+spold(18,j)+spold(
     +20,j)+spold(21,j)+spold(22,j)+spold(23,j))/delt
          funcs(10) = fnitra/fdocn*fdocs/(exp((so4min-kmax3*sp(15,j)/(ks
     +so4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kin
     +no3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso
     +4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno
     +3/(kinno3+sp(10,j))+fnitra/fdocn*fdocs/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(
     +ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*k
     +inno3/(kinno3+sp(10,j))+fnitra/fdocn/fsulph*fdocs*sp(15,j)/delt-fn
     +itra/fdocn/fsulph*fdocs*spold(15,j)/delt
          funcs(11) = -kl*sp(11,j)*vel**0.58D0-kdet*sp(11,j)/(exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(12,j)
     +-kdeac*sp(9,j)*(1-1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(
     +10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1))+kr
     +eac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/(
     +ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)-km*sp(20,j
     +)-km*sp(21,j)-0.1D-1*km*sp(22,j)-0.1D-1*km*sp(23,j)+1/fcn*faa*foo/
     +fdoco/foa*kpd*sp(18,j)+kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0
     +.1D0*amming-sp(5,j))/st/amming)+1)-1/fcn*faa/foa*(250-sp(4,j))*(1-
     +1/(exp(0.25D1*(0.4D0-ws)/st)+1))-(sp(1,j)+sp(2,j)+1/fcn*faa*foo/fd
     +oco/foa*sp(3,j)-1/fcn*faa/foa*sp(4,j)+1/fcn*sp(5,j)+sp(6,j)+sp(7,j
     +)+sp(8,j)+sp(9,j)-1/fcn*faa*foo/fdoco/foa/fnitra*fdocn*sp(10,j)+sp
     +(11,j)+sp(13,j)+sp(14,j)-1/fcn*faa/fsulph*fdocs*foo/fdoco/foa*sp(1
     +5,j)+sp(16,j)+sp(17,j)+sp(18,j)+sp(20,j)+sp(21,j)+sp(22,j)+sp(23,j
     +))/delt+(spold(1,j)+spold(2,j)+1/fcn*faa*foo/fdoco/foa*spold(3,j)-
     +1/fcn*faa/foa*spold(4,j)+1/fcn*spold(5,j)+spold(6,j)+spold(7,j)+sp
     +old(8,j)+spold(9,j)-1/fcn*faa*foo/fdoco/foa/fnitra*fdocn*spold(10,
     +j)+spold(11,j)+spold(13,j)+spold(14,j)-1/fcn*faa/fsulph*fdocs*foo/
     +fdoco/foa*spold(15,j)+spold(16,j)+spold(17,j)+spold(18,j)+spold(20
     +,j)+spold(21,j)+spold(22,j)+spold(23,j))/delt
          funcs(12) = -km*sp(8,j)-1/(exp((amming-sp(5,j))/st/amming)+1)*
     +Yo*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(
     +exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))
     +/st/doxmin)+1)-km*sp(9,j)-0.1D-1*km*sp(11,j)-1/(exp((amming-sp(5,j
     +))/st/amming)+1)*Yo*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)
     +/(ksox+sp(4,j))/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j
     +)/(ksox+sp(4,j)))/st/doxmin)+1)-km*sp(13,j)-km*sp(14,j)-0.1D-1*km*
     +sp(16,j)-0.1D-1*km*sp(17,j)+kpd*sp(18,j)-kmax5*sp(5,j)/(ksamm+sp(5
     +,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))/st/amming)+1)-(-sp(1,j)-sp(
     +2,j)-sp(6,j)-sp(7,j)-sp(18,j))/delt+(-spold(1,j)-spold(2,j)-spold(
     +6,j)-spold(7,j)-spold(18,j))/delt
          funcs(13) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(
     +13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-kl*sp(13,j)*vel**0.58D0-
     +kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1))*sp(14,j)+kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))-kreac*sp(17,j)/(e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-
     +km*sp(13,j)-0.1D-1*km*sp(16,j)-0.1D-1*km*sp(17,j)-(sp(13,j)+sp(16,
     +j)+sp(17,j))/delt+(spold(13,j)+spold(16,j)+spold(17,j))/delt
          funcs(14) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(
     +14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)+kl*sp(13,j)*vel**0.58D0+
     +kdet*sp(13,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j
     +)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)-katt*(1-1/(exp((Bfmax-s
     +p(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j
     +))/st/Bfmax)+1))*sp(14,j)-kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*s
     +p(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+s
     +p(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))+kreac*sp(17,j)/(e
     +xp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j)
     +)*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)-
     +km*sp(14,j)-sp(14,j)/delt+spold(14,j)/delt
          funcs(15) = -fsulph/fdocs/foo*fdoco*foa*sp(20,j)*kmax4*sp(5,j)
     +/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/
     +(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-fsulph/fdocs/
     +foo*fdoco*foa*sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox
     ++sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+
     +sp(4,j)))/st/ammin)+1)-fsulph/fdocs*kpd*sp(18,j)+fsulph/fdocs/foo*
     +fdoco*(250-sp(4,j))*(1-1/(exp(0.25D1*(0.4D0-ws)/st)+1))-(-fsulph/f
     +docs*sp(3,j)+fsulph/fdocs/foo*fdoco*sp(4,j)+fsulph/fnitra*fdocn/fd
     +ocs*sp(10,j)+sp(15,j))/delt+(-fsulph/fdocs*spold(3,j)+fsulph/fdocs
     +/foo*fdoco*spold(4,j)+fsulph/fnitra*fdocn/fdocs*spold(10,j)+spold(
     +15,j))/delt
          funcs(16) = -kl*sp(16,j)*vel**0.58D0-kdet*sp(16,j)/(exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(17,j)
     +-kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*
     +sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp
     +(10,j)))/st/so4min)+1))+kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)
     +/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))
     +*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)+0.1D-1*km*sp(17,j)+sp(17,
     +j)/delt-spold(17,j)/delt
          funcs(17) = kdeac*sp(13,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1))+kdeac*sp(14,j)*(1-1/(exp((s
     +o4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kin
     +dox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/st/so4min)+1))-krea
     +c*sp(16,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(k
     +ssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j)))/s
     +t/so4min)+1)-kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(
     +15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(ki
     +nno3+sp(10,j)))/st/so4min)+1)-0.1D-1*km*sp(16,j)-0.1D-1*km*sp(17,j
     +)-(sp(16,j)+sp(17,j))/delt+(spold(16,j)+spold(17,j))/delt
          funcs(18) = -0.1D-1*km*sp(12,j)-km*sp(20,j)-km*sp(21,j)-0.1D-1
     +*km*sp(22,j)-0.1D-1*km*sp(23,j)+1/fcn*faa*foo/fdoco/foa*kpd*sp(18,
     +j)+kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1D2*(0.1D0*amming-sp(5,j))
     +/st/amming)+1)-1/fcn*faa/foa*(250-sp(4,j))*(1-1/(exp(0.25D1*(0.4D0
     +-ws)/st)+1))-(sp(1,j)+sp(2,j)+1/fcn*faa*foo/fdoco/foa*sp(3,j)-1/fc
     +n*faa/foa*sp(4,j)+1/fcn*sp(5,j)+sp(6,j)+sp(7,j)+sp(8,j)+sp(9,j)-1/
     +fcn*faa*foo/fdoco/foa/fnitra*fdocn*sp(10,j)+sp(11,j)+sp(12,j)+sp(1
     +3,j)+sp(14,j)-1/fcn*faa/fsulph*fdocs*foo/fdoco/foa*sp(15,j)+sp(16,
     +j)+sp(17,j)+sp(18,j)+sp(20,j)+sp(21,j)+sp(22,j)+sp(23,j))/delt+(sp
     +old(1,j)+spold(2,j)+1/fcn*faa*foo/fdoco/foa*spold(3,j)-1/fcn*faa/f
     +oa*spold(4,j)+1/fcn*spold(5,j)+spold(6,j)+spold(7,j)+spold(8,j)+sp
     +old(9,j)-1/fcn*faa*foo/fdoco/foa/fnitra*fdocn*spold(10,j)+spold(11
     +,j)+spold(12,j)+spold(13,j)+spold(14,j)-1/fcn*faa/fsulph*fdocs*foo
     +/fdoco/foa*spold(15,j)+spold(16,j)+spold(17,j)+spold(18,j)+spold(2
     +0,j)+spold(21,j)+spold(22,j)+spold(23,j))/delt
          funcs(19) = sp(19,j)-spold(19,j)
          funcs(20) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-
     +kl*sp(20,j)*vel**0.58D0-kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)
     ++katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)+kdeac*sp(21,j)*(
     +1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/ammin)+1))-kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+
     +sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-km*sp(20,j)-0.1D-1*k
     +m*sp(22,j)-0.1D-1*km*sp(23,j)-(sp(20,j)+sp(22,j)+sp(23,j))/delt+(s
     +pold(20,j)+spold(22,j)+spold(23,j))/delt
          funcs(21) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+
     +kl*sp(20,j)*vel**0.58D0+kdet*sp(20,j)/(exp((Bfmax-sp(1,j)-sp(8,j)-
     +sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1)
     +-katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-s
     +p(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)-kdeac*sp(21,j)*(
     +1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/ammin)+1))+kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+
     +sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-km*sp(21,j)-sp(21,j)
     +/delt+spold(21,j)/delt
          funcs(22) = -kl*sp(22,j)*vel**0.58D0-kdet*sp(22,j)/(exp((Bfmax
     +-sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22
     +,j))/st/Bfmax)+1)+katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-s
     +p(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(23,j)
     +-kdeac*sp(21,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(
     +4,j)/(ksox+sp(4,j)))/st/ammin)+1))+kreac*sp(23,j)/(exp((ammin-kmax
     +4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)+0.1
     +D-1*km*sp(23,j)+sp(23,j)/delt-spold(23,j)/delt
          funcs(23) = kdeac*sp(20,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))+kdeac*sp(21,j)*(
     +1-1/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j
     +)))/st/ammin)+1))-kreac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+
     +sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)-kreac*sp(23,j)/(exp(
     +(ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/am
     +min)+1)-0.1D-1*km*sp(22,j)-0.1D-1*km*sp(23,j)-(sp(22,j)+sp(23,j))/
     +delt+(spold(22,j)+spold(23,j))/delt
      end
