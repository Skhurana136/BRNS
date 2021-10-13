c      
c     SUBROUTINE rates
c      
      subroutine rates(j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            r(1,j) = 1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(1,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))
            r(2,j) = 1/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,
     +j)/(ksox+sp(4,j)))/st/doxmin)+1)*kmax1*sp(2,j)*sp(3,j)/(ksodoc+sp(
     +3,j))*sp(4,j)/(ksox+sp(4,j))
            r(3,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(1
     +,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-k
     +max1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1
     +)
            r(4,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yo*kmax1*sp(2
     +,j)*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j))/(exp((doxmin-k
     +max1*sp(3,j)/(ksodoc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1
     +)
            r(5,j) = kl*sp(1,j)*vel**0.58E0+kdet*sp(1,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)
            r(6,j) = kl*sp(6,j)*vel**0.58E0+kdet*sp(6,j)/(exp((Bfmax-sp(
     +1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))
     +/st/Bfmax)+1)
            r(7,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(2,j)
            r(8,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(2
     +0,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(7,j)
            r(9,j) = kdeac*sp(1,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(ksod
     +oc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))
            r(10,j) = kdeac*sp(2,j)*(1-1/(exp((doxmin-kmax1*sp(3,j)/(kso
     +doc+sp(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1))
            r(11,j) = kreac*sp(6,j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
            r(12,j) = kreac*sp(7,j)/(exp((doxmin-kmax1*sp(3,j)/(ksodoc+s
     +p(3,j))*sp(4,j)/(ksox+sp(4,j)))/st/doxmin)+1)
            r(13,j) = km*sp(1,j)
            r(14,j) = km*sp(2,j)
            r(15,j) = 0.1E-1*km*sp(6,j)
            r(16,j) = 0.1E-1*km*sp(7,j)
            r(17,j) = 1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2
     +*sp(8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j))
            r(18,j) = 1/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10
     +,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)*kmax2
     +*sp(9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)
     +/(ksndoc+sp(3,j))
            r(19,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(
     +8,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
            r(20,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Yn*kmax2*sp(
     +9,j)*kindox/(kindox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ks
     +ndoc+sp(3,j))/(exp((no3min-kmax2*kindox/(kindox+sp(4,j))*sp(10,j)/
     +(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3min)+1)
            r(21,j) = kl*sp(8,j)*vel**0.58E0+kdet*sp(8,j)/(exp((Bfmax-sp
     +(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j)
     +)/st/Bfmax)+1)
            r(22,j) = kl*sp(11,j)*vel**0.58E0+kdet*sp(11,j)/(exp((Bfmax-
     +sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,
     +j))/st/Bfmax)+1)
            r(23,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(9,j)
            r(24,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(12,j)
            r(25,j) = kdeac*sp(8,j)*(1-1/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1))
            r(26,j) = kdeac*sp(9,j)*(1-1/(exp((no3min-kmax2*kindox/(kind
     +ox+sp(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st
     +/no3min)+1))
            r(27,j) = kreac*sp(11,j)/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1)
            r(28,j) = kreac*sp(12,j)/(exp((no3min-kmax2*kindox/(kindox+s
     +p(4,j))*sp(10,j)/(ksno3+sp(10,j))*sp(3,j)/(ksndoc+sp(3,j)))/st/no3
     +min)+1)
            r(29,j) = km*sp(8,j)
            r(30,j) = km*sp(9,j)
            r(31,j) = 0.1E-1*km*sp(11,j)
            r(32,j) = 0.1E-1*km*sp(12,j)
            r(33,j) = 1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(
     +3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j)))/st/so4min)+1)*kmax3*sp(13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))
            r(34,j) = 1/(exp((so4min-kmax3*sp(15,j)/(ksso4+sp(15,j))*sp(
     +3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10
     +,j)))/st/so4min)+1)*kmax3*sp(14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,
     +j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/(kinno3+sp(10,j
     +))
            r(35,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(
     +13,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)
            r(36,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ys*kmax3*sp(
     +14,j)*sp(15,j)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(k
     +indox+sp(4,j))*kinno3/(kinno3+sp(10,j))/(exp((so4min-kmax3*sp(15,j
     +)/(ksso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j)
     +)*kinno3/(kinno3+sp(10,j)))/st/so4min)+1)
            r(37,j) = kl*sp(13,j)*vel**0.58E0+kdet*sp(13,j)/(exp((Bfmax-
     +sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,
     +j))/st/Bfmax)+1)
            r(38,j) = kl*sp(16,j)*vel**0.58E0+kdet*sp(16,j)/(exp((Bfmax-
     +sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,
     +j))/st/Bfmax)+1)
            r(39,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(14,j)
            r(40,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(17,j)
            r(41,j) = kdeac*sp(13,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1))
            r(42,j) = kdeac*sp(14,j)*(1-1/(exp((so4min-kmax3*sp(15,j)/(k
     +sso4+sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*ki
     +nno3/(kinno3+sp(10,j)))/st/so4min)+1))
            r(43,j) = kreac*sp(16,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)
            r(44,j) = kreac*sp(17,j)/(exp((so4min-kmax3*sp(15,j)/(ksso4+
     +sp(15,j))*sp(3,j)/(kssdoc+sp(3,j))*kindox/(kindox+sp(4,j))*kinno3/
     +(kinno3+sp(10,j)))/st/so4min)+1)
            r(45,j) = km*sp(13,j)
            r(46,j) = km*sp(14,j)
            r(47,j) = 0.1E-1*km*sp(16,j)
            r(48,j) = 0.1E-1*km*sp(17,j)
            r(49,j) = sp(20,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +ox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(kso
     +x+sp(4,j)))/st/ammin)+1)
            r(50,j) = sp(21,j)*kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ks
     +ox+sp(4,j))/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(kso
     +x+sp(4,j)))/st/ammin)+1)
            r(51,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(20,j)*
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
            r(52,j) = 1/(exp((amming-sp(5,j))/st/amming)+1)*Ya*sp(21,j)*
     +kmax4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j))/(exp((ammin-k
     +max4*sp(5,j)/(ksamm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
            r(53,j) = kl*sp(20,j)*vel**0.58E0+kdet*sp(20,j)/(exp((Bfmax-
     +sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,
     +j))/st/Bfmax)+1)
            r(54,j) = kl*sp(22,j)*vel**0.58E0+kdet*sp(22,j)/(exp((Bfmax-
     +sp(1,j)-sp(8,j)-sp(13,j)-sp(20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,
     +j))/st/Bfmax)+1)
            r(55,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(21,j)
            r(56,j) = katt*(1-1/(exp((Bfmax-sp(1,j)-sp(8,j)-sp(13,j)-sp(
     +20,j)-sp(6,j)-sp(11,j)-sp(16,j)-sp(22,j))/st/Bfmax)+1))*sp(23,j)
            r(57,j) = kdeac*sp(20,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
            r(58,j) = kdeac*sp(21,j)*(1-1/(exp((ammin-kmax4*sp(5,j)/(ksa
     +mm+sp(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1))
            r(59,j) = kreac*sp(22,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
            r(60,j) = kreac*sp(23,j)/(exp((ammin-kmax4*sp(5,j)/(ksamm+sp
     +(5,j))*sp(4,j)/(ksox+sp(4,j)))/st/ammin)+1)
            r(61,j) = km*sp(20,j)
            r(62,j) = km*sp(21,j)
            r(63,j) = 0.1E-1*km*sp(22,j)
            r(64,j) = 0.1E-1*km*sp(23,j)
            r(65,j) = kpd*sp(18,j)
            r(66,j) = kmax5*sp(5,j)/(ksamm+sp(5,j))/(exp(0.1E2*(0.1E0*am
     +ming-sp(5,j))/st/amming)+1)
            r(67,j) = (250-sp(4,j))*(1-1/(exp(0.25E1*(0.4E0-ws)/st)+1))
      end
