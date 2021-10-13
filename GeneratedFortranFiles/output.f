c      
c     SUBROUTINE out
c      
      subroutine out(j,nt,time,depth,v_out,v_int)
        include 'common_geo.inc'
        include 'common.inc'
        real*8 time
        if (nt.eq.1.and.j.eq.1) then
          open(unit=11,file='Bfo1.dat',status='replace')
          close(11)
          open(unit=11,file='Bmo1.dat',status='replace')
          close(11)
          open(unit=11,file='doc1.dat',status='replace')
          close(11)
          open(unit=11,file='dox1.dat',status='replace')
          close(11)
          open(unit=11,file='Amm1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifo1.dat',status='replace')
          close(11)
          open(unit=11,file='Bimo1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfn1.dat',status='replace')
          close(11)
          open(unit=11,file='Bmn1.dat',status='replace')
          close(11)
          open(unit=11,file='nitra1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifn1.dat',status='replace')
          close(11)
          open(unit=11,file='Bimn1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfs1.dat',status='replace')
          close(11)
          open(unit=11,file='Bms1.dat',status='replace')
          close(11)
          open(unit=11,file='sulpha1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifs1.dat',status='replace')
          close(11)
          open(unit=11,file='Bims1.dat',status='replace')
          close(11)
          open(unit=11,file='POM1.dat',status='replace')
          close(11)
          open(unit=11,file='tr1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfa1.dat',status='replace')
          close(11)
          open(unit=11,file='Bma1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifa1.dat',status='replace')
          close(11)
          open(unit=11,file='Bima1.dat',status='replace')
          close(11)
          open(unit=11,file='xrate1.dat',status='replace')
          close(11)
          open(unit=11,file='xrate2.dat',status='replace')
          close(11)
          open(unit=11,file='xrate3.dat',status='replace')
          close(11)
          open(unit=11,file='xrate4.dat',status='replace')
          close(11)
          open(unit=11,file='xrate5.dat',status='replace')
          close(11)
          open(unit=11,file='xrate6.dat',status='replace')
          close(11)
          open(unit=11,file='xrate7.dat',status='replace')
          close(11)
          open(unit=11,file='xrate8.dat',status='replace')
          close(11)
          open(unit=11,file='xrate9.dat',status='replace')
          close(11)
          open(unit=11,file='xrate10.dat',status='replace')
          close(11)
          open(unit=11,file='xrate11.dat',status='replace')
          close(11)
          open(unit=11,file='xrate12.dat',status='replace')
          close(11)
          open(unit=11,file='xrate13.dat',status='replace')
          close(11)
          open(unit=11,file='xrate14.dat',status='replace')
          close(11)
          open(unit=11,file='xrate15.dat',status='replace')
          close(11)
          open(unit=11,file='xrate16.dat',status='replace')
          close(11)
          open(unit=11,file='xrate17.dat',status='replace')
          close(11)
          open(unit=11,file='xrate18.dat',status='replace')
          close(11)
          open(unit=11,file='xrate19.dat',status='replace')
          close(11)
          open(unit=11,file='xrate20.dat',status='replace')
          close(11)
          open(unit=11,file='xrate21.dat',status='replace')
          close(11)
          open(unit=11,file='xrate22.dat',status='replace')
          close(11)
          open(unit=11,file='xrate23.dat',status='replace')
          close(11)
          open(unit=11,file='xrate24.dat',status='replace')
          close(11)
          open(unit=11,file='xrate25.dat',status='replace')
          close(11)
          open(unit=11,file='xrate26.dat',status='replace')
          close(11)
          open(unit=11,file='xrate27.dat',status='replace')
          close(11)
          open(unit=11,file='xrate28.dat',status='replace')
          close(11)
          open(unit=11,file='xrate29.dat',status='replace')
          close(11)
          open(unit=11,file='xrate30.dat',status='replace')
          close(11)
          open(unit=11,file='xrate31.dat',status='replace')
          close(11)
          open(unit=11,file='xrate32.dat',status='replace')
          close(11)
          open(unit=11,file='xrate33.dat',status='replace')
          close(11)
          open(unit=11,file='xrate34.dat',status='replace')
          close(11)
          open(unit=11,file='xrate35.dat',status='replace')
          close(11)
          open(unit=11,file='xrate36.dat',status='replace')
          close(11)
          open(unit=11,file='xrate37.dat',status='replace')
          close(11)
          open(unit=11,file='xrate38.dat',status='replace')
          close(11)
          open(unit=11,file='xrate39.dat',status='replace')
          close(11)
          open(unit=11,file='xrate40.dat',status='replace')
          close(11)
          open(unit=11,file='xrate41.dat',status='replace')
          close(11)
          open(unit=11,file='xrate42.dat',status='replace')
          close(11)
          open(unit=11,file='xrate43.dat',status='replace')
          close(11)
          open(unit=11,file='xrate44.dat',status='replace')
          close(11)
          open(unit=11,file='xrate45.dat',status='replace')
          close(11)
          open(unit=11,file='xrate46.dat',status='replace')
          close(11)
          open(unit=11,file='xrate47.dat',status='replace')
          close(11)
          open(unit=11,file='xrate48.dat',status='replace')
          close(11)
          open(unit=11,file='xrate49.dat',status='replace')
          close(11)
          open(unit=11,file='xrate50.dat',status='replace')
          close(11)
          open(unit=11,file='xrate51.dat',status='replace')
          close(11)
          open(unit=11,file='xrate52.dat',status='replace')
          close(11)
          open(unit=11,file='xrate53.dat',status='replace')
          close(11)
          open(unit=11,file='xrate54.dat',status='replace')
          close(11)
          open(unit=11,file='xrate55.dat',status='replace')
          close(11)
          open(unit=11,file='xrate56.dat',status='replace')
          close(11)
          open(unit=11,file='xrate57.dat',status='replace')
          close(11)
          open(unit=11,file='xrate58.dat',status='replace')
          close(11)
          open(unit=11,file='xrate59.dat',status='replace')
          close(11)
          open(unit=11,file='xrate60.dat',status='replace')
          close(11)
          open(unit=11,file='xrate61.dat',status='replace')
          close(11)
          open(unit=11,file='xrate62.dat',status='replace')
          close(11)
          open(unit=11,file='xrate63.dat',status='replace')
          close(11)
          open(unit=11,file='xrate64.dat',status='replace')
          close(11)
          open(unit=11,file='xrate65.dat',status='replace')
          close(11)
          open(unit=11,file='xrate66.dat',status='replace')
          close(11)
          open(unit=11,file='Bfo1.inp',status='replace')
          close(11)
          open(unit=11,file='Bmo1.inp',status='replace')
          close(11)
          open(unit=11,file='doc1.inp',status='replace')
          close(11)
          open(unit=11,file='dox1.inp',status='replace')
          close(11)
          open(unit=11,file='Amm1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifo1.inp',status='replace')
          close(11)
          open(unit=11,file='Bimo1.inp',status='replace')
          close(11)
          open(unit=11,file='Bfn1.inp',status='replace')
          close(11)
          open(unit=11,file='Bmn1.inp',status='replace')
          close(11)
          open(unit=11,file='nitra1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifn1.inp',status='replace')
          close(11)
          open(unit=11,file='Bimn1.inp',status='replace')
          close(11)
          open(unit=11,file='Bfs1.inp',status='replace')
          close(11)
          open(unit=11,file='Bms1.inp',status='replace')
          close(11)
          open(unit=11,file='sulpha1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifs1.inp',status='replace')
          close(11)
          open(unit=11,file='Bims1.inp',status='replace')
          close(11)
          open(unit=11,file='POM1.inp',status='replace')
          close(11)
          open(unit=11,file='tr1.inp',status='replace')
          close(11)
          open(unit=11,file='Bfa1.inp',status='replace')
          close(11)
          open(unit=11,file='Bma1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifa1.inp',status='replace')
          close(11)
          open(unit=11,file='Bima1.inp',status='replace')
          close(11)
        v_out = 0.1D-2
        v_int = 100
        endif
        if (time.le.v_out.and.v_out.lt.time+delt) then
          open(unit=11,file='Bfo1.dat',
     +      status='old',access='append')
          write(11,2000) sp(1,j),depth
 2000     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmo1.dat',
     +      status='old',access='append')
          write(11,2001) sp(2,j),depth
 2001     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='doc1.dat',
     +      status='old',access='append')
          write(11,2002) sp(3,j),depth
 2002     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='dox1.dat',
     +      status='old',access='append')
          write(11,2003) sp(4,j),depth
 2003     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Amm1.dat',
     +      status='old',access='append')
          write(11,2004) sp(5,j),depth
 2004     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifo1.dat',
     +      status='old',access='append')
          write(11,2005) sp(6,j),depth
 2005     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimo1.dat',
     +      status='old',access='append')
          write(11,2006) sp(7,j),depth
 2006     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfn1.dat',
     +      status='old',access='append')
          write(11,2007) sp(8,j),depth
 2007     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmn1.dat',
     +      status='old',access='append')
          write(11,2008) sp(9,j),depth
 2008     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='nitra1.dat',
     +      status='old',access='append')
          write(11,2009) sp(10,j),depth
 2009     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifn1.dat',
     +      status='old',access='append')
          write(11,2010) sp(11,j),depth
 2010     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimn1.dat',
     +      status='old',access='append')
          write(11,2011) sp(12,j),depth
 2011     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfs1.dat',
     +      status='old',access='append')
          write(11,2012) sp(13,j),depth
 2012     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bms1.dat',
     +      status='old',access='append')
          write(11,2013) sp(14,j),depth
 2013     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='sulpha1.dat',
     +      status='old',access='append')
          write(11,2014) sp(15,j),depth
 2014     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifs1.dat',
     +      status='old',access='append')
          write(11,2015) sp(16,j),depth
 2015     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bims1.dat',
     +      status='old',access='append')
          write(11,2016) sp(17,j),depth
 2016     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='POM1.dat',
     +      status='old',access='append')
          write(11,2017) sp(18,j),depth
 2017     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='tr1.dat',
     +      status='old',access='append')
          write(11,2018) sp(19,j),depth
 2018     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfa1.dat',
     +      status='old',access='append')
          write(11,2019) sp(20,j),depth
 2019     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bma1.dat',
     +      status='old',access='append')
          write(11,2020) sp(21,j),depth
 2020     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifa1.dat',
     +      status='old',access='append')
          write(11,2021) sp(22,j),depth
 2021     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bima1.dat',
     +      status='old',access='append')
          write(11,2022) sp(23,j),depth
 2022     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate1.dat',
     +      status='old',access='append')
          write(11,2023) r(1,j),depth
 2023     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate2.dat',
     +      status='old',access='append')
          write(11,2024) r(2,j),depth
 2024     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate3.dat',
     +      status='old',access='append')
          write(11,2025) r(3,j),depth
 2025     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate4.dat',
     +      status='old',access='append')
          write(11,2026) r(4,j),depth
 2026     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate5.dat',
     +      status='old',access='append')
          write(11,2027) r(5,j),depth
 2027     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate6.dat',
     +      status='old',access='append')
          write(11,2028) r(6,j),depth
 2028     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate7.dat',
     +      status='old',access='append')
          write(11,2029) r(7,j),depth
 2029     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate8.dat',
     +      status='old',access='append')
          write(11,2030) r(8,j),depth
 2030     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate9.dat',
     +      status='old',access='append')
          write(11,2031) r(9,j),depth
 2031     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate10.dat',
     +      status='old',access='append')
          write(11,2032) r(10,j),depth
 2032     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate11.dat',
     +      status='old',access='append')
          write(11,2033) r(11,j),depth
 2033     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate12.dat',
     +      status='old',access='append')
          write(11,2034) r(12,j),depth
 2034     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate13.dat',
     +      status='old',access='append')
          write(11,2035) r(13,j),depth
 2035     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate14.dat',
     +      status='old',access='append')
          write(11,2036) r(14,j),depth
 2036     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate15.dat',
     +      status='old',access='append')
          write(11,2037) r(15,j),depth
 2037     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate16.dat',
     +      status='old',access='append')
          write(11,2038) r(16,j),depth
 2038     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate17.dat',
     +      status='old',access='append')
          write(11,2039) r(17,j),depth
 2039     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate18.dat',
     +      status='old',access='append')
          write(11,2040) r(18,j),depth
 2040     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate19.dat',
     +      status='old',access='append')
          write(11,2041) r(19,j),depth
 2041     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate20.dat',
     +      status='old',access='append')
          write(11,2042) r(20,j),depth
 2042     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate21.dat',
     +      status='old',access='append')
          write(11,2043) r(21,j),depth
 2043     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate22.dat',
     +      status='old',access='append')
          write(11,2044) r(22,j),depth
 2044     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate23.dat',
     +      status='old',access='append')
          write(11,2045) r(23,j),depth
 2045     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate24.dat',
     +      status='old',access='append')
          write(11,2046) r(24,j),depth
 2046     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate25.dat',
     +      status='old',access='append')
          write(11,2047) r(25,j),depth
 2047     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate26.dat',
     +      status='old',access='append')
          write(11,2048) r(26,j),depth
 2048     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate27.dat',
     +      status='old',access='append')
          write(11,2049) r(27,j),depth
 2049     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate28.dat',
     +      status='old',access='append')
          write(11,2050) r(28,j),depth
 2050     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate29.dat',
     +      status='old',access='append')
          write(11,2051) r(29,j),depth
 2051     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate30.dat',
     +      status='old',access='append')
          write(11,2052) r(30,j),depth
 2052     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate31.dat',
     +      status='old',access='append')
          write(11,2053) r(31,j),depth
 2053     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate32.dat',
     +      status='old',access='append')
          write(11,2054) r(32,j),depth
 2054     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate33.dat',
     +      status='old',access='append')
          write(11,2055) r(33,j),depth
 2055     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate34.dat',
     +      status='old',access='append')
          write(11,2056) r(34,j),depth
 2056     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate35.dat',
     +      status='old',access='append')
          write(11,2057) r(35,j),depth
 2057     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate36.dat',
     +      status='old',access='append')
          write(11,2058) r(36,j),depth
 2058     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate37.dat',
     +      status='old',access='append')
          write(11,2059) r(37,j),depth
 2059     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate38.dat',
     +      status='old',access='append')
          write(11,2060) r(38,j),depth
 2060     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate39.dat',
     +      status='old',access='append')
          write(11,2061) r(39,j),depth
 2061     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate40.dat',
     +      status='old',access='append')
          write(11,2062) r(40,j),depth
 2062     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate41.dat',
     +      status='old',access='append')
          write(11,2063) r(41,j),depth
 2063     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate42.dat',
     +      status='old',access='append')
          write(11,2064) r(42,j),depth
 2064     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate43.dat',
     +      status='old',access='append')
          write(11,2065) r(43,j),depth
 2065     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate44.dat',
     +      status='old',access='append')
          write(11,2066) r(44,j),depth
 2066     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate45.dat',
     +      status='old',access='append')
          write(11,2067) r(45,j),depth
 2067     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate46.dat',
     +      status='old',access='append')
          write(11,2068) r(46,j),depth
 2068     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate47.dat',
     +      status='old',access='append')
          write(11,2069) r(47,j),depth
 2069     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate48.dat',
     +      status='old',access='append')
          write(11,2070) r(48,j),depth
 2070     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate49.dat',
     +      status='old',access='append')
          write(11,2071) r(49,j),depth
 2071     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate50.dat',
     +      status='old',access='append')
          write(11,2072) r(50,j),depth
 2072     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate51.dat',
     +      status='old',access='append')
          write(11,2073) r(51,j),depth
 2073     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate52.dat',
     +      status='old',access='append')
          write(11,2074) r(52,j),depth
 2074     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate53.dat',
     +      status='old',access='append')
          write(11,2075) r(53,j),depth
 2075     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate54.dat',
     +      status='old',access='append')
          write(11,2076) r(54,j),depth
 2076     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate55.dat',
     +      status='old',access='append')
          write(11,2077) r(55,j),depth
 2077     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate56.dat',
     +      status='old',access='append')
          write(11,2078) r(56,j),depth
 2078     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate57.dat',
     +      status='old',access='append')
          write(11,2079) r(57,j),depth
 2079     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate58.dat',
     +      status='old',access='append')
          write(11,2080) r(58,j),depth
 2080     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate59.dat',
     +      status='old',access='append')
          write(11,2081) r(59,j),depth
 2081     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate60.dat',
     +      status='old',access='append')
          write(11,2082) r(60,j),depth
 2082     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate61.dat',
     +      status='old',access='append')
          write(11,2083) r(61,j),depth
 2083     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate62.dat',
     +      status='old',access='append')
          write(11,2084) r(62,j),depth
 2084     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate63.dat',
     +      status='old',access='append')
          write(11,2085) r(63,j),depth
 2085     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate64.dat',
     +      status='old',access='append')
          write(11,2086) r(64,j),depth
 2086     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate65.dat',
     +      status='old',access='append')
          write(11,2087) r(65,j),depth
 2087     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='xrate66.dat',
     +      status='old',access='append')
          write(11,2088) r(66,j),depth
 2088     format(1x,e14.7,2x,f12.4)
          close(11)
        if (j.eq.nx) then
        v_out = v_out+v_int
        endif
        endif
        if (time.eq.endt) then
          open(unit=11,file='Bfo1.inp',
     +      status='old',access='append')
          write(11,2089) sp(1,j),depth
 2089     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmo1.inp',
     +      status='old',access='append')
          write(11,2090) sp(2,j),depth
 2090     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='doc1.inp',
     +      status='old',access='append')
          write(11,2091) sp(3,j),depth
 2091     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='dox1.inp',
     +      status='old',access='append')
          write(11,2092) sp(4,j),depth
 2092     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Amm1.inp',
     +      status='old',access='append')
          write(11,2093) sp(5,j),depth
 2093     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifo1.inp',
     +      status='old',access='append')
          write(11,2094) sp(6,j),depth
 2094     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimo1.inp',
     +      status='old',access='append')
          write(11,2095) sp(7,j),depth
 2095     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfn1.inp',
     +      status='old',access='append')
          write(11,2096) sp(8,j),depth
 2096     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmn1.inp',
     +      status='old',access='append')
          write(11,2097) sp(9,j),depth
 2097     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='nitra1.inp',
     +      status='old',access='append')
          write(11,2098) sp(10,j),depth
 2098     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifn1.inp',
     +      status='old',access='append')
          write(11,2099) sp(11,j),depth
 2099     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimn1.inp',
     +      status='old',access='append')
          write(11,2100) sp(12,j),depth
 2100     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfs1.inp',
     +      status='old',access='append')
          write(11,2101) sp(13,j),depth
 2101     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bms1.inp',
     +      status='old',access='append')
          write(11,2102) sp(14,j),depth
 2102     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='sulpha1.inp',
     +      status='old',access='append')
          write(11,2103) sp(15,j),depth
 2103     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifs1.inp',
     +      status='old',access='append')
          write(11,2104) sp(16,j),depth
 2104     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bims1.inp',
     +      status='old',access='append')
          write(11,2105) sp(17,j),depth
 2105     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='POM1.inp',
     +      status='old',access='append')
          write(11,2106) sp(18,j),depth
 2106     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='tr1.inp',
     +      status='old',access='append')
          write(11,2107) sp(19,j),depth
 2107     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfa1.inp',
     +      status='old',access='append')
          write(11,2108) sp(20,j),depth
 2108     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bma1.inp',
     +      status='old',access='append')
          write(11,2109) sp(21,j),depth
 2109     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifa1.inp',
     +      status='old',access='append')
          write(11,2110) sp(22,j),depth
 2110     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bima1.inp',
     +      status='old',access='append')
          write(11,2111) sp(23,j),depth
 2111     format(1x,e14.7,2x,f12.4)
          close(11)
        endif
      end
