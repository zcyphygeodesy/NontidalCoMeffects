!  NontidalCoMeffects.f90 
!
!  FUNCTIONS:
!  NontidalCoMeffects - Entry point of console application.
!
!****************************************************************************

      program NontidalCoMeffects
      implicit none
	character*800::line,str,astr,stm
	integer i,j,k,nn,IY,IM,ID,ih,imn,kk,sn,kln
      real*8 mjd,mjd1,mjd2,sec,rec(800),pi,RAD,mdl,xx,yy,zz,rr
	real*8 GRS(6),BLH(3),rln(3),r0,ae,gm,fh1,fl1,tdn(14),tdn0(14),stmjd
	real*8 tm0,tm1,tm2,tm,GC(5000,4),mjd01,mjd02
	real*8 cosst,sinst,cosla,sinla,NFD(5),gr,tmp
	real*8 plon, plat, phgt, bgntm, endtm, tmdlt,ltm
	integer::status=0
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
      pi=datan(1.d0)*4.d0; RAD=pi/180.d0;r0=6371.d3;gm=GRS(1);ae=GRS(2)
      fh1=-0.2871129880d0;fl1=0.1045044062d0
      !Input longitude (degree decimal), latitude (degree decimal), ellipsoidal height (m) 
      !输入经纬度大地高
	plon=107.23d0; plat=29.91d0; phgt=72.4d0
      !Input starting time (long integer time), ending time, time interval (minute) 输入起止时间与时间间隔
      bgntm=20140101;endtm=20160101; tmdlt=720.d0/1440.d0
 	!Read the Earth's mass centric variation time series file from the measured SLR
      open(unit=8,file="Monthly_geocenter_MK.txt",status="old",iostat=status)
      if(status/=0) goto 902 
      do i=1,15
         read(8,'(a)') line
      enddo
	i=1
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickRecord(line,kln,rec,sn)
         if(sn<4)goto 407
         if(rec(2)>299.d0.or.rec(3)>299.d0.or.rec(4)>299.d0)goto 407
         GC(i,1:3)=rec(2:4)*1.d-3
         IY=floor(rec(1))
         call CAL2JD (IY,1,1,mjd,j) 
         GC(i,4)=mjd+(rec(1)-dble(IY))*365.25d0+15.d0
         i=i+1
407      continue
	enddo
903   close(8)
      nn=i-1
      if(nn<5) goto 902 
      mdl=(GC(nn-1,4)-GC(nn-5,4))/4.d0;kk=0
      tm1=GC(1,4);tm2=GC(nn,4)
      !Transform the long integer time (date) agreed by ETideLoad to year, month, day, hour, minute and second
      !ETideLoad格式日期tm转年月日时分秒。
      call tmcnt(bgntm,IY,IM,ID,ih,imn,sec)
      !Gregorian Calendar to Julian Date.
      call CAL2JD (IY,IM,ID,mjd,j)
      mjd01=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
      call tmcnt(endtm,IY,IM,ID,ih,imn,sec)
      call CAL2JD (IY,IM,ID,mjd,j)
      mjd02=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
      if(mjd02<mjd01)goto 902
      open(unit=10,file='reslt.txt',status="replace") !Output file 输出文件
      BLH(1)=plat;BLH(2)=plon;BLH(3)=phgt
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      write(10,'(a8,2F12.6,F10.3,F15.6)')'calcpnt',plon,plat,phgt,mjd01
      mjd=mjd01;mjd02=mjd02+1.d-6;k=0
      do while(mjd<mjd02)
         !Interpolate the Earth's centric variation of mass at mjd epoch time
         !内插mjd历元时刻地球质心
         call spline1d(mjd,xx,GC(1:nn,4),GC(1:nn,1),nn)
         call spline1d(mjd,yy,GC(1:nn,4),GC(1:nn,2),nn)
         call spline1d(mjd,zz,GC(1:nn,4),GC(1:nn,3),nn)
         cosst=dsin(rln(2)*RAD);sinst=dcos(rln(2)*RAD)
         cosla=dcos(rln(3)*RAD);sinla=dsin(rln(3)*RAD)
         !计算各种大地测量要素的地球质心变化效应
         !Calculate the non-tidal Earth's centric effects of mass on all-element geodetic variations
         tdn(1:14)=0.d0;tmp=gm/r0*ae
         tdn(1)=tmp/rr**2/gr*(xx*cosla*sinst+yy*sinla*sinst+zz*cosst)
         tdn(2)=2.d0*tmp/rr**3*(1.d0+2.d0*fh1)*(xx*cosla*sinst+yy*sinla*sinst+zz*cosst)
         tdn(3)=2.d0*tmp/rr**3*(xx*cosla*sinst+yy*sinla*sinst+zz*cosst)
         tdn(4)=tmp/rr**3/gr*sinst*(1.d0-fh1)*(xx*cosla*cosst+yy*sinla*sinst-zz*sinst)
         tdn(5)=tmp/rr**3/gr*(1.d0-fh1)*(xx*sinla-yy*cosla)
         tdn(6)=tmp/rr**3/gr*sinst*(xx*cosla*cosst+yy*sinla*sinst-zz*sinst)
         tdn(7)=tmp/rr**3/gr*(xx*sinla-yy*cosla)
         tdn(8)=-tmp/rr**2/gr*fl1*(xx*sinla-yy*cosla)
         tdn(9)=-tmp/rr**2/gr*sinst*fl1*(xx*cosla*cosst+yy*sinla*sinst-zz*sinst)
         tdn(10)=tmp/rr**2/gr*fh1*(xx*cosla*sinst+yy*sinla*sinst+zz*cosst)
         tdn(12)=6.d0*tmp/rr**4*(xx*cosla*sinst+yy*sinla*sinst+zz*cosst)
         tdn(13)=tmp/rr**4*(xx*cosla*sinst+yy*sinla*sinst+zz*cosst)
         tdn(14)=tmp/rr**4*(xx*sinla+yy*cosla)
         tdn(11)=tdn(10)-tdn(1)
	   tdn(1)=tdn(1)*1.d3
	   tdn(2:3)=tdn(2:3)*1.0e8
	   tdn(4:7)=tdn(4:7)/RAD*36.d5
	   tdn(8:11)=-tdn(8:11)*1.d3
	   tdn(12:14)=tdn(12:14)*1.d14 !0.01mE
         if(k<1)tdn0(1:14)=tdn(1:14)
         call mjdtotm(mjd,ltm)
         call tmtostr(ltm,stm)
         xx=xx*1.d3;yy=yy*1.d3;zz=zz*1.d3
         write(10,'(a15,F12.6,40F12.4)')adjustl(stm), mjd-mjd01, xx,yy,zz,(tdn(i),i=1,14)!(tdn(i)-tdn0(i),i=1,14)
         mjd=mjd+tmdlt;k=k+1
906      continue
	enddo
      close(10)
902	continue
101   format(a,40F12.4)
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
!
!*******************************************************
!
      subroutine spline1d(x,y,xa,ya,n)
        !不规则采样的1维样条插值,返回y,xa单调
        integer n
	  real*8 x,y,xa(n),ya(n),sx(1),sy(1),dsy(1)
!-------------------------------------------------------------
        sx(1)=x;sy(1)=0.d0
        call spline3(xa,ya,n,sx,sy,dsy,1)
        y=sy(1)
        if(dabs(y)<1.d-16)y=99999.d0
        end
!
!*******************************************************
!
      subroutine spline3( x, y, n, sx, f, f1, m)  !三次样条插值程序
      implicit none
      integer:: i, j, k, m, n, n1, n2
      integer, parameter :: dp = 8
      Real(dp):: x(n), y(n), sx(m), f(m), f1(m)
      Real(dp):: s2(n), h(n-1), dy(n-1), s(n-1), e(n-2)
      Real(dp):: z, h1, h2, h3, h4
!------------------------------------------------------------- 
      n1=n-1;n2=n-2
      do i = 1, n1
        h(i)  = x(i+1) - x(i)
        dy(i) = (y(i+1) - y(i))/h(i)
      end do
      s2(1) = 0.d0; s2(n) = 0.d0
      do i = 2, n1
        s2(i) = 6.d0 * (dy(i) - dy(i-1) )
      end do
      z = 0.5d0 / ( h(1) + h(2) )
      s(1) = -h(2) * z
      e(1) = s2(2) * z
      do i = 2, n2
        k = i - 1
        j = i + 1
        z = 1.d0 / ( 2.d0*( h(i)+h(j) ) + h(i)*s(k) )
        s(i) = -h(j) * z
        e(i) = ( s2(j)-h(i)*e(k) ) * z
      end do
      s2(n1) = e(n2)
      do i = n2, 2, -1
        k = i - 1
        s2(i) = s(k)*s2(i+1) + e(k)
      end do
      do i = 1, n1
        s(i) = ( s2(i+1) - s2(i) ) / h(i)
      end do
      i = 2; k = 1
      do j = 1, m
        do 
          if ( i<n) then
            if ( sx(j) > x(i) ) then
              k = i; i = i + 1
            else
              exit
            endif
          else
            goto 333
          end if
        end do
        h1 = sx(j) - x(k); h2 = sx(j) - x(i)
        h3 = h1 * h2; h4 = s2(k) + h1*s(k)
        z = ( s2(i) + s2(k) + h4 ) / 6.d0
        f(j) = y(k) + h1*dy(k) + h3*Z
        f1(j) = dy(k) + z*( h1+h2 ) + h3 * s(k) / 6.d0
333     continue
      end do
      end
