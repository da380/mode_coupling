      subroutine rdsource(lu,iyr,imo,ida,iho,imi,fsec,
     #           tcmt,ityphdur,hdur,epla,eplo,dept,itypso,xm,ierr)
      dimension xm(6)
      character*24 ireg
      character*128 string
c
c---- first read hypocenter info
c
      read(lu,"(a4,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)") 
     #   isour,iyr,imo,ida,iho,imi,fsec,eplat,eplong,depth,xmb,xms,ireg
      write(6,"(a4,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)") 
     #   isour,iyr,imo,ida,iho,imi,fsec,eplat,eplong,depth,xmb,xms,ireg
c
      ierr=0
      ityphdur=2
      itypso=0
      do i=1,6
	xm(i)=0.
      enddo
c
      ios=0
      do while(ios.eq.0) 
      read(lu,"(a)",iostat=ios) string
      lstr=lnblnk(string)
      if(string(1:10).eq.'event name') then
      else if(string(1:10).eq.'time shift') then
        read(string(12:lstr),*) tcmt
      else if(string(1:13).eq.'half duration') then
        read(string(15:lstr),*) hdur
      else if(string(1:8).eq.'latitude') then
        read(string(10:lstr),*) epla
      else if(string(1:9).eq.'longitude') then
        read(string(11:lstr),*) eplo
      else if(string(1:5).eq.'depth') then
        read(string(7:lstr),*) dept
      else if(string(1:3).eq.'Mrr') then
        read(string(5:lstr),*) xm(1)
      else if(string(1:3).eq.'Mtt') then
        read(string(5:lstr),*) xm(2)
      else if(string(1:3).eq.'Mpp') then
        read(string(5:lstr),*) xm(3)
      else if(string(1:3).eq.'Mrt') then
        read(string(5:lstr),*) xm(4)
      else if(string(1:3).eq.'Mrp') then
        read(string(5:lstr),*) xm(5)
      else if(string(1:3).eq.'Mtp') then
        read(string(5:lstr),*) xm(6)
      else if(string(1:3).eq.'CFr') then
        itypso=1
        read(string(5:lstr),*) xm(1)
      else if(string(1:3).eq.'CFt') then
        itypso=1
        read(string(5:lstr),*) xm(2)
      else if(string(1:3).eq.'CFp') then
        itypso=1
        read(string(5:lstr),*) xm(3)
      else if(string(1:3).eq.'SFr') then
        itypso=2
        read(string(5:lstr),*) xm(1)
      else if(string(1:3).eq.'SFt') then
        itypso=2
        read(string(5:lstr),*) xm(2)
      else if(string(1:3).eq.'SFp') then
        itypso=2
        read(string(5:lstr),*) xm(3)
      endif
      enddo
      if(itypso.gt.1) then
        xm(4)=0.
        xm(5)=0.
        xm(6)=0.
      endif
      return
      end
