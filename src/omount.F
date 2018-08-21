c---------------------------------------------------------------------
        subroutine omount(name,inew,istat,idev,path,lpath)
        character*(*) name,path
        character*100 patht
        character*1 quote
        parameter (quote='"')
        character*24 tapnam
        character*20 value,volnam,back
        character*4 vers
        character*60 fulnam
        character*200 mess
        character*1 str1
        dimension ihms(3)
#if ( ! defined(Machinel) )

        include 'oasstat.h'
   80   istat=0
        lname=istlen(name)
        nbys=1
        do while (nbys.le.lname.and.name(nbys:nbys).ne.'/')
          nbys=1+nbys
        enddo
        if(nbys.le.lname.and.nbys.gt.1) then
          volnam=name(1:nbys-1)
          lvolnam=nbys-1
          tapnam=name(nbys+1:lname)
          ltapnam=lname-nbys
        else
          lvolnam=0
          if(inew.ne.0) then
            write(0,*) 'No volume specified for new tape image'
            call exit(2)
          endif
          if(nbys.gt.1) then
            tapnam=name(1:lname)
            ltapnam=lname
          else
            tapnam=name(2:lname)
            ltapnam=lname-1
          endif
          knt=0
   42     call oassend('stat/tape/type=1 '
     1         //quote//tapnam(1:ltapnam)//quote,mess,lmess,-1)
          if(lmess.ge.2.and.mess(1:2).eq.'E-') then
            if(knt.gt.3) then
              write(0,*) mess(1:lmess)
              istat=4
              return
            else
              knt=knt+1
              call wait(300)
              goto 42
            endif
          endif
          if(lmess.le.0) then
            write(0,*) 'omount: Volume not found'
            call exit(2)
          endif
          volnam=mess(1:lmess)
          lvolnam=lmess
        endif
        knt=0
   41   call oassend('stat/vol/typ=9 '
     1     //quote//volnam(1:lvolnam)//quote//';',mess,lmess,-1)
        if(lmess.ge.2.and.mess(1:2).eq.'E-') then
          if(knt.gt.3) then
            write(0,*) mess(1:lmess)
            istat=4
            return
          else
            knt=knt+1
            call wait(300)
            goto 41
          endif
        endif
        back=mess(1:lmess)
        lback=lmess
        do i=1,MXJTAP
          if(oasvols(i).eq.back(1:lback)) then
            if(itchan(i).ge.-1) then
              if(itchan(i).ge.0) call cclose(itchan(i),ires,ierrno)
              itchan(i)=-1
              write(str1,'(i1)') i-1
              call oassend('user/def='//str1,mess,lmess,-1)
              call oassend('offline',mess,lmess,-1)
              call oassend('dismount',mess,lmess,-1)
              call oassend('dismount/v',mess,lmess,-1)
              oasvols(i)=char(0)
            else
              istat=-2
              write(0,*) 'omount: Back of volume in use; waiting ...'
              goto 88
            endif
          else if(oasvols(i).eq.volnam(1:lvolnam)) then
            if(itchan(i).ge.-1) then
              if(i.gt.noasdevs) then
                write(str1,'(i1)') i-1
                call oassend('user/def='//str1,mess,lmess,-1)
                call oassend('offline',mess,lmess,-1)
                call oassend('dismount',mess,lmess,-1)
                call oassend('dismount/v',mess,lmess,-1)
                oasvols(i)=char(0)
              endif
            else
              write(0,*) 'omount: Volume in use; waiting ...'
              istat=-1
              goto 88
            endif
          endif
        enddo
        iuser=-1
        iuser1=-1
        ismounted=0
        call itime(ihms)
        iup=mod(ihms(3),noasdevs)
        do i=1,noasdevs
          if(itchan(i).ge.-1) then
            if(iuser1.lt.0.or.iuser1.ne.iup)  iuser1=i-1
            if(oasvols(i).eq.volnam(1:lvolnam)) then
              iuser=i-1
              ismounted=1
            endif
          endif
        enddo
        if(iuser.lt.0) iuser=iuser1
        if(iuser.lt.0) then
          istat=-3
          write(0,*) 'omount: No available unit; waiting ...'
          goto 88
        endif
        write(str1,'(i1)') iuser
        call oassend('user/def='//str1,mess,lmess,-1)
        call oassend('offline',mess,lmess,-1)
        call oassend('dismount',mess,lmess,-1)
        call wait(700)
        if(tapnam(ltapnam:ltapnam).eq.'>') then
          vers=tapnam(ltapnam-3:ltapnam)
          lvers=4
          ltapnam=ltapnam-4
          if(vers(1:1).ne.'<') stop 'Invalid tape name verion'
        else 
          lvers=0
        endif

c       if(ismounted.eq.0) then
c         call oassend('dismount/vol',mess,lmess,-1)
c         write(0,*) 'mount '//quote//volnam(1:lvolnam)//quote//';'
c         call oassend('mount '//quote//volnam(1:lvolnam)//quote//';'
c    1            ,mess,lmess,-1)
c         if(lmess.ne.0) then
c           write(0,*) 'From OAS: ',mess(1:lmess)
c           call exit(2)
c         endif
c       endif

        if(lvers.gt.0) then
          fulnam=quote//volnam(1:lvolnam)//quote//';'//quote//tapnam(1:ltapnam)//quote//vers
        else
          fulnam=quote//volnam(1:lvolnam)//quote//';'//quote//tapnam(1:ltapnam)//quote
        endif
        lfulnam=istlen(fulnam)
          

        if(inew.eq.0) then
c         call glinkd(tapnam(1:ltapnam),value,lvalue,ilstat)
c         if(ilstat.ne.1) pause 'omount: unexpected status'
c         tapnam=value
c         ltapnam=lvalue
         

          write(0,*) 'mount '//fulnam(1:lfulnam)
          call oassend('mount '//fulnam(1:lfulnam),mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:21).eq.'E-DMGR-TAPE NOT FOUND') then
              istat=4
              goto 99
            else
              write(0,*) 'From OAS: ',mess(1:lmess)
              call exit(2)
            endif
          endif
c         call oassend('rewind d:',mess,lmess,-1)
c         if(lmess.ne.0) then
c           write(0,*) 'omount: rewinding:'//mess(1:lmess)
c           call exit(2)
c         endif
          call wait(500)


        else if(inew.eq.1) then

          write(0,'(a)') 'mount '//quote//volnam(1:lvolnam)//quote//';'
          call oassend('mount '//quote//volnam(1:lvolnam)//quote//';',mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS: ',mess(1:lmess)
            call exit(2)
          endif
          write(0,'(a)') 'mount/blank '//quote//tapnam(1:ltapnam)//quote

          call wait(1000)
          call oassend('mount/blank '//quote//tapnam(1:ltapnam)//quote,mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS: ',mess(1:lmess)
            call exit(2)
          endif



c         write(0,*) 'mount/blank '//fulnam(1:lfulnam)
c         call oassend('mount/blank '//fulnam(1:lfulnam),mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS: ',mess(1:lmess)
            if(mess(1:26).eq.'E-OMSG-NAME ALREADY EXISTS') then
              istat=3
              goto 99
            else
              call exit(2)
            endif
          endif
        else 
          write(0,*) 'opnflc: illegal inew for a tape image'
          call exit(2)
        endif
        call oassend('online',mess,lmess,-1)
        call wait(1000)
        if(loashost.gt.0) then
          path=oashost(1:loashost)//'/dev/nrst'//tapeid(iuser+1)
        else
          path='/dev/nrst'//tapeid(iuser+1)
        endif
        lpath=istlen(path)
        idev=iuser+1
        patht=path(1:lpath)
        write(0,*) 'Device: '//patht(1:lpath)
        goto 99
   88   continue
        do i=1,MXJTAP
          if(itchan(i).ge.0) then
            call cclose(itchan(i),ires,ierrno)
          endif
          call oasclose()
        enddo
        call wait(20000)
        call oasopen()
        goto 80
   99   continue
#else
      stop 'omount: not implemented in gfortran'
#endif
        return
        end
