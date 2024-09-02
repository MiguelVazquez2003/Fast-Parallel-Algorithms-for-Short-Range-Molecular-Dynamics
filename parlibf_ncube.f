c Parallel Molecular Dynamics Algorithms
c
c Authored by Steve Plimpton
c   (505) 845-7873, sjplimp@cs.sandia.gov
c   Dept 1421, MS 1111, Sandia National Labs, Albuquerque, NM  87185-1111
c
c Last modified:  10 March 1995
c
c See the README file for more information
c

      subroutine whatami(node,nprocs,ncube)

      call whoami(node,i,j,ncube)
      nprocs = 2**ncube

      return
      end


      subroutine broadcast(data,n,node,nprocs,ncube)
      
      iflag = 0
      if (node.eq.0) iflag = 1
      do i = 0,ncube-1
        if (iflag.eq.1) then
          idest = node + 2**i
          call big_write(data,n,idest)
        else if (node.lt.2**(i+1)) then
          iflag = 1
          isrc = node - 2**i
          call big_read(data,n,isrc)
        endif
      enddo

      return
      end


      subroutine big_read(data,n,isrc)
      integer data(*)
      parameter (max_size=16384)
      
      itype = 0
      nleft = n
      ipnt = 1
 10   isize = min(nleft,max_size)
      nleft = nleft - isize
      ierr = nread(data(ipnt),isize,isrc,itype,0)
      if (nleft.eq.0) return
      ierr = nwrite(tmp,0,isrc,itype,null)
      ipnt = ipnt + max_size/4
      goto 10

      return
      end


      subroutine big_write(data,n,idest)
      integer data(*)
      parameter (max_size=16384)
      
      itype = 0
      nleft = n
      ipnt = 1
 10   isize = min(nleft,max_size)
      nleft = nleft - isize
      ierr = nwrite(data(ipnt),isize,idest,itype,null)
      if (nleft.eq.0) return
      ierr = nread(tmp,0,idest,itype,null)
      ipnt = ipnt + max_size/4
      goto 10

      return
      end
      

      subroutine merge_i(ivalue,node,nprocs,ncube)
      integer ivalue,itmp
      parameter (ibyte=4)

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
        ierr = nread(itmp,ibyte,ipartner,itype,null)
        ivalue = ivalue + itmp
      enddo

      return
      end


      subroutine merge_imax(ivalue,node,nprocs,ncube)
      parameter (ibyte=4)
      integer ivalue,itmp

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
        ierr = nread(itmp,ibyte,ipartner,itype,null)
        if (itmp.gt.ivalue) ivalue = itmp
      enddo

      return
      end


      subroutine merge_iv(vec,vectmp,n,node,nprocs,ncube)
      integer vec(*),vectmp(*)
      parameter (ibyte=4)

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(vec,n*ibyte,ipartner,itype,null)
        ierr = nread(vectmp,n*ibyte,ipartner,itype,null)
	do j = 1,n
          vec(j) = vec(j) + vectmp(j)
        enddo
      enddo

      return
      end


      subroutine merge_r(value,node,nprocs,ncube)
      real*4 value,tmp
      parameter (nbyte=4)

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(value,nbyte,ipartner,itype,null)
        ierr = nread(tmp,nbyte,ipartner,itype,null)
        value = value + tmp
      enddo

      return
      end


      subroutine merge_r8(value,node,nprocs,ncube)
      real*8 value,tmp
      parameter (nbyte=8)

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(value,nbyte,ipartner,itype,null)
        ierr = nread(tmp,nbyte,ipartner,itype,null)
        value = value + tmp
      enddo

      return
      end


      subroutine merge_r8max(value,node,nprocs,ncube)
      real*8 value,tmp
      parameter (nbyte=8)

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(value,nbyte,ipartner,itype,null)
        ierr = nread(tmp,nbyte,ipartner,itype,null)
        if (tmp.gt.value) value = tmp
      enddo

      return
      end


      subroutine merge_r8min(value,node,nprocs,ncube)
      real*8 value,tmp
      parameter (nbyte=8)

      itype = 0
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(value,nbyte,ipartner,itype,null)
        ierr = nread(tmp,nbyte,ipartner,itype,null)
        if (tmp.lt.value) value = tmp
      enddo

      return
      end


      subroutine stats_r8(data,n,ave,xmax,xmin,
     $     histo,histotmp,nhisto,node,nprocs,ncube)
      implicit real*8 (a-h,o-z)
      real*8 data(*),ave,xmax,xmin
      integer histo(*),histotmp(*)
      
      xmin = 1.0E20
      xmax = -1.0E20
      ave = 0.0
      do i = 1,n
        ave = ave + data(i)
        if (data(i).lt.xmin) xmin = data(i)
        if (data(i).gt.xmax) xmax = data(i)
      enddo

      ntot = n
      call merge_i(ntot,node,nprocs,ncube)
      ave = ave/ntot
      call merge_r8(ave,node,nprocs,ncube)
      call merge_r8max(xmax,node,nprocs,ncube)
      call merge_r8min(xmin,node,nprocs,ncube)

      do i = 1,nhisto
        histo(i) = 0
      enddo
      
      del = xmax-xmin
      do i = 1,n
        if (del.eq.0.0) then
          j = 1
        else
          j = (data(i)-xmin)/del * nhisto + 1
          if (j.gt.nhisto) j = nhisto
        endif
        histo(j) = histo(j) + 1
      enddo
      
      call merge_iv(histo,histotmp,nhisto,node,nprocs,ncube)

      return
      end


      subroutine split_scalar(n,nlocal,nstart,nend,node,nprocs)

      nstart = nint(float(node)/nprocs*n) + 1
      nend = nint(float(node+1)/nprocs*n)
      nlocal = nend - nstart + 1

      return
      end


      subroutine expand_setup(nlocal,ntotal,nstart,nend,nbyte,ndim,
     $     vec1,vec2,vec3,vec4,vec5,node,nprocs,ncube)
      integer vec1(*),vec2(*),vec3(*),vec4(*),vec5(*)
      parameter (ibyte=4)
      
      ilocal = nlocal
      istart = 1
      iend = nlocal
      itype = 0
      
      do i = 1,ncube
        
	ipartner = ieor(node,2**(i-1))

        vec1(i) = istart
        vec2(i) = ilocal
        
        ierr = nwrite(ilocal,ibyte,vec5(i),itype,null)
        ierr = nread(jlocal,ibyte,vec5(i),itype,null)

        if (ipartner.lt.node) then
          istart = istart - jlocal
          vec3(i) = istart
        else
          vec3(i) = iend + 1
          iend = iend + jlocal
        endif

        vec4(i) = jlocal
        ilocal = ilocal + jlocal

      enddo
      
      ntotal = ilocal
      nshift = 1 - istart
      nstart = 1 + nshift
      nend = nlocal + nshift

      do i = 1,ncube
        vec1(i) = (vec1(i)+nshift-1)*ndim + 1
        vec2(i) = vec2(i) * ndim*nbyte
        vec3(i) = (vec3(i)+nshift-1)*ndim + 1
        vec4(i) = vec4(i) * ndim*nbyte
      enddo

      return
      end


      subroutine expand(vec,vec1,vec2,vec3,vec4,vec5,
     $     node,nprocs,ncube)
      real*4 vec(*)
      integer vec1(*),vec2(*),vec3(*),vec4(*),vec5(*)
      parameter (max_size=16384)
      
      itype = 0
      do i = 1,ncube
        iwrite = vec1(i)
        islen = vec2(i)
        iread = vec3(i)
        irlen = vec4(i)
        ipartner = vec5(i)
        if (islen.le.max_size.and.irlen.le.max_size) then
          ierr = nwrite(vec(iwrite),islen,ipartner,itype,null)
          ierr = nread(vec(iread),irlen,ipartner,itype,null)
        else
          call swap(-1,vec(iwrite),islen,ipartner,
     $         vec(iread),irlen,ipartner,0,1,0,0)
        endif
      enddo

      return
      end


      subroutine fold(vec,vectmp,vec1,vec2,vec3,vec4,vec5,
     $     node,nprocs,ncube)
      real*4 vec(*),vectmp(*)
      integer vec1(*),vec2(*),vec3(*),vec4(*),vec5(*)
      parameter (max_size=16384)
      parameter (nbyte=4)

      itype = 0
      do i = ncube,1,-1
        iread = vec1(i)
        irlen = vec2(i)
        iwrite = vec3(i)
        islen = vec4(i)
        ipartner = vec5(i)
        if (islen.le.max_size.and.irlen.le.max_size) then
          ierr = nwrite(vec(iwrite),islen,ipartner,itype,null)
          ierr = nread(vectmp,irlen,ipartner,itype,null)
        else
          call swap(-1,vec(iwrite),islen,ipartner,
     $         vectmp,irlen,ipartner,0,1,0,0)
        endif
        do j = 1,irlen/nbyte
          vec(iread) = vec(iread) + vectmp(j)
          iread = iread + 1
        enddo
      enddo

      return
      end


      subroutine swap(node,sbuf,islen,isnode,rbuf,irlen,irnode,
     $     iunknownflag,ihandflag,isyncflag,icopyflag)
      integer sbuf(*),rbuf(*)
      parameter (max_size=16384)
      parameter (ibyte=4)

      if (isnode.ne.node) then

        if (ihandflag.eq.0) then

          itype = 0
          ierr = nwrite(sbuf,islen,isnode,itype,null)
          if (ierr.ne.0) write (6,*) 'Error: Swap overflow',node,
     $         isnode,islen,ierr
          inum = nread(rbuf,irlen,irnode,itype,null)
          if (isyncflag.eq.1) then
            ierr = nwrite(tmp,0,irnode,itype,null)
            ierr = nread(tmp,0,isnode,itype,null)
          endif
          
        else

          itype = 0
          nsleft = islen
          nrleft = irlen
          inum = 0
          isflag = 1
          irflag = 1
          ipnt = 1
 10       if (isflag.eq.1) then
            isize = min(nsleft,max_size)
            ierr = nwrite(sbuf(ipnt),isize,isnode,itype,null)
            nsleft = nsleft - isize
            if (isize.lt.max_size) isflag = 0
          endif
          if (irflag.eq.1) then
            isize = min(nrleft,max_size)
            iactual = nread(rbuf(ipnt),isize,irnode,itype,null)
            nrleft = nrleft - isize
            if (iactual.lt.max_size) irflag = 0
            inum = inum + iactual
          endif
          ipnt = ipnt + max_size/ibyte
          if (isnode.ne.irnode.or.isflag.ne.irflag) then
            if (irflag.eq.1) ierr = nwrite(tmp,0,irnode,itype,null)
            if (isflag.eq.1) ierr = nread(tmp,0,isnode,itype,null)
          endif
          if (isflag.eq.1.or.irflag.eq.1) goto 10

        endif

      else if (icopyflag.eq.1) then

        isize = islen
        if (iunknownflag.eq.1) isize = min(islen,irlen)
        do i = 1,isize/ibyte
          rbuf(i) = sbuf(i)
        enddo
        inum = isize

      else

        inum = 0

      endif

      if (iunknownflag.eq.1) irlen = inum

      return
      end


      subroutine synchro(node,nprocs,ncube)
      
      itype = 32767
      do i = 0,ncube-1
        ipartner = ieor(node,2**i)
        ierr = nwrite(tmp,0,ipartner,itype,null)
        ierr = nread(tmp,0,ipartner,itype,null)
      enddo
      
      return
      end


      subroutine timer(time)
      real*8 time
      
      nstart = n_time(nstart_micro)
      time = dble(nstart) + dble(nstart_micro)/1.0D6
      
      return
      end


      subroutine exit(i)
      
      stop
      end
