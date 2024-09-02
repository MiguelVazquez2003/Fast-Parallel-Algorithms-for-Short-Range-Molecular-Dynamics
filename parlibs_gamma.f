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

      node = mynode()
      ncube = nodedim()
      nprocs = 2**ncube

      return
      end


      subroutine broadcast(data,n,node,nprocs,ncube)
      
      itype = 0
      iflag = 0
      if (node.eq.0) iflag = 1
      do i = 0,ncube-1
        if (iflag.eq.1) then
          idest = node + 2**i
          ierr = nwrite(data,n,idest,itype,null)
        else if (node.lt.2**(i+1)) then
          iflag = 1
          isrc = node - 2**i
          ierr = nread(data,n,isrc,itype,null)
        endif
      enddo

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


      subroutine mesh_3d(nx,ny,nz,node,ix,iy,iz,
     $     iwest,ieast,isouth,inorth,idown,iup)
      integer seq_to_gray,gray_to_seq
      
      ixseq = mod(node,nx)
      iyseq = mod(node/nx,ny)
      izseq = node/nx/ny

      ix = gray_to_seq(ixseq)
      iy = gray_to_seq(iyseq)
      iz = gray_to_seq(izseq)

      iwest = ix - 1
      if (iwest.eq.-1) iwest = nx - 1
      ieast = ix + 1
      if (ieast.eq.nx) ieast = 0
      isouth = iy - 1
      if (isouth.eq.-1) isouth = ny - 1
      inorth = iy + 1
      if (inorth.eq.ny) inorth = 0
      idown = iz - 1
      if (idown.eq.-1) idown = nz - 1
      iup = iz + 1
      if (iup.eq.nz) iup = 0

      iwest = nx*ny*izseq + nx*iyseq + seq_to_gray(iwest)
      ieast = nx*ny*izseq + nx*iyseq + seq_to_gray(ieast)
      isouth = nx*ny*izseq + nx*seq_to_gray(isouth) + ixseq
      inorth = nx*ny*izseq + nx*seq_to_gray(inorth) + ixseq
      idown = nx*ny*seq_to_gray(idown) + nx*iyseq + ixseq
      iup = nx*ny*seq_to_gray(iup) + nx*iyseq + ixseq

      return
      end


      integer function seq_to_gray(seq)
      integer seq
      
      seq_to_gray = ieor(seq,seq/2)

      return
      end
      
      
      integer function gray_to_seq(gray)
      integer gray,rsgray
      
      gray_to_seq = gray
      rsgray = gray/2
 10   gray_to_seq = ieor(gray_to_seq,rsgray)
      rsgray = rsgray/2
      if (rsgray.ne.0) goto 10

      return
      end


      subroutine swap(node,sbuf,islen,isnode,rbuf,irlen,irnode,
     $     iunknownflag,ihandflag,isyncflag,icopyflag)
      integer sbuf(*),rbuf(*)
      parameter (ibyte=4)

      if (isnode.ne.node) then

        itype = 0
        ierr = nwrite(sbuf,islen,isnode,itype,null)
        inum = nread(rbuf,irlen,irnode,itype,null)
        if (isyncflag.eq.1) then
          ierr = nwrite(tmp,0,irnode,itype,null)
          ierr = nread(tmp,0,isnode,itype,null)
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
      real*8 time,dclock
      
      time = dclock()
      
      return
      end


      integer function nread(buf,isize,ipartner,itype,null)
      parameter (max_type=32768)

      iflag = ipartner*max_type + itype
      call crecv(iflag,buf,isize) 
      nread = infocount()

      return
      end


      integer function nwrite(buf,isize,ipartner,itype,null)
      parameter (max_type=32768)

      iflag = mynode()*max_type + itype
      call csend(iflag,buf,isize,ipartner,0)
      nwrite = 0

      return
      end
