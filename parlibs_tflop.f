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
      include "mpif.h"

      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,node,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)
      ncube = 0
 10   if (2**ncube.lt.nprocs) then
        ncube = ncube + 1
        goto 10
      endif

      return
      end


      subroutine broadcast(data,n,node,nprocs,ncube)
      
      itype = 0
      iflag = 0
      if (node.eq.0) iflag = 1
      do i = 0,ncube-1
        if (iflag.eq.1) then
          idest = node + 2**i
          if (idest.lt.nprocs) ierr = nwrite(data,n,idest,itype,null)
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
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
          ierr = nread(itmp,ibyte,ipartner,itype,null)
          ivalue = ivalue + itmp
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(itmp,ibyte,ipartner,itype,null)
          ivalue = ivalue + itmp
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
            ierr = nread(itmp,ibyte,ipartner,itype,null)
            ivalue = ivalue + itmp
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(ivalue,ibyte,ipartner,itype,null)
        endif
        
      endif

      return
      end


      subroutine merge_imax(ivalue,node,nprocs,ncube)
      parameter (ibyte=4)
      integer ivalue,itmp

      itype = 0
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
          ierr = nread(itmp,ibyte,ipartner,itype,null)
          if (itmp.gt.ivalue) ivalue = itmp
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(itmp,ibyte,ipartner,itype,null)
          if (itmp.gt.ivalue) ivalue = itmp
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
            ierr = nread(itmp,ibyte,ipartner,itype,null)
            if (itmp.gt.ivalue) ivalue = itmp
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(ivalue,ibyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(ivalue,ibyte,ipartner,itype,null)
        endif
        
      endif

      return
      end


      subroutine merge_iv(vec,vectmp,n,node,nprocs,ncube)
      integer vec(*),vectmp(*)
      parameter (ibyte=4)

      itype = 0
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(vec,n*ibyte,ipartner,itype,null)
          ierr = nread(vectmp,n*ibyte,ipartner,itype,null)
          do j = 1,n
            vec(j) = vec(j) + vectmp(j)
          enddo
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(vec,n*ibyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(vectmp,n*ibyte,ipartner,itype,null)
          do j = 1,n
            vec(j) = vec(j) + vectmp(j)
          enddo
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(vec,n*ibyte,ipartner,itype,null)
            ierr = nread(vectmp,n*ibyte,ipartner,itype,null)
            do j = 1,n
              vec(j) = vec(j) + vectmp(j)
            enddo
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(vec,n*ibyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(vec,n*ibyte,ipartner,itype,null)
        endif
        
      endif

      return
      end


      subroutine merge_r(value,node,nprocs,ncube)
      real*4 value,tmp
      parameter (nbyte=4)

      itype = 0
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(value,nbyte,ipartner,itype,null)
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          value = value + tmp
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          value = value + tmp
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(value,nbyte,ipartner,itype,null)
            ierr = nread(tmp,nbyte,ipartner,itype,null)
            value = value + tmp
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(value,nbyte,ipartner,itype,null)
        endif
        
      endif

      return
      end


      subroutine merge_r8(value,node,nprocs,ncube)
      real*8 value,tmp
      parameter (nbyte=8)

      itype = 0
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(value,nbyte,ipartner,itype,null)
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          value = value + tmp
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          value = value + tmp
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(value,nbyte,ipartner,itype,null)
            ierr = nread(tmp,nbyte,ipartner,itype,null)
            value = value + tmp
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(value,nbyte,ipartner,itype,null)
        endif
        
      endif

      return
      end


      subroutine merge_r8max(value,node,nprocs,ncube)
      real*8 value,tmp
      parameter (nbyte=8)

      itype = 0
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(value,nbyte,ipartner,itype,null)
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          if (tmp.gt.value) value = tmp
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          if (tmp.gt.value) value = tmp
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(value,nbyte,ipartner,itype,null)
            ierr = nread(tmp,nbyte,ipartner,itype,null)
            if (tmp.gt.value) value = tmp
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(value,nbyte,ipartner,itype,null)
        endif
        
      endif

      return
      end


      subroutine merge_r8min(value,node,nprocs,ncube)
      real*8 value,tmp
      parameter (nbyte=8)

      itype = 0
      if (nprocs.eq.2**ncube) then

        do i = 0,ncube-1
          ipartner = ieor(node,2**i)
          ierr = nwrite(value,nbyte,ipartner,itype,null)
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          if (tmp.lt.value) value = tmp
        enddo

      else

        nhalf = 2**(ncube-1)
        nextra = nprocs - nhalf
        if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nread(tmp,nbyte,ipartner,itype,null)
          if (tmp.lt.value) value = tmp
        endif
        if (node.lt.nhalf) then
          do i = 0,ncube-2
            ipartner = ieor(node,2**i)
            ierr = nwrite(value,nbyte,ipartner,itype,null)
            ierr = nread(tmp,nbyte,ipartner,itype,null)
            if (tmp.lt.value) value = tmp
          enddo
        endif
        if (node.lt.nextra) then
          ipartner = node + nhalf
          ierr = nwrite(value,nbyte,ipartner,itype,null)
        else if (node.ge.nhalf) then
          ipartner = node - nhalf
          ierr = nread(value,nbyte,ipartner,itype,null)
        endif
        
      endif

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
      
      ix = mod(node,nx)
      iy = mod(node/nx,ny)
      iz = node/nx/ny

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

      iwest = nx*ny*iz + nx*iy + iwest
      ieast = nx*ny*iz + nx*iy + ieast
      isouth = nx*ny*iz + nx*isouth + ix
      inorth = nx*ny*iz + nx*inorth + ix
      idown = nx*ny*idown + nx*iy + ix
      iup = nx*ny*iup + nx*iy + ix

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
        if (ipartner.lt.nprocs) then
          ierr = nwrite(tmp,0,ipartner,itype,null)
          ierr = nread(tmp,0,ipartner,itype,null)
        endif
      enddo
      
      return
      end


      subroutine timer(time)
      real*8 time,dclock
      
      time = dclock()
      
      return
      end
