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

      node = 0
      nprocs = 1
      ncube = 0
     
      return
      end


      subroutine broadcast(data,n,node,nprocs,ncube)

      return
      end


      subroutine merge_i(ivalue,node,nprocs,ncube)

      return
      end


      subroutine merge_imax(ivalue,node,nprocs,ncube)

      return
      end


      subroutine merge_iv(vec,vectmp,n,node,nprocs,ncube)
      integer vec(*),vectmp(*)

      return
      end


      subroutine merge_r(value,node,nprocs,ncube)

      return
      end


      subroutine merge_r8(value,node,nprocs,ncube)
      real*8 value

      return
      end


      subroutine merge_r8max(value,node,nprocs,ncube)
      real*8 value

      return
      end


      subroutine merge_r8min(value,node,nprocs,ncube)
      real*8 value

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

      ix = 0
      iy = 0
      iz = 0

      iwest = node
      ieast = node
      isouth = node
      inorth = node
      idown = node
      iup = node

      return
      end


      subroutine swap(node,sbuf,islen,isnode,rbuf,irlen,irnode,
     $     iunknownflag,ihandflag,isyncflag,icopyflag)
      integer sbuf(*),rbuf(*)
      parameter (ibyte=4)

      if (icopyflag.eq.1) then

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
      
      return
      end


      subroutine timer(time)
      real*8 time
      real array(2)
      
      r = etime(array)
      time = array(1)
      
      return
      end
