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


      subroutine split_scalar(n,nlocal,nstart,nend,node,nprocs)

      nstart = nint(float(node)/nprocs*n) + 1
      nend = nint(float(node+1)/nprocs*n)
      nlocal = nend - nstart + 1

      return
      end


      subroutine expand_setup(nlocal,ntotal,nstart,nend,nbyte,ndim,
     $     vec1,vec2,vec3,vec4,vec5,node,nprocs,ncube)
      integer vec1(*),vec2(*),vec3(*),vec4(*),vec5(*)

      ntotal = nlocal
      nstart = 1
      nend = nlocal

      return
      end


      subroutine expand(vec,vec1,vec2,vec3,vec4,vec5,
     $     node,nprocs,ncube)
      real*4 vec(*)
      integer vec1(*),vec2(*),vec3(*),vec4(*),vec5(*)

      return
      end


      subroutine fold(vec,vectmp,vec1,vec2,vec3,vec4,vec5,
     $     node,nprocs,ncube)
      real*4 vec(*),vectmp(*)
      integer vec1(*),vec2(*),vec3(*),vec4(*),vec5(*)

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
