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

c force-decomposition algorithm

      program ljf

      include 'ljf.h'
      real*8 time1,time2

 10   call whatami(node,nprocs,ncube)
      
      call input

      if (node.eq.0) write (6,*) 'Setting up ...'

      call setup_general

      call setup_parallel
      if (ierrorflag.lt.0) goto 20

      call setup_atom

      call neighbor
      call status

      if (node.eq.0) write (6,*) 'Starting dynamics ...'

      call synchro(node,nprocs,ncube)
      call timer(time1)
      
      do i = 1,ntimes
        call integrate(i)
        if (mod(i,nstat).eq.0) call status
      enddo
      
      call synchro(node,nprocs,ncube)
      call timer(time2)
      
      time_total = time2 - time1
      call output
      
 20   continue

      if (node.eq.0) then
        write (6,*) '(0) Stop, (1) Continue'
        read (5,*) iflag
      endif
      
      call broadcast(iflag,ibyte,node,nprocs,ncube)
      if (iflag.eq.1) goto 10
      
      end
      
      
c *************************
c Subroutines
c *************************
      
c input from file and prompts
      
      subroutine input
      include 'ljf.h'

      if (node.eq.0) then

        open (unit=1,file='lj.in',status='old')
        read (1,*)
        read (1,*)
        read (1,*) dtstar
        read (1,*) tstar
        read (1,*) rhostar
        read (1,*) nneigh
        read (1,*) rc,rs
        read (1,*) nstat
        read (1,*) iseed
        close (1)

        write (6,*) 'Atoms: nx,ny,nz'
        read (5,*) nxsize,nysize,nzsize
        write (6,*) '# of Timesteps'
        read (5,*) ntimes
        write (6,*) 'Neighboring: (0) N^2 (1) Binned'
        read (5,*) ineigh
        if (ineigh.eq.1) then
          write (6,*) 'Specify binning: (0) No (1) Yes'
          read (5,*) ibin
          if (ibin.eq.1) then
            write (6,*) 'Bins in each direction: x,y,z'
            read (5,*) nbinx,nbiny,nbinz
          endif
        endif
        write (6,*) 'Newtons 3rd law: (0) No (1) Yes'
        read (5,*) inewt
        write (6,*) 'Transpose: (0) No (1) Yes'
        read (5,*) itrans
        if (nprocs.eq.2**ncube) then
          write (6,*) 'Specify processor grid: (0) No (1) Yes'
          read (5,*) igrid
        else
          igrid = 1
        endif
        if (igrid.eq.1) then
          write (6,*) 'Processor grid: numrow,numcol'
          read (5,*) numrow,numcol
        endif

      endif

      inum = 17*ibyte + 5*nbyte
      call broadcast(dtstar,inum,node,nprocs,ncube)

      return
      end
        

c setup program units

      subroutine setup_general
      include 'ljf.h'

      natoms = nxsize*nysize*nzsize*4
      xprd = 1.0
      xprd2 = 0.5
      yprd = float(nysize)/nxsize
      yprd2 = yprd/2.0
      zprd = float(nzsize)/nxsize
      zprd2 = zprd/2.0
      alat = 1.0 / nxsize
      sigma = alat / (4.0/rhostar)**(1.0/3.0)
      dt = dtstar*sigma
      dtforce = 48.0*dtstar*sigma
      sigsq = sigma*sigma
   
      cutsq1 = (rc*sigma)*(rc*sigma)
      cutsq2 = (rs*sigma)*(rs*sigma)
      xmc = xprd - rc*sigma
      ymc = yprd - rc*sigma
      zmc = zprd - rc*sigma
      xms = xprd - rs*sigma
      yms = yprd - rs*sigma
      zms = zprd - rs*sigma
      
      time_force = 0.0D0
      time_neigh = 0.0D0
      time_comm = 0.0D0

      mstat = 0
      mneigh = 0
      neighmax = 0

      return
      end
      
      
c setup neighbor bins, force-decomposition
c  error check

      subroutine setup_parallel
      include 'ljf.h'
      parameter (small=1.0E-6)
      
c set neighbor bin size
c  addding small -> bins slightly larger
c  prevents round-off error (bin # ix = nbinx) when atoms are binned
      
      if (ineigh.eq.1) then
	if (ibin.eq.0) then
          nbinx = xprd / (rs*sigma)
          nbiny = yprd / (rs*sigma)
          nbinz = zprd / (rs*sigma)
        endif
        binsizex = (xprd+small*xprd) / nbinx
        binsizey = (yprd+small*yprd) / nbiny
        binsizez = (zprd+small*zprd) / nbinz
      endif
      
c map processors to grid in calendar order
c  reset nprocs,ncube to only used processors, map out unused processors
c  allows one to run on less nodes than allocated (so can still allocate a 
c  rectangular region)

      if (igrid.eq.0) then
        numrow = 2 ** ((ncube+1)/2)
        numcol = 2 ** ((ncube+0)/2)
      endif
      
      if (numrow*numcol.gt.nprocs) then
        call error('Bad grid specification')
        return
      endif
      
      nprocs = numrow*numcol
      ncube = 0
 10   if (2**ncube.lt.nprocs) then
        ncube = ncube + 1
        goto 10
      endif
      if (node.ge.nprocs) then
        ierrorflag = -1
        return
      endif
      
      irow = node/numcol
      icol = mod(node,numcol)

      numrowdim = 0
 20   if (2**numrowdim.lt.numrow) then
        numrowdim = numrowdim + 1
        goto 20
      endif
      numcoldim = 0
 30   if (2**numcoldim.lt.numcol) then
        numcoldim = numcoldim + 1
        goto 30
      endif
      
c setup # of atoms on each processor
      
      call split_scalar(natoms,nlocal,i,j,node,nprocs)
      
c setup communication within each row (across columns)
      
      do i = 1,numcoldim
        icolnew = ieor(icol,2**(i-1))
        vec5(i) = irow*numcol + icolnew
      enddo
      call expand_setup(nlocal,nglocal,nstart,nend,nbyte,3,
     $     vec1,vec2,vec3,vec4,vec5,icol,numcol,numcoldim)
      
c setup transpose if necessary
      
      if (itrans.eq.0) then
        itranrecv = node
        itransend = node
        nolocal = nlocal
      else
        itranrecv = icol*numrow + irow
        irowtmp = mod(node,numrow)
        icoltmp = node / numrow
        itransend = irowtmp*numcol + icoltmp
        call split_scalar(natoms,nolocal,i,j,itranrecv,nprocs)
      endif
      
c setup communication within each column (across rows)
      
      do i = 1,numrowdim
        irownew = ieor(irow,2**(i-1))
        vec10(i) = irownew*numcol + icol
      enddo
      call expand_setup(nolocal,ngolocal,nostart,noend,nbyte,3,
     $     vec6,vec7,vec8,vec9,vec10,irow,numrow,numrowdim)

c error check
      
      ierrorflag = 0
      
      if (npmax.gt.(namax+1)/2)
     $     call error('Ftmp possibly too small - boost namax')
      if (rs*sigma.ge.xprd2.or.rs*sigma.ge.yprd2.or.rs*sigma.ge.zprd2)
     $     call error('Outer cutoff >= 1/2 box size')
      if (ineigh.eq.1.and.nbinx*nbiny*nbinz.gt.nbmax)
     $     call error('Not enough bins - boost nbmax')
      if (ineigh.eq.1.and.(nbinx.le.2.or.nbiny.le.2.or.nbinz.le.2))
     $     call error('Two or less bins in a dimension')
      if (ntimes/nstat+1.gt.ntmax)
     $     call error('Stats array too small - boost ntmax')
      
      iflag = 0
      if (nglocal.gt.namax.or.ngolocal.gt.namax) iflag = -1
      call merge_i(iflag,node,nprocs,ncube)
      if (iflag.lt.0)
     $     call error('Too many row/col atoms - boost namax')
      iflag = 0
      if (nlocal.gt.npmax) iflag = -1
      call merge_i(iflag,node,nprocs,ncube)
      if (iflag.lt.0)
     $     call error('Too many atoms/processor - boost npmax')

      if (ierrorflag.lt.0) return
      
c setup ngbefore and pntcopy
      
      call split_scalar(natoms,i,ngbefore,j,irow*numcol,nprocs)
      ngbefore = ngbefore - 1
      
      do i = nstart,nend
        x(1,i) = ngbefore + i
      enddo
      call swap(node,x(1,nstart),3*nbyte*nlocal,itransend,
     $     xcopy(1,nostart),3*nbyte*nolocal,itranrecv,0,0,0,1)
      call expand(xcopy,vec6,vec7,vec8,vec9,vec10,
     $     irow,numrow,numrowdim)
      do i = 1,ngolocal
        pntcopy(i) = xcopy(1,i)
      enddo

      return
      end
      
      
c initialize atoms on fcc lattice
c  initialize velocities so that each node generates all velocities, only
c   saves its own
c  this is slow, but nice for debugging

      subroutine setup_atom
      include 'ljf.h'
      
c create fcc lattice
c  all nodes generate entire lattice, only save theirs
c  generate atom positions sequentially, use RNG to permute order
c  works because RNG is full period (im) by Knuth p. 16,20
c  throw out numbers for ivalue >= natoms
c  uncomment "n = m" lines (2 places) to use standard ordering
      
c im = power of 2 >= natoms
c ia,ic (prime) = RNG constants, see Knuth
c nlo/nhi = global bounds of this node's atoms
      
      im = 1
 10   im = im * 2
      if (im.lt.natoms) goto 10
      ia = 365
      ic = 8161
      
      nlo = ngbefore + nstart
      nhi = ngbefore + nend
      
      ivalue = 0
      m = 0
      do k = 1,nzsize*2
        do j = 1,nysize*2
          do i = 1,nxsize*2
            if (mod(i+j+k,2).eq.1) then
              m = m + 1
 20           ivalue = mod(ia*ivalue+ic,im)
              if (ivalue.ge.natoms) goto 20
              n = ivalue + 1
c              n = m
              if (n.ge.nlo.and.n.le.nhi) then
                ntmp = n - ngbefore
                x(1,ntmp) = (i-1)*alat/2.0 - xprd2
                x(2,ntmp) = (j-1)*alat/2.0 - yprd2
                x(3,ntmp) = (k-1)*alat/2.0 - zprd2
              endif
            endif
          enddo
        enddo
      enddo
      
      call expand(x,vec1,vec2,vec3,vec4,vec5,
     $     icol,numcol,numcoldim)
      call swap(node,x(1,nstart),3*nbyte*nlocal,itransend,
     $     xcopy(1,nostart),3*nbyte*nolocal,itranrecv,0,0,0,1)
      call expand(xcopy,vec6,vec7,vec8,vec9,vec10,
     $     irow,numrow,numrowdim)
      
c initialize (in permuted order) and scale velocities
      
      ivalue = 0
      vxtot = 0.0
      vytot = 0.0
      vztot = 0.0
      do m = 1,natoms
        vxtmp = random(iseed)*dt
        vytmp = random(iseed)*dt
        vztmp = random(iseed)*dt
 30     ivalue = mod(ia*ivalue+ic,im)
        if (ivalue.ge.natoms) goto 30
        n = ivalue + 1
c        n = m
        if (n.ge.nlo.and.n.le.nhi) then
          ntmp = n - nlo + 1
          v(1,ntmp) = vxtmp
          v(2,ntmp) = vytmp
          v(3,ntmp) = vztmp
        endif
        vxtot = vxtot + vxtmp
        vytot = vytot + vytmp
        vztot = vztot + vztmp
      enddo

      do i = 1,nlocal
        v(1,i) = v(1,i) - vxtot/natoms
        v(2,i) = v(2,i) - vytot/natoms
        v(3,i) = v(3,i) - vztot/natoms
        vold(1,i) = v(1,i)
        vold(2,i) = v(2,i)
        vold(3,i) = v(3,i)
      enddo

      call temperature(t)
      factor = sqrt(tstar/t)
      do i = 1,nlocal
        v(1,i) = v(1,i) * factor
        v(2,i) = v(2,i) * factor
        v(3,i) = v(3,i) * factor
        vold(1,i) = v(1,i)
        vold(2,i) = v(2,i)
        vold(3,i) = v(3,i)
      enddo
      
      return
      end


c leapfrog integrator

      subroutine integrate(m)
      include 'ljf.h'
      real*8 time1,time2

      j = 0
      do i = nstart,nend
        j = j + 1
        x(1,i) = x(1,i) + dt*v(1,j)
        if (x(1,i).ge.xprd2) x(1,i) = x(1,i) - xprd
        if (x(1,i).lt.-xprd2) x(1,i) = x(1,i) + xprd
        x(2,i) = x(2,i) + dt*v(2,j)
        if (x(2,i).ge.yprd2) x(2,i) = x(2,i) - yprd
        if (x(2,i).lt.-yprd2) x(2,i) = x(2,i) + yprd
        x(3,i) = x(3,i) + dt*v(3,j)
        if (x(3,i).ge.zprd2) x(3,i) = x(3,i) - zprd
        if (x(3,i).lt.-zprd2) x(3,i) = x(3,i) + zprd
      enddo
        
      call timer(time1)
      call expand(x,vec1,vec2,vec3,vec4,vec5,
     $     icol,numcol,numcoldim)
      call swap(node,x(1,nstart),3*nbyte*nlocal,itransend,
     $     xcopy(1,nostart),3*nbyte*nolocal,itranrecv,0,0,0,1)
      call expand(xcopy,vec6,vec7,vec8,vec9,vec10,
     $     irow,numrow,numrowdim)
      call timer(time2)
      time_comm = time_comm + time2-time1

      if (mod(m,nneigh).eq.0) then
        call timer(time1)
        call neighbor
        call timer(time2)
        time_neigh = time_neigh + time2-time1
      endif

      call timer(time1)
      if (inewt.eq.0) then
        call force0
      else
        call force1
      endif
      call timer(time2)
      time_force = time_force + time2-time1

      call timer(time1)
      if (inewt.ne.0) then
        call fold(fcopy,ftmp,vec6,vec7,vec8,vec9,vec10,
     $       irow,numrow,numrowdim)
        call swap(node,fcopy(1,nostart),3*nbyte*nolocal,itranrecv,
     $     ftmp,3*nbyte*nlocal,itransend,0,0,0,1)
        j = 1
        do i = nstart,nend
          f(1,i) = f(1,i) + ftmp(1,j)
          f(2,i) = f(2,i) + ftmp(2,j)
          f(3,i) = f(3,i) + ftmp(3,j)
          j = j + 1
        enddo
      
      endif
      call fold(f,ftmp,vec1,vec2,vec3,vec4,vec5,
     $     icol,numcol,numcoldim)
      call timer(time2)
      time_comm = time_comm + time2-time1

      j = 0
      do i = nstart,nend
        j = j + 1
        vold(1,j) = v(1,j)
        vold(2,j) = v(2,j)
        vold(3,j) = v(3,j)
        v(1,j) = v(1,j) + dtforce*f(1,i)
        v(2,j) = v(2,j) + dtforce*f(2,i)
        v(3,j) = v(3,j) + dtforce*f(3,i)
      enddo

      return
      end
      
      
c driver for neighbor-list creation

      subroutine neighbor
      include 'ljf.h'
      
      if (inewt.eq.0) then
        if (ineigh.eq.0) then
          call neighbor0
        else
          call neighbor1
        endif
      else
        if (ineigh.eq.0) then
          call neighbor2
        else
          call neighbor3
        endif
      endif
      
      mneigh = mneigh + 1
      neighmax = max(neighmax,nnlist(nglocal+1)-1)
      
      return
      end
      
      
c no Newton's 3rd law, no binning
c  each group atom looks at all atoms in other group
      
      subroutine neighbor0
      include 'ljf.h'
      
      ii = ngbefore
      npnt = 1
      do i = 1,nglocal
        ii = ii + 1
        nnlist(i) = npnt
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do j = 1,ngolocal
          jj = pntcopy(j)
          delx = xtmp - xcopy(1,j)
          dely = ytmp - xcopy(2,j)
          delz = ztmp - xcopy(3,j)
          if (abs(delx).gt.xms) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.yms) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zms) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          if (rsq.le.cutsq2.and.ii.ne.jj) then
            nlist(npnt) = j
            npnt = npnt + 1
          endif
        enddo
        if (npnt.gt.npmax*nnmax) then
          write (6,*) 'Neighbor list too big:',node,npnt,i,nglocal
          call exit(0)
        endif
      enddo
      nnlist(nglocal+1) = npnt
      
      return
      end
      
      
c no Newton's 3rd law, binning
c  all other group atoms binned, each group atom looks at other groups bins
      
      subroutine neighbor1
      include 'ljf.h'
      
      do i = 1,nbinx*nbiny*nbinz
        binpnt(i) = 0
      enddo
      
      do i = 1,ngolocal
        ix = (xcopy(1,i) + xprd2) / binsizex
        iy = (xcopy(2,i) + yprd2) / binsizey
        iz = (xcopy(3,i) + zprd2) / binsizez
        ib = iz*nbiny*nbinx + iy*nbinx + ix + 1
        bin(i) = binpnt(ib)
        binpnt(ib) = i
      enddo
      
      ii = ngbefore
      npnt = 1
      do i = 1,nglocal
        ii = ii + 1
        nnlist(i) = npnt
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        ixx = (xtmp + xprd2) / binsizex
        iyy = (ytmp + yprd2) / binsizey
        izz = (ztmp + zprd2) / binsizez
        do k = 0,26
          ix = ixx + mod(k,3) - 1
          if (ix.lt.0) ix = nbinx - 1
          if (ix.eq.nbinx) ix = 0
          iy = iyy + mod(k/3,3) - 1
          if (iy.lt.0) iy = nbiny - 1
          if (iy.eq.nbiny) iy = 0
          iz = izz + k/9 - 1
          if (iz.lt.0) iz = nbinz - 1
          if (iz.eq.nbinz) iz = 0
          ib = iz*nbiny*nbinx + iy*nbinx + ix + 1
          j = binpnt(ib)
 30       if (j.ne.0) then
            jj = pntcopy(j)
            delx = xtmp - xcopy(1,j)
            dely = ytmp - xcopy(2,j)
            delz = ztmp - xcopy(3,j)
            if (abs(delx).gt.xms) then
              if (delx.lt.0.0) then
                delx = delx + xprd
              else
                delx = delx - xprd
              endif
            endif
            if (abs(dely).gt.yms) then
              if (dely.lt.0.0) then
                dely = dely + yprd
              else
                dely = dely - yprd
              endif
            endif
            if (abs(delz).gt.zms) then
              if (delz.lt.0.0) then
                delz = delz + zprd
              else
                delz = delz - zprd
              endif
            endif
            rsq = delx*delx + dely*dely + delz*delz
            if (rsq.le.cutsq2.and.ii.ne.jj) then
              nlist(npnt) = j
              npnt = npnt + 1
            endif
            j = bin(j)
            goto 30
          endif
        enddo
        if (npnt.gt.npmax*nnmax) then
          write (6,*) 'Neighbor list too big:',node,npnt,i,nglocal
          call exit(0)
        endif
      enddo
      nnlist(nglocal+1) = npnt
      
      return
      end
      
      
c Newton's 3rd law, no binning
c  each group atom looks at all atoms in other group, stores 1/2
      
      subroutine neighbor2
      include 'ljf.h'
      
      ii = ngbefore
      npnt = 1
      do i = 1,nglocal
        ii = ii + 1
        nnlist(i) = npnt
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do j = 1,ngolocal
          jj = pntcopy(j)
          if (ii.le.jj) then
            if (mod(ii+jj,2).eq.0) goto 10
          else
            if (mod(ii+jj,2).eq.1) goto 10
          endif
          delx = xtmp - xcopy(1,j)
          dely = ytmp - xcopy(2,j)
          delz = ztmp - xcopy(3,j)
          if (abs(delx).gt.xms) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.yms) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zms) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          if (rsq.le.cutsq2) then
            nlist(npnt) = j
            npnt = npnt + 1
          endif
 10       continue
        enddo
        if (npnt.gt.npmax*nnmax) then
          write (6,*) 'Neighbor list too big:',node,npnt,i,nglocal
          call exit(0)
        endif
      enddo
      nnlist(nglocal+1) = npnt
      
      return
      end
      
c Newton's 3rd law, binning
c  all other group atoms binned, each group atom looks at other groups bins,
c  stores 1/2
      
      subroutine neighbor3
      include 'ljf.h'
      
      do i = 1,nbinx*nbiny*nbinz
        binpnt(i) = 0
      enddo
      
      do i = 1,ngolocal
        ix = (xcopy(1,i) + xprd2) / binsizex
        iy = (xcopy(2,i) + yprd2) / binsizey
        iz = (xcopy(3,i) + zprd2) / binsizez
        ib = iz*nbiny*nbinx + iy*nbinx + ix + 1
        bin(i) = binpnt(ib)
        binpnt(ib) = i
      enddo
      
      ii = ngbefore
      npnt = 1
      do i = 1,nglocal
        ii = ii + 1
        nnlist(i) = npnt
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        ixx = (xtmp + xprd2) / binsizex
        iyy = (ytmp + yprd2) / binsizey
        izz = (ztmp + zprd2) / binsizez
        do k = 0,26
          ix = ixx + mod(k,3) - 1
          if (ix.lt.0) ix = nbinx - 1
          if (ix.eq.nbinx) ix = 0
          iy = iyy + mod(k/3,3) - 1
          if (iy.lt.0) iy = nbiny - 1
          if (iy.eq.nbiny) iy = 0
          iz = izz + k/9 - 1
          if (iz.lt.0) iz = nbinz - 1
          if (iz.eq.nbinz) iz = 0
          ib = iz*nbiny*nbinx + iy*nbinx + ix + 1
          j = binpnt(ib)
 30       if (j.ne.0) then
            jj = pntcopy(j)
            if (ii.le.jj) then
              if (mod(ii+jj,2).eq.0) goto 35
            else
              if (mod(ii+jj,2).eq.1) goto 35
            endif
            delx = xtmp - xcopy(1,j)
            dely = ytmp - xcopy(2,j)
            delz = ztmp - xcopy(3,j)
            if (abs(delx).gt.xms) then
              if (delx.lt.0.0) then
                delx = delx + xprd
              else
                delx = delx - xprd
              endif
            endif
            if (abs(dely).gt.yms) then
              if (dely.lt.0.0) then
                dely = dely + yprd
              else
                dely = dely - yprd
              endif
            endif
            if (abs(delz).gt.zms) then
              if (delz.lt.0.0) then
                delz = delz + zprd
              else
                delz = delz - zprd
              endif
            endif
            rsq = delx*delx + dely*dely + delz*delz
            if (rsq.le.cutsq2) then
              nlist(npnt) = j
              npnt = npnt + 1
            endif
 35         j = bin(j)
            goto 30
          endif
        enddo
        if (npnt.gt.npmax*nnmax) then
          write (6,*) 'Neighbor list too big:',node,npnt,i,nglocal
          call exit(0)
        endif
      enddo
      nnlist(nglocal+1) = npnt
      
      return
      end
      
      
c no Newton's 3rd law
c  store computed force only for atoms in group
      
      subroutine force0
      include 'ljf.h'
      
      do i = 1,nglocal
        f(1,i) = 0.0
        f(2,i) = 0.0
        f(3,i) = 0.0
      enddo
      
      do i = 1,nglocal
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnlist(i),nnlist(i+1)-1
          j = nlist(k)
          delx = xtmp - xcopy(1,j)
          dely = ytmp - xcopy(2,j)
          delz = ztmp - xcopy(3,j)
          if (abs(delx).gt.xmc) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          if (rsq.lt.cutsq1) then
            sr2 = sigsq/rsq
            sr6 = sr2*sr2*sr2
            tmp = sr6*(sr6-0.5)/rsq
            f(1,i) = f(1,i) + delx*tmp
            f(2,i) = f(2,i) + dely*tmp
            f(3,i) = f(3,i) + delz*tmp
          endif
        enddo
      enddo
      
      return
      end
      
      
c Newton's 3rd law
c  store computed force for atoms in own group and other group
      
      subroutine force1
      include 'ljf.h'
      
      do i = 1,nglocal
        f(1,i) = 0.0
        f(2,i) = 0.0
        f(3,i) = 0.0
      enddo
      
      do i = 1,ngolocal
        fcopy(1,i) = 0.0
        fcopy(2,i) = 0.0
        fcopy(3,i) = 0.0
      enddo
      
      do i = 1,nglocal
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnlist(i),nnlist(i+1)-1
          j = nlist(k)
          delx = xtmp - xcopy(1,j)
          dely = ytmp - xcopy(2,j)
          delz = ztmp - xcopy(3,j)
          if (abs(delx).gt.xmc) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          if (rsq.lt.cutsq1) then
            sr2 = sigsq/rsq
            sr6 = sr2*sr2*sr2
            tmp = sr6*(sr6-0.5)/rsq
            f(1,i) = f(1,i) + delx*tmp
            f(2,i) = f(2,i) + dely*tmp
            f(3,i) = f(3,i) + delz*tmp
            fcopy(1,j) = fcopy(1,j) - delx*tmp
            fcopy(2,j) = fcopy(2,j) - dely*tmp
            fcopy(3,j) = fcopy(3,j) - delz*tmp
          endif
        enddo
      enddo
      
      return
      end


c thermodynamic computations

      subroutine status
      include 'ljf.h'
      
      mstat = mstat + 1
      call temperature(t)
      call energy(eng)
      call pressure(p,t)
      tmparr(mstat) = t
      engarr(mstat) = eng
      prsarr(mstat) = p
      if (mstat.eq.1) then
        enginit = 1.5*t + eng
        conarr(mstat) = 1.0
      else
        conarr(mstat) = (1.5*t+eng)/enginit
      endif

      return
      end
      
      
c reduced temperature

      subroutine temperature(t)
      include 'ljf.h'
      
      t = 0.0
      do i = 1,nlocal
        vx = (v(1,i)+vold(1,i))/2.0
        vy = (v(2,i)+vold(2,i))/2.0
        vz = (v(3,i)+vold(3,i))/2.0
        t = t + vx*vx + vy*vy + vz*vz
      enddo
      call merge_r(t,node,nprocs,ncube)
      t = t / (3.0*natoms)
      
      return
      end
      
      
c reduced potential energy

      subroutine energy(eng)
      include 'ljf.h'
      
      eng = 0.0
      do i = 1,nglocal
        do j = nnlist(i),nnlist(i+1)-1
          rsq = rdstsq(i,nlist(j))
          if (rsq.lt.cutsq1) eng = eng + phi(rsq)
        enddo
      enddo
      call merge_r(eng,node,nprocs,ncube)
      if (inewt.eq.0) eng = eng/2.0
      eng = eng/natoms
      
      return
      end
      
      
c reduced pressure from virial

      subroutine pressure(p,t)
      include 'ljf.h'
      
      virial = 0.0
      do i = 1,nglocal
        do j = nnlist(i),nnlist(i+1)-1
          rsq = rdstsq(i,nlist(j))
          if (rsq.lt.cutsq1) virial = virial + fphi(rsq)
        enddo
      enddo
      call merge_r(virial,node,nprocs,ncube)
      if (inewt.eq.0) virial = virial/2.0
      p = t*rhostar + rhostar/3.0/natoms*virial
      
      return
      end
      
      
c distance between two atoms

      real*4 function rdstsq(i,j)
      include 'ljf.h'
      
      delx = x(1,i) - xcopy(1,j)
      dely = x(2,i) - xcopy(2,j)
      delz = x(3,i) - xcopy(3,j)
      if (abs(delx).gt.xmc) then
        if (delx.lt.0.0) then
          delx = delx + xprd
        else
          delx = delx - xprd
        endif
      endif
      if (abs(dely).gt.ymc) then
        if (dely.lt.0.0) then
          dely = dely + yprd
        else
          dely = dely - yprd
        endif
      endif
      if (abs(delz).gt.zmc) then
        if (delz.lt.0.0) then
          delz = delz + zprd
        else
          delz = delz - zprd
        endif
      endif
      rdstsq = delx*delx + dely*dely + delz*delz

      return
      end
      
      
c output to screen and file
c  add long-range correction to energy and pressure

      subroutine output
      include 'ljf.h'
      integer ihisto(10),ihistotmp(10)
      real*8 rtmp,ttmp,ave,xmax,xmin
 900  format(' ',a,f15.6,f13.4)
 901  format(' ',a,f13.4,a,f13.4,a,f13.4,a)
 902  format(' ',a,10i5)

      engcorr = 8.0*3.1415926*rhostar *
     $     (1.0/(9.0*rc**9) - 1.0/(3.0*rc**3))
      prscorr = 8.0*3.1415926*rhostar*rhostar *
     $     (4.0/(9.0*rc**9) - 2.0/(3.0*rc**3))
      
      if (node.eq.0) then
        
        open (unit=1,file='lj.out')
        close (1,status='delete')
        open (unit=1,file='lj.out')
        m = ntimes/nstat + 1
        tmpave = -tmparr(1)
        engave = -engarr(1) - engcorr
        prsave = -prsarr(1) - prscorr
        conave = -conarr(1)
        write (1,*) 'Timestep, T*, U*, P*, Conservation:'
        do i = 1,m
          engarr(i) = engarr(i) + engcorr
          prsarr(i) = prsarr(i) + prscorr
          write (1,*) (i-1)*nstat,
     $         tmparr(i),engarr(i),prsarr(i),conarr(i)
          tmpave = tmpave + tmparr(i)
          engave = engave + engarr(i)
          prsave = prsave + prsarr(i)
          conave = conave + conarr(i)
        enddo
        if (m.gt.1) then
          write (1,*) 'Averages: (without timestep 0)' 
          write (1,*) '     ',tmpave/(m-1),engave/(m-1),
     $         prsave/(m-1),conave/(m-1)
        endif

      endif
      
      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      rtmp = time_total
      call merge_r8(rtmp,node,nprocs,ncube)
      ttmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,*) 'Total time:',ttmp,
     $       ' on',nprocs,' procs for',natoms,' atoms'
        write (1,*) 'Total time:',ttmp,
     $       ' on',nprocs,' procs for',natoms,' atoms'
      endif

      if (ttmp.eq.0.0) ttmp = 1.0

      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      rtmp = time_force
      call merge_r8(rtmp,node,nprocs,ncube)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Force time/%:',rtmp,rtmp/ttmp*100
        write (1,900) 'Force time/%:',rtmp,rtmp/ttmp*100
      endif
      
      rtmp = time_neigh
      call merge_r8(rtmp,node,nprocs,ncube)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Neigh time/%:',rtmp,rtmp/ttmp*100
        write (1,900) 'Neigh time/%:',rtmp,rtmp/ttmp*100
      endif
      
      rtmp = time_comm
      call merge_r8(rtmp,node,nprocs,ncube)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Comm  time/%:',rtmp,rtmp/ttmp*100
        write (1,900) 'Comm  time/%:',rtmp,rtmp/ttmp*100
      endif
      
      rtmp = time_total - (time_force + time_neigh + time_comm)
      call merge_r8(rtmp,node,nprocs,ncube)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Other time/%:',rtmp,rtmp/ttmp*100
        write (1,900) 'Other time/%:',rtmp,rtmp/ttmp*100
      endif
      
      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      call stats_r8(time_force,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Force time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Force time:',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      call stats_r8(time_neigh,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Neigh time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Neigh time:',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      call stats_r8(time_comm,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Comm  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Comm  time:',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      rtmp = time_total - (time_force + time_neigh + time_comm)
      call stats_r8(rtmp,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Other time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Other time:',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      rtmp = nnlist(nglocal+1) - 1
      call stats_r8(rtmp,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      call merge_imax(neighmax,node,nprocs,ncube)
      if (node.eq.0) then
        write (1,*) 'Max # of neighbors =',neighmax,
     $       ' out of',nnmax*npmax
      endif
      
      if (node.eq.0) close (1)

      return
      end
      

c derivative of LJ energy (force)

      real*4 function fphi(rsq)
      include 'ljf.h'
      
      sr2 = sigsq/rsq
      sr6 = sr2*sr2*sr2
      fphi = 48.0*sr6*(sr6-0.5)

      return
      end
      
      
c LJ energy

      real*4 function phi(rsq)
      include 'ljf.h'
      
      sr2 = sigsq/rsq
      sr6 = sr2*sr2*sr2
      sr12 = sr6*sr6
      phi = 4.0*(sr12-sr6)

      return
      end
      
      
c Park/Miller RNG

      real*4 function random(iseed)
      real*8 aa,mm,sseed
      parameter (aa=16807.0D0,mm=2147483647.0D0)
      
      sseed = iseed
      sseed = mod(aa*sseed,mm)
      random = sseed/mm
      iseed = sseed

      return
      end


c set error flag, print message

      subroutine error(str)
      include 'ljf.h'
      character*(*) str

      ierrorflag = -1
      if (node.eq.0) write (6,*) 'Error: ',str

      return
      end
