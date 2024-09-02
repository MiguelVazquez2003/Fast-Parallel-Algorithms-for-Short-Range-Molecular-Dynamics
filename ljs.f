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

c spatial-decomposition algorithm
      
      program ljs
      
      include 'ljs.h'
      real*8 time1,time2
      
 10   call whatami(node,nprocs,ncube)
      
      call input
      
      if (node.eq.0) write (6,*) 'Setting up ...'

      call setup_general

      call setup_parallel
      if (ierrorflag.lt.0) goto 20

      call setup_atom
      if (ierrorflag.lt.0) goto 20

      call exchange
      call borders
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
      include 'ljs.h'

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
        if (nprocs.eq.2**ncube) then
          write (6,*) 'Specify processor grid: (0) No (1) Yes'
          read (5,*) igrid
        else
          igrid = 1
        endif
        if (igrid.eq.1) then
          write (6,*) 'Processor grid: x,y,z'
          read (5,*) npdim(1),npdim(2),npdim(3)
        endif

      endif

      inum = 16*ibyte + 5*nbyte
      call broadcast(dtstar,inum,node,nprocs,ncube)

      return
      end
        

c setup program units

      subroutine setup_general
      include 'ljs.h'

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
      
      time_force = 0.0
      time_neigh = 0.0
      time_comm = 0.0

      mstat = 0
      mneigh = 0

      nlocalmax = 0
      nothermax = 0
      neighmax = 0
      nslistmax = 0
      nexcmax = 0
      nswpmax = 0

      return
      end
      
      
c setup neighbor bins, spatial-decomposition communication patterns
c  error check

      subroutine setup_parallel
      include 'ljs.h'
      parameter (small=1.0E-6)
      
      if (igrid.eq.0) then
        npdim(1) = 2 ** ((ncube+0)/3)
        npdim(2) = 2 ** ((ncube+1)/3)
        npdim(3) = 2 ** ((ncube+2)/3)
      endif

      if (npdim(1)*npdim(2)*npdim(3).ne.nprocs) then
        call error('Bad grid of processors')
        return
      endif

      call mesh_3d(npdim(1),npdim(2),npdim(3),node,me(1),me(2),me(3),
     $     mpart(1,1),mpart(2,1),mpart(1,2),mpart(2,2),
     $     mpart(1,3),mpart(2,3))

      border(1,1) = float(me(1))/npdim(1) * xprd - xprd2
      border(2,1) = float(me(1)+1)/npdim(1) * xprd - xprd2
      border(1,2) = float(me(2))/npdim(2) * yprd - yprd2
      border(2,2) = float(me(2)+1)/npdim(2) * yprd - yprd2
      border(1,3) = float(me(3))/npdim(3) * zprd - zprd2
      border(2,3) = float(me(3)+1)/npdim(3) * zprd - zprd2
      
      need(1) = (rs*sigma) / (xprd/npdim(1)) + 1
      need(2) = (rs*sigma) / (yprd/npdim(2)) + 1
      need(3) = (rs*sigma) / (zprd/npdim(3)) + 1

c don't exchange if only 1 box in a dimension

      if (npdim(1).eq.1) need(1) = 0
      if (npdim(2).eq.1) need(2) = 0
      if (npdim(3).eq.1) need(3) = 0

c don't exchange more than 1/2 way over (e.g. 3 boxes away when npdim = 5)

      if (2*need(1).gt.npdim(1)) need(1) = need(1) - 1
      if (2*need(2).gt.npdim(2)) need(2) = need(2) - 1
      if (2*need(3).gt.npdim(3)) need(3) = need(3) - 1
      
c setup 4 parameters for each exchange: (spart,rpart,boundlo,boundhi)
c  (1,2) nodes to swap with
c  (3,4) slab boundaries (in correct dimension) of atoms that will be sent
c 1st part of if is sending to the west (south,down)
c 2nd part of if is sending to the east (north,up)
c nbox = box (in this dimension) who originally owned the atoms 
c   I will be sending in this swap
c kk/2+1 vs npdim/2 test is to make sure a node doesn't get more than 1/2
c   of a box from both directions (else would get some atoms twice)

      prd(1) = xprd
      prd(2) = yprd
      prd(3) = zprd
      nswap = 0
      
      do k = 1,3
        do kk = 0,2*need(k)-1
          nswap = nswap + 1
          if (nswap.le.nsmax) then
            if (mod(kk,2).eq.0) then
              spart(nswap) = mpart(1,k)
              rpart(nswap) = mpart(2,k)
              nbox = me(k) + kk/2
              if (nbox.ge.npdim(k)) nbox = nbox - npdim(k)
              blo = float(nbox)/npdim(k) * prd(k) - prd(k)/2.0
              if (kk/2+1.lt.need(k)) then
                bhi = float(nbox+1)/npdim(k) * prd(k) - prd(k)/2.0
              else
                bhi = border(1,k) + rs*sigma
                if (bhi.ge.prd(k)/2.0) bhi = bhi - prd(k)
              endif
              if (kk/2+1.eq.npdim(k)/2) then
                btmp = float(nbox+1)/npdim(k) * prd(k) - prd(k)/2.0
                bmid = (blo+btmp) / 2.0
                bhi = min(bhi,bmid)
              endif
            else
              spart(nswap) = mpart(2,k)
              rpart(nswap) = mpart(1,k)
              nbox = me(k) - kk/2
              if (nbox.lt.0) nbox = nbox + npdim(k)
              bhi = float(nbox+1)/npdim(k) * prd(k) - prd(k)/2.0
              if (kk/2+1.lt.need(k)) then
                blo = float(nbox)/npdim(k) * prd(k) - prd(k)/2.0
              else
                blo = border(2,k) - rs*sigma
                if (blo.lt.-prd(k)/2.0) blo = blo + prd(k)
              endif
              if (kk/2+1.eq.npdim(k)/2) then
                btmp = float(nbox)/npdim(k) * prd(k) - prd(k)/2.0
                bmid = (btmp+bhi) / 2.0
                blo = max(blo,bmid)
              endif
            endif
            boundlo(nswap) = blo
            boundhi(nswap) = bhi
          endif
        enddo
      enddo
      
c setup neighbor binning parameters in box owned by each processor
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
        
        mbinxlo = int((border(1,1) + xprd2) / binsizex) - 1
        mbinylo = int((border(1,2) + yprd2) / binsizey) - 1
        mbinzlo = int((border(1,3) + zprd2) / binsizez) - 1
        
        ix = int((border(2,1) + xprd2) / binsizex) + 1
        iy = int((border(2,2) + yprd2) / binsizey) + 1
        iz = int((border(2,3) + zprd2) / binsizez) + 1
        if (ix.gt.nbinx) ix = nbinx
        if (iy.gt.nbiny) iy = nbiny
        if (iz.gt.nbinz) iz = nbinz
        mbinx = ix - mbinxlo + 1
        mbiny = iy - mbinylo + 1
        mbinz = iz - mbinzlo + 1
        mbinx = min(mbinx,nbinx)
        mbiny = min(mbiny,nbiny)
        mbinz = min(mbinz,nbinz)
        
        if (mbinxlo.lt.0) mbinxlo = nbinx - 1
        if (mbinylo.lt.0) mbinylo = nbiny - 1
        if (mbinzlo.lt.0) mbinzlo = nbinz - 1
        
      endif
      
      ierrorflag = 0
      if (rs*sigma.ge.xprd2.or.rs*sigma.ge.yprd2.or.rs*sigma.ge.zprd2)
     $     call error('Outer cutoff >= 1/2 box size')
      if (ineigh.eq.1.and.(nbinx.le.2.or.nbiny.le.2.or.nbinz.le.2))
     $     call error('Two or less bins in a dimension')
      if (ntimes/nstat+1.gt.ntmax) 
     $     call error('Stats array too small - boost ntmax')
      if (nswap.gt.nsmax)
     $     call error('Swap array too small - boost nsmax')

      iflag = 0
      if (ineigh.eq.1.and.mbinx*mbiny*mbinz.gt.nbmax) iflag = -1
      call merge_i(iflag,node,nprocs,ncube)
      if (iflag.lt.0)
     $     call error('Too many local bins - boost nbmax')
      
      return
      end
      
      
c initialize atoms on fcc lattice
c  initialize velocities so that each node generates all velocities, only
c   saves its own
c  this is slow, but nice for debugging

      subroutine setup_atom
      include 'ljs.h'
      
      freepnt = 1
      do i = 1,npmax
        list(i) = i + 1
      enddo
      list(npmax) = 0
      nlocal = 0
      atompnt = npmax + 1
      
c each processor generates positions and velocities of all atoms
c  only stores ones in its box
      
      vxtot = 0.0
      vytot = 0.0
      vztot = 0.0
      do k = 1,nzsize*2
        do j = 1,nysize*2
          do i = 1,nxsize*2
            if (mod(i+j+k,2).eq.1) then
              xtmp = (i-1)*alat/2.0 - xprd2
              ytmp = (j-1)*alat/2.0 - yprd2
              ztmp = (k-1)*alat/2.0 - zprd2
              vxtmp = random(iseed)*dt
              vytmp = random(iseed)*dt
              vztmp = random(iseed)*dt
              vxtot = vxtot + vxtmp
              vytot = vytot + vytmp
              vztot = vztot + vztmp
              if (xtmp.ge.border(1,1).and.xtmp.lt.border(2,1).and.
     $             ytmp.ge.border(1,2).and.ytmp.lt.border(2,2).and.
     $             ztmp.ge.border(1,3).and.ztmp.lt.border(2,3)) then
                nlocal = nlocal + 1
                if (freepnt.ne.0) then
                  itmp = atompnt
                  atompnt = freepnt
                  freepnt = list(freepnt)
                  list(atompnt) = itmp
                  x(1,atompnt) = xtmp
                  x(2,atompnt) = ytmp
                  x(3,atompnt) = ztmp
                  v(1,atompnt) = vxtmp
                  v(2,atompnt) = vytmp
                  v(3,atompnt) = vztmp
                endif
              endif
            endif
          enddo
        enddo
        if (natoms.gt.1000000.and.mod(k,10).eq.0) then
          call synchro(node,nprocs,ncube)
          if (node.eq.0) write (6,*) 'Created z-layer',k/2
          call synchro(node,nprocs,ncube)
        endif
      enddo
      
      iflag = 0
      if (nlocal.gt.npmax) iflag = -1
      call merge_i(iflag,node,nprocs,ncube)
      if (iflag.lt.0)
     $     call error('Too many atoms/processor - boost npmax')
      if (ierrorflag.lt.0) return

      i = atompnt
      do ii = 1,nlocal
        v(1,i) = v(1,i) - vxtot/natoms
        v(2,i) = v(2,i) - vytot/natoms
        v(3,i) = v(3,i) - vztot/natoms
        vold(1,i) = v(1,i)
        vold(2,i) = v(2,i)
        vold(3,i) = v(3,i)
        i = list(i)
      enddo

      call temperature(t)
      factor = sqrt(tstar/t)
      i = atompnt
      do ii = 1,nlocal
        v(1,i) = v(1,i) * factor
        v(2,i) = v(2,i) * factor
        v(3,i) = v(3,i) * factor
        vold(1,i) = v(1,i)
        vold(2,i) = v(2,i)
        vold(3,i) = v(3,i)
        i = list(i)
      enddo
      
      return
      end
      
      
c leapfrog integrator

      subroutine integrate(m)
      include 'ljs.h'
      real*8 time1,time2
      
      i = atompnt
      do ii = 1,nlocal
        x(1,i) = x(1,i) + dt*v(1,i)
        if (x(1,i).ge.xprd2) x(1,i) = x(1,i) - xprd
        if (x(1,i).lt.-xprd2) x(1,i) = x(1,i) + xprd
        x(2,i) = x(2,i) + dt*v(2,i)
        if (x(2,i).ge.yprd2) x(2,i) = x(2,i) - yprd
        if (x(2,i).lt.-yprd2) x(2,i) = x(2,i) + yprd
        x(3,i) = x(3,i) + dt*v(3,i)
        if (x(3,i).ge.zprd2) x(3,i) = x(3,i) - zprd
        if (x(3,i).lt.-zprd2) x(3,i) = x(3,i) + zprd
        i = list(i)
      enddo
      
      if (mod(m,nneigh).ne.0) then
        call timer(time1)
        call communicate
        call timer(time2)
        time_comm = time_comm + time2-time1
      else
        call timer(time1)
        call exchange
        call borders
        call timer(time2)
        time_comm = time_comm + time2-time1
        call timer(time1)
        call neighbor
        call timer(time2)
        time_neigh = time_neigh + time2-time1
      endif
      
      call timer(time1)
      call force
      call timer(time2)
      time_force = time_force + time2-time1
      
      i = atompnt
      do ii = 1,nlocal
        vold(1,i) = v(1,i)
        vold(2,i) = v(2,i)
        vold(3,i) = v(3,i)
        v(1,i) = v(1,i) + dtforce*f(1,i)
        v(2,i) = v(2,i) + dtforce*f(2,i)
        v(3,i) = v(3,i) + dtforce*f(3,i)
        i = list(i)
      enddo
      
      return
      end
      
      
c swap slabs of atoms with other processors in all 3 directions
c  done every timestep
      
      subroutine communicate
      include 'ljs.h'
      
      ipnt = npmax
      do k = 1,nswap
        j = 0
        do ii = nslist(k),nslist(k+1)-1
          i = slist(ii)
          buf1(j+1) = x(1,i)
          buf1(j+2) = x(2,i)
          buf1(j+3) = x(3,i)
          j = j + 3
        enddo
        icnt = (namax-ipnt)*3*nbyte
        call swap(node,buf1,j*nbyte,spart(k),
     $       x(1,ipnt+1),icnt,rpart(k),1,0,0,0)
        ipnt = ipnt + icnt/3/nbyte
      enddo
      
      return
      end
      
      
c send out atoms that have left my box, receive ones entering my box
c  since last reneighboring
c  done in all 3 directions
c  3 positions and 3 velocities go with each atom
c  done before every reneighboring
      
      subroutine exchange
      include 'ljs.h'
      
      do k = 1,3
        
        if (npdim(k).gt.1) then
          
          blo = border(1,k)
          bhi = border(2,k)
          ndelete = 0
          iprev = 0
          j = 0
          i = atompnt
          
c fill buffer with atoms leaving my box, update local list
          
          do ii = 1,nlocal
            if (x(k,i).lt.blo.or.x(k,i).ge.bhi) then
              if (j.le.nfmax) then
                buf1(j+1) = x(1,i)
                buf1(j+2) = x(2,i)
                buf1(j+3) = x(3,i)
                buf1(j+4) = v(1,i)
                buf1(j+5) = v(2,i)
                buf1(j+6) = v(3,i)
              endif
              j = j + 6
              ndelete = ndelete + 1
              if (iprev.eq.0) then
                atompnt = list(i)
              else
                list(iprev) = list(i)
              endif
              itmp = list(i)
              list(i) = freepnt
              freepnt = i
              i = itmp
            else
              iprev = i
              i = list(i)
            endif
          enddo
          
          nlocal = nlocal - ndelete
          nexcmax = max(nexcmax,j/6)
          
          if (j.gt.nfmax) then
            write (6,*) 'Sending too many exchange atoms:',
     $           node,mpart(1,k),k
            call exit(0)
          endif
          
c send them out in both directions (if neighboring nodes are different)
          
          icnt = nfmax*nbyte
          call swap(node,buf1,j*nbyte,mpart(1,k),
     $         buf2,icnt,mpart(2,k),1,0,0,0)
          if (npdim(k).gt.2) then
            itmp = 2*nfmax*nbyte - icnt
            call swap(node,buf1,j*nbyte,mpart(2,k),
     $           buf2(icnt/nbyte+1),itmp,mpart(1,k),1,0,0,0)
            icnt = icnt + itmp
          endif
          
c check incoming atoms to see if they are in my box (could be in node's
c  box on other side of the sender)
          
          do j = 0,icnt/nbyte-1,6
            tmp = buf2(j+k)
            if (tmp.ge.blo.and.tmp.lt.bhi) then
              nlocal = nlocal + 1
              if (freepnt.ne.0) then
                itmp = atompnt
                atompnt = freepnt
                freepnt = list(freepnt)
                list(atompnt) = itmp
                x(1,atompnt) = buf2(j+1)
                x(2,atompnt) = buf2(j+2)
                x(3,atompnt) = buf2(j+3)
                v(1,atompnt) = buf2(j+4)
                v(2,atompnt) = buf2(j+5)
                v(3,atompnt) = buf2(j+6)
              endif
            endif
          enddo
          
          if (nlocal.gt.npmax) then
            write (6,*) 'Received too many exchange atoms:',
     $           node,nlocal,k
            call exit(0)
          endif
          
        endif
        
      enddo

      return
      end
      
      
c make lists of nearby atoms to send to neighboring nodes at every timestep
c  one list made for every swap that will be made
c  as list is made, actually do swaps
c  this does equivalent of a communicate (so don't need to explicitly
c   call communicate routine on reneighboring timestep)
c  3 atom positions are what is swapped
c  this routine called before every reneighboring
      
      subroutine borders
      include 'ljs.h'
      
      nswap = 0
      nother = 0
      npnt = 1
      
      do k = 1,3
        do kk = 0,2*need(k)-1
          
c buffer up all atoms I own (plus those previously received) that are
c  inside slab boundaries
c also store pointers to those atoms in slist for communicate routine
c  to use in future timesteps
          
          nswap = nswap + 1
          nslist(nswap) = npnt
          blo = boundlo(nswap)
          bhi = boundhi(nswap)
          
          j = 0
          i = atompnt
          do ii = 1,nlocal+nother
            if (x(k,i).ge.blo.and.x(k,i).lt.bhi) then
              if (npnt.le.nemax) slist(npnt) = i
              npnt = npnt + 1
              if (j.le.nfmax) then
                buf1(j+1) = x(1,i)
                buf1(j+2) = x(2,i)
                buf1(j+3) = x(3,i)
              endif
              j = j + 3
            endif
            if (ii.le.nlocal) then
              i = list(i)
            else
              i = i + 1
            endif
          enddo
          
          nswpmax = max(nswpmax,j/3)
          
          if (npnt.gt.nemax) then
            write (6,*) 'Too many atoms in border list:',
     $           node,npnt,k,kk
            call exit(0)
          endif
          
          if (j.gt.nfmax) then
            write (6,*) 'Sending too many border atoms:',
     $           node,j,k,kk
            call exit(0)
          endif
          
c swap atoms, put incoming ones at end of my position array
          
          icnt = (nomax-nother)*3*nbyte
          call swap(node,buf1,j*nbyte,spart(nswap),
     $         x(1,npmax+nother+1),icnt,rpart(nswap),1,0,0,0)
          nother = nother + icnt/3/nbyte
          
          if (nother.ge.nomax) then
            write (6,*) 'Received too many border atoms:',
     $           node,nother,k,kk
            call exit(0)
          endif
          
        enddo
      enddo
      
      nslist(nswap+1) = npnt
      
      return
      end


c driver for neighbor-list creation

      subroutine neighbor
      include 'ljs.h'

      if (ineigh.eq.0) then
        call neighbor0
      else
        call neighbor1
      endif
      
      mneigh = mneigh + 1
      nlocalmax = max(nlocalmax,nlocal)
      nothermax = max(nothermax,nother)
      neighmax = max(neighmax,nnlist(nlocal+1)-1)
      nslistmax = max(nslistmax,nslist(nswap+1)-1)

      return
      end
      
      
c no binning, Newton's 3rd law
c  N^2 / 2 search for neighbor pairs in my box
c  pair stored once in list if atoms i AND j are in my box (and i < j)
c  pair stored by me if j is NOT in my box (also stored by node owning j)
      
      subroutine neighbor0
      include 'ljs.h'
      
      npnt = 1
      i = atompnt
      do ii = 1,nlocal
        nnlist(ii) = npnt
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        j = list(i)
        do jj = ii+1,nlocal+nother
          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
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
          if (jj.le.nlocal) then
            j = list(j)
          else
            j = j + 1
          endif
        enddo
        if (npnt.gt.npmax*nnmax) then
          write (6,*) 'Neighbor list too big:',node,npnt,ii,nlocal
          call exit(0)
        endif
        i = list(i)
      enddo
      nnlist(nlocal+1) = npnt
      
      return
      end
      
      
c binning, Newton's 3rd law
c  all of mine and nearby atoms binned once
c  each owned atom i checks 27 surrounding boxes
c  pair stored once in list if atoms i AND j are in my box (and i < j)
c  pair stored by me if j is NOT in my box (also stored by node owning j)
      
      subroutine neighbor1
      include 'ljs.h'
      
      do i = 1,mbinx*mbiny*mbinz
        binpnt(i) = 0
      enddo
      
      i = atompnt
      do ii = 1,nlocal+nother
        ix = (x(1,i) + xprd2) / binsizex
        iy = (x(2,i) + yprd2) / binsizey
        iz = (x(3,i) + zprd2) / binsizez
        ix = ix - mbinxlo
        if (ix.lt.0) ix = ix + nbinx
        iy = iy - mbinylo
        if (iy.lt.0) iy = iy + nbiny
        iz = iz - mbinzlo
        if (iz.lt.0) iz = iz + nbinz
        if (ix.le.mbinx.and.iy.le.mbiny.and.iz.le.mbinz) then
          ib = iz*mbiny*mbinx + iy*mbinx + ix + 1
          bin(i) = binpnt(ib)
          binpnt(ib) = i
        endif
        if (ii.le.nlocal) then
          i = list(i)
        else
          i = i + 1
        endif
      enddo
      
      npnt = 1
      i = atompnt
      do ii = 1,nlocal
        nnlist(ii) = npnt
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        ixx = (xtmp + xprd2) / binsizex
        iyy = (ytmp + yprd2) / binsizey
        izz = (ztmp + zprd2) / binsizez
        do k = 0,26
          ix = ixx + mod(k,3) - 1
          iy = iyy + mod(k/3,3) - 1
          iz = izz + k/9 - 1
          ix = ix - mbinxlo
          if (ix.lt.0) ix = ix + nbinx
          if (ix.eq.nbinx) ix = 0
          iy = iy - mbinylo
          if (iy.lt.0) iy = iy + nbiny
          if (iy.eq.nbiny) iy = 0
          iz = iz - mbinzlo
          if (iz.lt.0) iz = iz + nbinz
          if (iz.eq.nbinz) iz = 0
          ib = iz*mbiny*mbinx + iy*mbinx + ix + 1
          j = binpnt(ib)
 30       if (j.ne.0) then
            if (j.le.i) goto 40
            delx = xtmp - x(1,j)
            dely = ytmp - x(2,j)
            delz = ztmp - x(3,j)
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
 40         j = bin(j)
            goto 30
          endif
        enddo
        if (npnt.gt.npmax*nnmax) then
          write (6,*) 'Neighbor list too big:',node,npnt,ii,nlocal
          call exit(0)
        endif
        i = list(i)
      enddo
      nnlist(nlocal+1) = npnt
      
      return
      end
      
      
c Newton's 3rd law
c  force components stored for atom i always and for atom j IF I own it
      
      subroutine force
      include 'ljs.h'
      
      i = atompnt
      do ii = 1,nlocal
        f(1,i) = 0.0
        f(2,i) = 0.0
        f(3,i) = 0.0
        i = list(i)
      enddo
      
      i = atompnt
      do ii = 1,nlocal
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnlist(ii),nnlist(ii+1)-1
          j = nlist(k)
          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
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
            if (j.le.npmax) then
              f(1,j) = f(1,j) - delx*tmp
              f(2,j) = f(2,j) - dely*tmp
              f(3,j) = f(3,j) - delz*tmp
            endif
          endif
        enddo
        i = list(i)
      enddo
      
      return
      end
      

c thermodynamic computations

      subroutine status
      include 'ljs.h'
      
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
      include 'ljs.h'
      
      t = 0.0
      i = atompnt
      do ii = 1,nlocal
        vx = (v(1,i)+vold(1,i))/2.0
        vy = (v(2,i)+vold(2,i))/2.0
        vz = (v(3,i)+vold(3,i))/2.0
        t = t + vx*vx + vy*vy + vz*vz
        i = list(i)
      enddo
      call merge_r(t,node,nprocs,ncube)
      t = t / (3.0*natoms)
      
      return
      end
      
      
c reduced potential energy

      subroutine energy(eng)
      include 'ljs.h'
      
      eng = 0.0
      i = atompnt
      do ii = 1,nlocal
        do j = nnlist(ii),nnlist(ii+1)-1
          rsq = rdstsq(i,nlist(j))
          if (rsq.lt.cutsq1) then
            factor = 1.0
            if (nlist(j).gt.npmax) factor = 0.5
            eng = eng + phi(rsq)*factor
          endif
        enddo
        i = list(i)
      enddo
      call merge_r(eng,node,nprocs,ncube)
      eng = eng/natoms
      
      return
      end
      
      
c reduced pressure from virial

      subroutine pressure(p,t)
      include 'ljs.h'
      
      virial = 0.0
      i = atompnt
      do ii = 1,nlocal
        do j = nnlist(ii),nnlist(ii+1)-1
          rsq = rdstsq(i,nlist(j))
          if (rsq.lt.cutsq1) then
            factor = 1.0
            if (nlist(j).gt.npmax) factor = 0.5
            virial = virial + fphi(rsq)*factor
          endif
        enddo
        i = list(i)
      enddo
      call merge_r(virial,node,nprocs,ncube)
      p = t*rhostar + rhostar/3.0/natoms*virial
      
      return
      end
      
      
c distance between two atoms

      real*4 function rdstsq(i,j)
      include 'ljs.h'
      
      delx = x(1,i) - x(1,j)
      dely = x(2,i) - x(2,j)
      delz = x(3,i) - x(3,j)
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
      include 'ljs.h'
      integer ihisto(10),ihistotmp(10)
      real*8 rtmp,ttmp,ave,xmax,xmin
 900  format(' ',a,f15.6,f13.4)
 901  format(' ',a,f13.4,a,f13.4,a,f13.4,a)
 902  format(' ',a,10i5)

c check for lost atoms

      nsum = nlocal
      call merge_i(nsum,node,nprocs,ncube)

      nlost = 0
      i = atompnt
      do ii = 1,nlocal
        if (x(1,i).lt.-xprd2.or.x(1,i).ge.xprd2.or.
     $       x(2,i).lt.-yprd2.or.x(2,i).ge.yprd2.or.
     $       x(3,i).lt.-zprd2.or.x(3,i).ge.zprd2)
     $       nlost = nlost + 1
        i = list(i)
      enddo
      call merge_i(nlost,node,nprocs,ncube)

      if (nsum.ne.natoms.or.nlost.gt.0) then
        if (node.eq.0) write (6,*) 'Atom counts =',nsum,nlost
        call error('Incorrect number of atoms')
      endif

c output

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

      rtmp = nnlist(nlocal+1) - 1
      call stats_r8(rtmp,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      rtmp = nlocal
      call stats_r8(rtmp,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Nlocal:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nlocal:    ',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      rtmp = nother
      call stats_r8(rtmp,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Nother:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nother:    ',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      rtmp = nslist(nswap+1) - 1
      call stats_r8(rtmp,1,ave,xmax,xmin,
     $     ihisto,ihistotmp,10,node,nprocs,ncube)
      if (node.eq.0) then
        write (6,901) 'Nswaps:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nswaps:    ',ave,' ave',xmax,' max',xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif
      
      call merge_imax(nlocalmax,node,nprocs,ncube)
      call merge_imax(nothermax,node,nprocs,ncube)
      call merge_imax(neighmax,node,nprocs,ncube)
      call merge_imax(nslistmax,node,nprocs,ncube)
      call merge_imax(nexcmax,node,nprocs,ncube)
      call merge_imax(nswpmax,node,nprocs,ncube)
      mbinmax = mbinx*mbiny*mbinz
      call merge_imax(mbinmax,node,nprocs,ncube)
      
      if (node.eq.0) then
        
        write (1,*) 'Max # of local atoms =',nlocalmax,
     $       ' out of',npmax
        write (1,*) 'Max # of other atoms =',nothermax,
     $       ' out of',nomax
        write (1,*) 'Max # of neighbors =',neighmax,
     $       ' out of',nnmax*npmax
        write (1,*) 'Max size of swap list =',nslistmax,
     $       ' out of',nemax
        write (1,*) 'Max used in exchange buffer =',nexcmax*6,
     $       ' out of',nfmax
        write (1,*) 'Max used in swap buffer =',nswpmax*3,
     $       ' out of',nfmax
        if (mbinmax.gt.0) write (1,*) 'Max # of bins =',mbinmax,
     $       ' out of',nbmax
        write (1,*)
        write (1,*) '# of swaps =',nswap,
     $       ' Needs =',need(1),need(2),need(3)
        write (1,*) 'Rs*sigma =',rs*sigma,' Cut/Box =',
     $       rs*sigma*npdim(1)/xprd,rs*sigma*npdim(2)/yprd,
     $       rs*sigma*npdim(3)/zprd
        
      endif

      if (node.eq.0) close (1)
      
      return
      end
      
      
c derivative of LJ energy (force)

      real*4 function fphi(rsq)
      include 'ljs.h'
      
      sr2 = sigsq/rsq
      sr6 = sr2*sr2*sr2
      fphi = 48.0*sr6*(sr6-0.5)

      return
      end
      
      
c LJ energy

      real*4 function phi(rsq)
      include 'ljs.h'
      
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


c set errorflag, print message

      subroutine error(str)
      include 'ljs.h'
      character*(*) str

      ierrorflag = -1
      if (node.eq.0) write (6,*) 'Error: ',str

      return
      end
