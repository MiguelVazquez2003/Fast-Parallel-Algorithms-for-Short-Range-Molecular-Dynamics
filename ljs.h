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

      implicit real*4 (a-h,o-z)
      parameter (nbyte=4,ibyte=4)

c nCUBE 10,000,000 run on 1024 procs
      parameter (npmax=11000,nomax=10000,namax=npmax+nomax)
      parameter (nnmax=50,nbmax=10000,ntmax=101)
      parameter (nemax=20000,nfmax=10000,nsmax=32)

c Delta 10,000,000 run or iPSC/860 1,000,000 run
c      parameter (npmax=22000,nomax=21000,namax=npmax+nomax)
c      parameter (nnmax=50,nbmax=10000,ntmax=101)
c      parameter (nemax=20000,nfmax=20000,nsmax=32)

c Paragon 10,000,000 run on 512 procs - run with yods (comm=1000000)
c      parameter (npmax=22000,nomax=21000,namax=npmax+nomax)
c      parameter (nnmax=50,nbmax=10000,ntmax=101)
c      parameter (nemax=20000,nfmax=20000,nsmax=32)

c Paragon 100,000,000 run on 1840 procs - run with default stack,heap,comm
c      parameter (npmax=58000,nomax=40000,namax=npmax+nomax)
c      parameter (nnmax=42,nbmax=6500,ntmax=101)
c      parameter (nemax=45000,nfmax=35000,nsmax=32)

      real*4 x(3,namax)
      real*4 v(3,npmax),vold(3,npmax)
      real*4 f(3,npmax)
      integer list(npmax),atompnt,freepnt

      integer nlist(npmax*nnmax+nnmax),nnlist(npmax+1)
      integer bin(namax),binpnt(nbmax)

      real*4 boundlo(nsmax),boundhi(nsmax),border(2,3),prd(3)
      integer slist(nemax),nslist(nsmax+1),spart(nsmax),rpart(nsmax)
      integer mpart(2,3),npdim(3),need(3),me(3)

      real*4 buf1(nfmax+6),buf2(2*nfmax)

      real*4 tmparr(ntmax),engarr(ntmax),prsarr(ntmax),conarr(ntmax)
      
      real*8 time_total,time_neigh,time_comm,time_force

      common /blk00/ dtstar,tstar,rhostar,rc,rs,
     $     nneigh,nstat,iseed,nxsize,nysize,nzsize,ntimes,
     $     ineigh,ibin,nbinx,nbiny,nbinz,igrid,npdim
      common /blk01/ alat,dt,dtforce,sigma,sigsq,cutsq1,cutsq2
      common /blk02/ xprd,yprd,zprd,xprd2,yprd2,zprd2
      common /blk03/ xmc,ymc,zmc,xms,yms,zms

      common /blk10/ x,v,vold,f
      common /blk11/ natoms,nlocal,nother
      common /blk12/ atompnt,freepnt,list

      common /blk20/ nlist,nnlist,bin,binpnt
      common /blk21/ binsizex,binsizey,binsizez
      common /blk22/ mbinx,mbiny,mbinz,mbinxlo,mbinylo,mbinzlo
      common /blk23/ border,boundlo,boundhi,prd
      common /blk24/ slist,nslist,spart,rpart,mpart,need,me,nswap

      common /blk30/ enginit,tmparr,engarr,prsarr,conarr
      common /blk31/ buf1,buf2
      common /blk32/ mstat,mneigh,nlocalmax,nothermax,neighmax
      common /blk33/ nslistmax,nexcmax,nswpmax
      common /blk34/ time_total,time_neigh,time_comm,time_force
      common /blk35/ node,nprocs,ncube,ierrorflag
