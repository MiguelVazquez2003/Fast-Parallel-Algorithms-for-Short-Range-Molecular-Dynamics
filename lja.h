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

c for no Newton's 3rd law, set nfmax=npmax, nsmax=1
c for Newton's 3rd law, set nfmax=namax, nsmax=namax

c nCUBE 100,000 run
      parameter (namax=100000,npmax=500,nfmax=npmax,nsmax=1)

c iPSC/860 or Delta 100,000 run on 32 or more processors
c      parameter (namax=100000,npmax=4000,nfmax=npmax,nsmax=1)

c Paragon SUNMOS 100,000 run on 32 or more processors (use comm = 2000000)
c      parameter (namax=100000,npmax=4000,nfmax=npmax,nsmax=1)

      parameter (nnmax=80,nbmax=5400,ntmax=101)

      real*4 x(3,namax)
      real*4 v(3,npmax),vold(3,npmax)
      real*4 f(3,nfmax),ftmp(3,(nsmax+1)/2)

      integer nlist(npmax*nnmax+nnmax),nnlist(npmax+1)
      integer bin(namax),binpnt(nbmax),offset(3,14)
      
      real*4 tmparr(ntmax),engarr(ntmax),prsarr(ntmax),conarr(ntmax)
      
      integer vec1(20),vec2(20),vec3(20),vec4(20),vec5(20)

      real*8 time_total,time_neigh,time_comm,time_force

      common /blk00/ dtstar,tstar,rhostar,rc,rs,
     $     nneigh,nstat,iseed,nxsize,nysize,nzsize,ntimes,
     $     ineigh,ibin,nbinx,nbiny,nbinz,inewt
      common /blk01/ alat,dt,dtforce,sigma,sigsq,cutsq1,cutsq2
      common /blk02/ xprd,yprd,zprd,xprd2,yprd2,zprd2
      common /blk03/ xmc,ymc,zmc,xms,yms,zms

      common /blk10/ x,v,vold,f,ftmp
      common /blk11/ natoms,nstart,nend,nlocal

      common /blk20/ nlist,nnlist,bin,binpnt
      common /blk21/ binsizex,binsizey,binsizez,offset

      common /blk30/ enginit,tmparr,engarr,prsarr,conarr
      common /blk31/ vec1,vec2,vec3,vec4,vec5
      common /blk32/ mstat,mneigh,neighmax
      common /blk33/ time_total,time_neigh,time_comm,time_force

      common /blk34/ node,nprocs,ncube,ierrorflag




