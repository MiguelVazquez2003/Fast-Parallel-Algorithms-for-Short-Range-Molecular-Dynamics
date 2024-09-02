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

c nCUBE 1,000,000 run
      parameter (namax=32000,npmax=2000)
      parameter (nnmax=80,nbmax=54000,ntmax=101)

c iPSC/860 500,000 run
c      parameter (namax=64000,npmax=8000)
c      parameter (nnmax=60,nbmax=27000,ntmax=101)

c Delta 1,000,000 run
c      parameter (namax=64000,npmax=4000)
c      parameter (nnmax=80,nbmax=54000,ntmax=101)

c Paragon 1,000,000 run on 128 or more processors
c      parameter (namax=128000,npmax=8000)
c      parameter (nnmax=80,nbmax=54000,ntmax=101)

      real*4 x(3,namax),xcopy(3,namax)
      real*4 v(3,npmax),vold(3,npmax)
      real*4 f(3,namax),fcopy(3,namax),ftmp(3,(namax+1)/2)
      integer pntcopy(namax)

      integer nlist(npmax*nnmax+nnmax),nnlist(namax+1)
      integer bin(namax),binpnt(nbmax)
      
      real*4 tmparr(ntmax),engarr(ntmax),prsarr(ntmax),conarr(ntmax)
      
      integer vec1(20),vec2(20),vec3(20),vec4(20),vec5(20)
      integer vec6(20),vec7(20),vec8(20),vec9(20),vec10(20)

      real*8 time_total,time_comm,time_neigh,time_force

      common /blk00/ dtstar,tstar,rhostar,rc,rs,
     $     nneigh,nstat,iseed,nxsize,nysize,nzsize,ntimes,
     $     ineigh,ibin,nbinx,nbiny,nbinz,inewt,itrans,
     $     igrid,numrow,numcol
      common /blk01/ alat,dt,dtforce,sigma,sigsq,cutsq1,cutsq2
      common /blk02/ xprd,yprd,zprd,xprd2,yprd2,zprd2
      common /blk03/ xmc,ymc,zmc,xms,yms,zms

      common /blk10/ x,xcopy,v,vold,f,fcopy,ftmp,pntcopy
      common /blk11/ natoms,nstart,nend,nlocal,nostart,noend,nolocal
      common /blk12/ nglocal,ngolocal,ngbefore,itransend,itranrecv
      common /blk13/ numrowdim,numcoldim,irow,icol

      common /blk20/ nlist,nnlist,bin,binpnt
      common /blk21/ binsizex,binsizey,binsizez

      common /blk30/ enginit,tmparr,engarr,prsarr,conarr
      common /blk31/ vec1,vec2,vec3,vec4,vec5
      common /blk32/ vec6,vec7,vec8,vec9,vec10
      common /blk33/ mstat,mneigh,neighmax
      common /blk34/ time_total,time_comm,time_neigh,time_force
      common /blk35/ node,nprocs,ncube,ierrorflag
