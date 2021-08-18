      subroutine msxfreq(xx0,mm,nclu,symb,el_zero,jobtype,icut,
     &                    es,emax,js,jmax,hso12,sc_qelec)

      implicit none

c     Passed into the subroutine
      double precision xx0(3,nclu),mm(nclu)
      integer nclu
      character*2 symb(nclu)
      character*5 jobtype
      logical icut

c     from input file
      double precision el_zero
      double precision es,emax
      integer js,jmax
      double precision hso12
      double precision sc_qelec

c     Phys constants and conversions
      double precision :: pi=dacos(-1.d0)
      double precision autoang,autokcal,autos
      double precision autoev,amutoau,autocmi
      double precision mu,kb
      parameter(autoang=0.52917706d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autos=2.4188843d-17)
      parameter(amutoau=1822.844987d0) 
      parameter(autoev=27.2113961d0)
      parameter(kb=3.166829d-6)     
      parameter(mu=1.d0*amutoau)

c     Grid and surf parameters
      integer mnsurf,mgrid
      parameter(mnsurf=3)
      parameter(mgrid=100000)
      integer nsurf
      integer mnclu

c     Variables Needed for msxfreq
      integer :: n1 = 1
      integer :: n2 = 2

c     Values initialized to Zero
      double precision :: ezero = 0.0
      double precision :: emin = 0.0
      integer :: jmin = 0
  
c     For linear checks      
      logical linear
      double precision linearinfinity

c     mass variables 
      double precision :: mtot = 0.0 
      double precision :: mx = 0.0 
      double precision :: my = 0.0 
      double precision :: mz = 0.0 

c     coordinates variables
      double precision xx(3,nclu)
      double precision x(nclu),y(nclu),z(nclu)
      double precision xxx
      double precision be

c     moment-of-inertia variables
      double precision rot(3,3),mom(3,3),momi(3),mom1,mom2
      double precision eig(3),ap(6)

c     gradients variables     
      double precision pemd(mnsurf,mnsurf)
      double precision gpemd(3,nclu,mnsurf,mnsurf)
      double precision gv1a(3,nclu),gv2a(3,nclu)
      double precision gv1b(3,nclu),gv2b(3,nclu)
      double precision gperp(3*nclu)
      double precision hh
      double precision g1g2normsum 
      double precision :: g1norm = 0.0
      double precision :: g2norm = 0.0

c     Hessian variables
      integer nmax,ndim
      double precision hess1(3*nclu,3*nclu),hess2(3*nclu,3*nclu)
      double precision mwhess1(3*nclu,3*nclu),mwhess2(3*nclu,3*nclu)
      double precision hessx(3*nclu,3*nclu)
      double precision hessval

c     Frequency variables
      integer nfreq
      double precision freq(3*nclu)

c     State-count variables
      double precision ee(50000)
      integer ne
      double precision fx(3*nclu)
      integer is,im
      double precision t(mgrid),at(mgrid),estep,ejk,djk,eee
      double precision t_lz(mgrid),t_ai(mgrid)
      integer imax,iejk,kmax,ifreq
      double precision etmp,ss
      double precision plz,rho,rholz,p2pass
      double precision a1,a2,a3,a4,e0,ex,tmpf,pairy,
     & airyarg,airypre
      double precision rholz_lz,rholz_ai,p2pass_lz,
     & p2pass_ai,rhox,rhox0
     
c     Variables for matrix diagonalization routines
      integer info
      integer lwork
      double precision work(9*nclu-1),work2(9)

c     Loop Variables
      integer i,j,ii,ij,jj,k,l,kl

c     Variables for Temp Storage
      double precision tmp,tmp1,tmp2,tmp3
      double precision temp1,temp2,temp3

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     Initialize certain variables as necessary

      lwork=9*nclu-1

      do i=1,3*nclu
        do j=1,3*nclu
          hess1(i,j) = 0.d0
          hess2(i,j) = 0.d0
         enddo
      enddo

      do i=1,nclu
        mtot=mtot+mm(i)
      enddo
      
      do i=1,3
        do j=1,nclu
          xx0(i,j)=xx0(i,j)/autoang
        enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate Moments-of-Inertia; Reorient Geometry

c     rotate geometry to align with principal axes of rotation
      write(6,*) "Rotating Geometry to Principal Axes..."
      write(6,*) 

c     center & reorient
      do i=1,nclu
        mx=mx+mm(i)*xx0(1,i)/mtot
        my=my+mm(i)*xx0(2,i)/mtot
        mz=mz+mm(i)*xx0(3,i)/mtot
      enddo
      do i=1,nclu
       xx0(1,i)=xx0(1,i)-mx
       xx0(2,i)=xx0(2,i)-my
       xx0(3,i)=xx0(3,i)-mz
      enddo

c     compute moment of intertia matrix (mom)
      do i=1,3
        do j=1,3
          mom(i,j) = 0.d0
        enddo
      enddo

      do i=1,nclu
         mom(1,1)=mom(1,1)+mm(i)*(xx0(2,i)**2+xx0(3,i)**2)
         mom(2,2)=mom(2,2)+mm(i)*(xx0(1,i)**2+xx0(3,i)**2)
         mom(3,3)=mom(3,3)+mm(i)*(xx0(1,i)**2+xx0(2,i)**2)
         mom(1,2)=mom(1,2)-mm(i)*(xx0(1,i)*xx0(2,i))
         mom(1,3)=mom(1,3)-mm(i)*(xx0(1,i)*xx0(3,i))
         mom(2,3)=mom(2,3)-mm(i)*(xx0(2,i)*xx0(3,i))
      enddo
      mom(2,1)=mom(1,2)
      mom(3,1)=mom(1,3)
      mom(3,2)=mom(2,3)

c     diagonalize the mom matrix
      do i=1,3
        do j=i,3
          ap(i+(j-1)*j/2)=mom(i,j)
        enddo
      enddo
      call dspev( 'v','u',3,ap,eig,rot,3,work2,info )

      do i=1,3
      momi(i)=0.5d0/eig(i)
      enddo

      write(6,'(a,3f15.5)')
     & " Moments of intertia (cm-1)",(autocmi*momi(i),i=1,3)
      tmp1=dabs(momi(1)-momi(2))
      tmp2=dabs(momi(2)-momi(3))
      tmp3=dabs(momi(1)-momi(3))
      if (tmp1.lt.tmp2.and.tmp1.lt.tmp3) then
      mom2=(momi(1)+momi(2))/2.d0
      mom1=momi(3)
      elseif (tmp2.lt.tmp3.and.tmp2.lt.tmp1) then
      mom2=(momi(2)+momi(3))/2.d0
      mom1=momi(1)
      else
      mom2=(momi(1)+momi(3))/2.d0
      mom1=momi(2)
      endif

c     check if the molecule is linear using mom eigenvalues 
      linearinfinity=1.d3
      if(mom1.gt.linearinfinity) then
      write(6,*)" Linear species found!"
      linear=.true.
      write(6,'(a,1f15.5,a)')" Symmetrized to (cm-1)     ",
     &    mom2*autocmi," (x2)"
      else
      linear=.false.
      write(6,'(a,2f15.5,a)')" Symmetrized to (cm-1)     ",
     &    mom1*autocmi,mom2*autocmi," (x2)"
      endif
      write(6,*)

c     rotate to diagonalize mom
      do i=1,nclu
         temp1 = xx0(1,i)
         temp2 = xx0(2,i)
         temp3 = xx0(3,i)
         xx0(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         xx0(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         xx0(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
      enddo

      do i=1,nclu
        x(i)=xx0(1,i)
        y(i)=xx0(2,i)
        z(i)=xx0(3,i)
      enddo

c     rotate geometry to align with principal axes of rotation
      write(6,*)"Rotated geometry (A)"
      do i=1,nclu
      write(6,'  (a,3f22.15)')symb(i),
     &      x(i)*autoang,y(i)*autoang,z(i)*autoang
      enddo
      write(6,*)
      write(6,*)

c     kill the calculation if only rotation is requested      
      if (jobtype.eq."ROTGM") then
        write(6,*) "Exiting Program..."
        stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the gradients
      write(6,*) "Calculating the Gradients..."
      write(6,*)

c     gradient of the gap
      do i=1,nclu
      x(i)=xx0(1,i)
      y(i)=xx0(2,i)
      z(i)=xx0(3,i)
      enddo
      nsurf=mnsurf
      mnclu=nclu
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
      pemd(n1,n1)=pemd(n1,n1)-el_zero
      pemd(n2,n2)=pemd(n2,n2)-el_zero
      write(6,*) "Energy of Rot'd MSX Geometry = ",
     &           pemd(n1,n1)*627.509,pemd(n2,n2)*627.509
      
      tmp1 = 0.0
      tmp2 = 0.0
      tmp3 = 0.0
      do i=1,3
        do j=1,nclu
          tmp1=tmp1+(gpemd(i,j,n1,n1)**2)*mu/mm(j)
          tmp2=tmp2+(gpemd(i,j,n2,n2)**2)*mu/mm(j)
          tmp3=tmp3+((gpemd(i,j,n1,n1)-gpemd(i,j,n2,n2))**2)*mu/mm(j)
        enddo
      enddo
      g1norm=dsqrt(tmp1)
      g2norm=dsqrt(tmp2)
      g1g2normsum=g1norm+g2norm
      tmp3=dsqrt(tmp3)
      do i=1,3
      do j=1,nclu
        ij = (i-1)*nclu + j
        gperp(ij)=(gpemd(i,j,n1,n1)-gpemd(i,j,n2,n2))
     &    *dsqrt(mu/mm(j))/tmp3
        write(6,'(3i5,3f15.5)')ij,i,j,gperp(ij)
      enddo
      enddo
      write(6,*) tmp3
      write(6,*) 

c     write gradient information to output file
      write(6, *) "|grad(n1)| / ( |grad(n1)| + |grad(n2)| ) =  ",
     &            g1norm / g1g2normsum
      write(6, *) "|grad(n2)| / ( |grad(n1)| + |grad(n2)| ) =  ",
     &            g2norm / g1g2normsum
      write(6,*) 
      write(6,*) 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       Obtain the Hessians

c     either calculates Hess with NST or reads them from hess.x file
c     procedure is determined based on user input
      write(6,*)"Obtaining the Hessians..."
      if (jobtype.eq."HESSC") then
        write(6,*)"Calculating Hessians for Each State with NST..."
        hh = 0.01
        do i=1,3
        do j=1,nclu
          write(6,'(4(a,i5))')" step (",i,",",j,") of ( 3,",nclu,")"
          ij = (i-1)*nclu + j
          do k=1,3
          do l=1,nclu
            xx(k,l) = xx0(k,l)
          enddo
          enddo
          xx(i,j) = xx0(i,j) + hh
        do k=1,nclu
        x(k)=xx(1,k)
        y(k)=xx(2,k)
        z(k)=xx(3,k)
        enddo
        call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
          do k=1,3
          do l=1,nclu
            gv1a(k,l) = gpemd(k,l,n1,n1)
            gv1b(k,l) = gpemd(k,l,n2,n2)
          enddo
          enddo
          xx(i,j) = xx0(i,j) - hh
        do k=1,nclu
        x(k)=xx(1,k)
        y(k)=xx(2,k)
        z(k)=xx(3,k)
        enddo
        call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
            do k=1,3
            do l=1,nclu
              gv2a(k,l) = gpemd(k,l,n1,n1)
              gv2b(k,l) = gpemd(k,l,n2,n2)
            enddo
            enddo
          do k=1,3
          do l=1,nclu
            kl = (k-1)*nclu + l
c           Hessian matrix, 3NCLU X 3NCLU matrix
c           data ordered (x1,x2,...,y1,...,z1,...,zNCLU)
            hess1(ij,kl) = (gv1a(k,l) - gv2a(k,l))/(2.d0*hh)
            hess2(ij,kl) = (gv1b(k,l) - gv2b(k,l))/(2.d0*hh)
c           mass-scale
            mwhess1(ij,kl) = hess1(ij,kl)*mu/dsqrt(mm(j)*mm(l))
            mwhess2(ij,kl) = hess2(ij,kl)*mu/dsqrt(mm(j)*mm(l))
          enddo
          enddo
        enddo
        enddo
      elseif (jobtype.eq."HESSR") then
        write(6,*)"Reading Hessians for Each State from Files..."
        open(30, file="hess.1")
        do i=1,3
          do j=1,nclu
            do k=1,3
              do l=1,nclu
                ij = nclu*(j-1)+i
                kl = nclu*(l-1)+k
                read(30,*) hessval
                mwhess1(ij,kl)=hessval*mu/dsqrt(mm(i)*mm(k))
              enddo
            enddo
          enddo
        enddo
        open(31, file="hess.2")
        do i=1,3
          do j=1,nclu
            do k=1,3
              do l=1,nclu
                ij=i*j
                kl=k*l
                ij = nclu*(j-1)+i
                kl = nclu*(l-1)+k
                read(31,*) hessval
                mwhess2(ij,kl)=hessval*mu/dsqrt(mm(i)*mm(k))
              enddo
            enddo
          enddo
        enddo
      endif
      write (6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the Hessians used to determine the state counts

      write (6,*)"Surface 1 grad(E1-E2)-projected"
      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
        do j=1,ndim
          hessx(i,j)=mwhess1(i,j)
        enddo
      enddo
      call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp,hessx)
      call dsyev('v','u',ndim,hessx,nmax,freq,work,lwork,info)

      write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)

      write (6,*)"Surface 2 grad(E1-E2)-projected"
      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessx(i,j)=mwhess2(i,j)
      enddo
      enddo
      call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp,hessx)
      call dsyev('v','u',ndim,hessx,nmax,freq,work,lwork,info)

      write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)
      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the Hessians used to determine the state counts
     
c       write(6,*) "Diagonalizing the original mw Hessians"
c       write(6,*) "Hess 1"
c       ndim = 3*nclu
c       nmax = 3*nclu
c       call dsyev('v','u',ndim,mwhess1,nmax,freq,work,lwork,info)
c       do k=1,ndim
c         if (freq(k).gt.0.d0) then
c             tmp=dsqrt(freq(k)/mu)
c         else
c             tmp=-dsqrt(-freq(k)/mu)
c         endif
c         write(6,150)k,freq(k),tmp*autocmi
c       enddo
c       write(6,*) "Hess 2"
c       ndim = 3*nclu
c       nmax = 3*nclu
c       call dsyev('v','u',ndim,mwhess2,nmax,freq,work,lwork,info)
c       do k=1,ndim
c         if (freq(k).gt.0.d0) then
c             tmp=dsqrt(freq(k)/mu)
c         else
c             tmp=-dsqrt(-freq(k)/mu)
c         endif
c         write(6,150)k,freq(k),tmp*autocmi
c       enddo

      write (6,*)"Calculating Effective Two-State Hessians..."
      write (6,*)

c     calculate the effective two-state Hessian from paper      
      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
        hessx(i,j)=(g1norm*mwhess1(i,j) + g2norm*mwhess2(i,j)) / 
     &  g1g2normsum
      enddo
      enddo
      call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp,hessx)
      call dsyev('v','u',ndim,hessx,nmax,freq,work,lwork,info)

      write(6,*)"   Index  Force Const (mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the state counts and write the ne* files for MESS

c numerical grid for nej.dat
      ne=int(emax/es)
      emin=0.  
      emax=emax/autocmi
      emin=emin/autocmi
      jmin=0.  
      hso12=hso12/autocmi

c MSX properties
      nfreq=(3*nclu-7)
      if (linear) nfreq=nfreq+1
      print *,"MSX frequencies from Harvey Effective Hessian"
      ii=0
      rho=1.d0  ! calc classical harmonic rho prefactor
      do k=ndim-nfreq+1,ndim
       tmp=dsqrt(freq(k)/mu)  
       ii=ii+1
       fx(ii)=tmp
       rho=rho/tmp
       print *,tmp*autocmi
      enddo
      do i=2,nfreq-1
      rho=rho/dble(i)
      enddo
      print *

c NEJ.DAT
      open(33,file="nej_lz.dat")
      open(44,file="nej_wc.dat")
c     initialize density of state arrays t(i) (and at(i))
c     t(i) corresponds to energy bin from E = ESTEP*(i-1) to ESTEP*i
      estep=es/autocmi
      ezero=0.d0
      imax = int(emax/estep)+1
      write(33,*)imax,jmax/js+1,1.d0
      write(44,*)imax,jmax/js+1,1.d0

      do jj=jmin,jmax,js

      do i=1,imax
      at(i) = 0.d0
      enddo
      if (linear) then
        djk = dble(2*jj+1)         ! degeneracy for tau
        ejk = mom2*dble(jj*(jj+1))
        iejk = int(ejk/estep)+1
        if (iejk.le.0) then
          write(6,*)"Rotational energy level found below energy grid"
          write(6,*)jj,i,iejk,ejk,estep
          stop
        endif
        if (iejk.le.imax) at(iejk) = at(iejk)+djk
      else
      do k=0,jj
        djk = dble(2*jj+1)        ! degeneracy for j
        if (k.ne.0) djk=djk*2.d0  ! degeneracy for k
        ejk = mom2*dble(jj*(jj+1))+(mom1-mom2)*dble(k**2)
        iejk = int(ejk/estep)+1
        if (iejk.le.0) then
          write(6,*)"Rotational energy level found below energy grid"
          write(6,*)jj,i,iejk,ejk,estep
          stop
        endif
        if (iejk.le.imax) at(iejk) = at(iejk)+djk
      enddo
      endif
      do i=1,imax
        t(i)=at(i)
      enddo
      do j=1,nfreq
        kmax = int(emax/fx(j))  ! maximum quanta w/ E < EMAX
        do k=1,kmax               ! loop over allowed quanta
          ifreq = int(dble(k)*fx(j)/estep) ! vib level spacing in grid units
          do i=1,imax-ifreq
          at(i+ifreq)=at(i+ifreq)+t(i)
          enddo
        enddo
        do i=1,imax
          t(i)=at(i)
        enddo
      enddo

      do i=1,imax
        t_lz(i)=at(i)
        t_ai(i)=at(i)
      enddo
      do i=1,imax
       ee(i)=dble(i)*estep
      enddo
      do i=1,imax
       rholz_lz=0.d0
       rholz_ai=0.d0
       do j=1,imax ! convolute
        if (i.eq.j) then
          etmp = estep/2.d0
        elseif (i.lt.j) then
          etmp = -ee(j-i)+estep/2.d0
        else
          etmp = ee(i-j)+estep/2.d0
        endif
c       Landau-Zener transition probability
        if (j.le.i) then
        plz=1.d0-dexp(-2.d0*pi*hso12**2/tmp3*dsqrt(0.5d0*mu/etmp))  
        p2pass_lz=plz+(1.d0-plz)*plz
        else
        plz=0.d0
        p2pass_lz=0.d0
        endif
c        p2pass_lz=1.d0  ! test
c       Weak-coupling transition probability
        tmpf=dsqrt(tmp1*tmp2)
        e0=(tmpf**4/(2.d0*mu*tmp3**2))**(1.d0/3.d0)
        be=(2.d0*hso12*tmpf/e0/tmp3)**(3.d0/2.d0)
        ex=etmp*tmp3/(2.d0*hso12*tmpf)
        airyarg=-ex*(be**(2.d0/3.d0))
        airypre=pi**2*be**(4.d0/3.d0)
        call airya(airyarg,a1,a2,a3,a4)
        pairy=airypre*a1**2
        rholz_lz=rholz_lz+at(j)*p2pass_lz
        rholz_ai=rholz_ai+at(j)*pairy
       enddo
       t_lz(i)=rholz_lz
       t_ai(i)=rholz_ai
      enddo

      eee=0.d0
      do while(eee.lt.ezero)
      write(33,133)eee*autocmi,jj,0.d0
      write(44,133)eee*autocmi,jj,0.d0
      eee=eee+estep
      enddo
      write(33,133)eee*autocmi,jj,0.d0
      write(44,133)eee*autocmi,jj,0.d0
      do i=1,imax
      eee=dble(i)*estep+ezero
      if (eee.le.emax) write(33,133)eee*autocmi,jj,t_lz(i)*sc_qelec
      if (eee.le.emax) write(44,133)eee*autocmi,jj,t_ai(i)*sc_qelec
      enddo

      enddo

c create NE.DAT from NEJ.DAT by just summing over J
      close(33)
      close(44)
      rewind(33)
      rewind(44)
      open(33,file="nej_lz.dat")
      read(33,*)
      open(44,file="nej_wc.dat")
      read(44,*)
      open(34,file="ne_lz.dat")
      open(45,file="ne_wc.dat")
      do i=1,imax
         t_lz(i)=0.d0
         t_ai(i)=0.d0
      enddo
      do j=0,jmax,js
      do i=1,imax
      read (33,*)eee,jj,xxx
      t_lz(i)=t_lz(i)+xxx
      read (44,*)eee,jj,xxx
      t_ai(i)=t_ai(i)+xxx
      if (j.eq.jmax) write(34,134)eee,t_lz(i)*js
      if (j.eq.jmax) write(45,134)eee,t_ai(i)*js
      enddo
      enddo

      write(6,*)"Writing the nej.dat and ne.dat files..."
      write(6,*)
      write(6,*)

 133  format(f12.2,i12,1pe20.8)
 134  format(f12.2,1pe20.8)
 150  format(i10,e15.5,f15.5)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate cuts along each normal mode if requested

      if (icut.eqv..true.) then
      write(6,*) "Calculating Cuts Along Each Normal Mode..."        
      do im=0,ndim
        if (im.gt.0) print *,"mode=",im,dsqrt(dabs(freq(k))/mu)*autocmi
        if (im.eq.0) print *,"gradient"
        do is=-10,10
          ss=dble(is)*.2d0
          do i=1,3
          do j=1,nclu
            ij = (i-1)*nclu + j
            if (im.gt.0) xx(i,j) = xx0(i,j) + hessx(ij,im)*ss
            if (im.eq.0) xx(i,j) = xx0(i,j) + gperp(ij)*ss
          enddo
          enddo
          do k=1,nclu
          x(k)=xx(1,k)
          y(k)=xx(2,k)
          z(k)=xx(3,k)
          enddo
          call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
          pemd(n1,n1)=pemd(n1,n1)-el_zero
          pemd(n2,n2)=pemd(n2,n2)-el_zero
          print *,ss,pemd(n1,n1)*autoev,pemd(n2,n2)*autoev
        enddo
      enddo
      endif
      write(6,*)
      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine msxfreq
