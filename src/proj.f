      subroutine proj(symb,xx0,mm,mu,natom,nsurf,zero,gvec,hess)

      implicit none

      double precision autoang,amutoau,mu,autocmi
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.52917706d0)
      parameter(amutoau=1822.844987d0)
      integer nsurf
      integer natom,i,j,info,k,ij,l,kl,i1,i2,j1,j2,ndim,nmax,lwork
      double precision x(natom),y(natom),z(natom),mm(natom),
     & pemd(nsurf,nsurf),gpemd(3,natom,nsurf,nsurf),zero,
     & hess(3*natom,3*natom),xx0(3,natom),
     * rot(3,3),eig(3),work2(9),temp1,temp2,temp3,
     & gvec(3*natom),
     & gperp(3*natom),aaa(3*natom,3*natom),freq(3*natom),
     & work(9*natom-1),tmph(3*natom,3*natom),hess2(3*natom,3*natom),
     & gv1a(3,natom),gv2a(3,natom),tmp,hh
      character*2 symb(natom)

c project out gradient
      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        gperp(ij)=gvec(ij)
      enddo
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif 
        enddo
        enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

c translation
      do i=1,3
      do j=1,natom
      ij = (i-1)*natom + j
          gperp(ij)=0.d0
      enddo
      enddo
      do j=1,natom
      ij = (1-1)*natom + j
          gperp(ij)=dsqrt(mm(j))
      enddo
      tmp=0.d0
      do i=1,3*natom
        tmp=tmp+gperp(i)**2
      enddo
      tmp=dsqrt(tmp)
      do i=1,3*natom
        gperp(i)=gperp(i)/tmp
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif
        enddo
        enddo
      enddo
      enddo

      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,natom
      ij = (i-1)*natom + j
          gperp(ij)=0.d0
      enddo
      enddo
      do j=1,natom
      ij = (2-1)*natom + j
          gperp(ij)=dsqrt(mm(j))
      enddo
      tmp=0.d0
      do i=1,3*natom
        tmp=tmp+gperp(i)**2
      enddo
      tmp=dsqrt(tmp)
      do i=1,3*natom
        gperp(i)=gperp(i)/tmp
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif
        enddo
        enddo
      enddo
      enddo

      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,natom
      ij = (i-1)*natom + j
          gperp(ij)=0.d0
      enddo
      enddo
      do j=1,natom
      ij = (3-1)*natom + j
          gperp(ij)=dsqrt(mm(j))
      enddo
      tmp=0.d0
      do i=1,3*natom
        tmp=tmp+gperp(i)**2
      enddo
      tmp=dsqrt(tmp)
      do i=1,3*natom
        gperp(i)=gperp(i)/tmp
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif
        enddo
        enddo
      enddo
      enddo

      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

c rotation
      do i=1,3
      do j=1,natom
      ij = (i-1)*natom + j
          gperp(ij)=0.d0
      enddo
      enddo
      do j=1,natom
      ij = (2-1)*natom + j
c          gperp(ij)=gperp(ij)-xx0(3,j)/dsqrt(mm(j))
          gperp(ij)=gperp(ij)-xx0(3,j)/dsqrt(mu)
      ij = (3-1)*natom + j
c          gperp(ij)=gperp(ij)+xx0(2,j)/dsqrt(mm(j))
          gperp(ij)=gperp(ij)+xx0(2,j)/dsqrt(mu)
      enddo
      tmp=0.d0
      do i=1,3*natom
        tmp=tmp+gperp(i)**2
      enddo
      tmp=dsqrt(tmp)
      do i=1,3*natom
        gperp(i)=gperp(i)/tmp
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif
        enddo
        enddo
      enddo
      enddo

      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,natom
      ij = (i-1)*natom + j
          gperp(ij)=0.d0
      enddo
      enddo
      do j=1,natom
      ij = (1-1)*natom + j
c          gperp(ij)=gperp(ij)+xx0(3,j)/dsqrt(mm(j))
          gperp(ij)=gperp(ij)+xx0(3,j)/dsqrt(mu)
      ij = (3-1)*natom + j
c          gperp(ij)=gperp(ij)-xx0(1,j)/dsqrt(mm(j))
          gperp(ij)=gperp(ij)-xx0(1,j)/dsqrt(mu)
      enddo
      tmp=0.d0
      do i=1,3*natom
        tmp=tmp+gperp(i)**2
      enddo
      tmp=dsqrt(tmp)
      do i=1,3*natom
        gperp(i)=gperp(i)/tmp
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif
        enddo
        enddo
      enddo
      enddo

      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,natom
      ij = (i-1)*natom + j
          gperp(ij)=0.d0
      enddo
      enddo
      do j=1,natom
      ij = (1-1)*natom + j
c          gperp(ij)=gperp(ij)-xx0(2,j)/dsqrt(mm(j))
          gperp(ij)=gperp(ij)-xx0(2,j)/dsqrt(mu)
      ij = (2-1)*natom + j
c          gperp(ij)=gperp(ij)+xx0(1,j)/dsqrt(mm(j))
          gperp(ij)=gperp(ij)+xx0(1,j)/dsqrt(mu)
      enddo
      tmp=0.d0
      do i=1,3*natom
        tmp=tmp+gperp(i)**2
      enddo
      tmp=dsqrt(tmp)
      do i=1,3*natom
        gperp(i)=gperp(i)/tmp
      enddo

      do i=1,3
      do j=1,natom
        ij = (i-1)*natom + j
        do k=1,3
        do l=1,natom
          kl = (k-1)*natom + l
          if (ij.eq.kl) then
             aaa(ij,kl)=1.d0-gperp(ij)*gperp(kl)
          else
             aaa(ij,kl)=-gperp(ij)*gperp(kl)
          endif
        enddo
        enddo
      enddo
      enddo

      do l=1,natom*3
      do k=1,natom*3
      tmph(k,l)=0.d0
      do j=1,natom*3
         tmph(k,l)=tmph(k,l)+aaa(k,j)*hess(j,l)
      enddo
      enddo
      enddo
      do l=1,natom*3
      do k=1,natom*3
      hess(k,l)=0.d0
      do j=1,natom*3
         hess(k,l)=hess(k,l)+tmph(k,j)*aaa(j,l)
      enddo
      enddo
      enddo

      end
