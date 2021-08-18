      program nst

      implicit none

c     Parameters for the Job      
      double precision autoang,amutoau
      parameter(autoang=0.52917706d0)
      parameter(amutoau=1822.844987d0) 
      integer mnat,mnsurf
      parameter(mnat=13)    ! maximum number of atoms
      parameter(mnsurf=3)   ! maximum number of electronic surfaces in the PES call

c     Variables for Job Control      
      character*5 jobtype
      logical :: barbor=.true.  ! Turn on Barzilai-Borwein step-size control by default
      integer :: qcprog=1       ! Gaussian as default program package
      integer nargs
c      character*4,allocatable::args(:)

c     Molecule Variables
      integer nclu
      character*2 symb(mnat)
      double precision xx(3,mnat)
      double precision x(mnat),y(mnat),z(mnat)
      double precision mm(mnat)

c     Variables read from input file
      double precision el_zero
      double precision es,emax
      integer js,jmax
      double precision hso12
      double precision sc_qelec
      logical icut

c     Loop Variables      
      integer i,j

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Obtain arguments from command line (Does this do anything?)

c     Sets the optimizer and the program to use
c      nargs=COMMAND_ARGUMENT_COUNT()
c      if (nargs .gt. 0) then
c        allocate(args(nargs))
c        do i=1,nargs
c          call GET_COMMAND_ARGUMENT (i, args(i))
c          args(i)=adjustl(args(i))
c          if (trim(args(i))=='nobb' .or. trim(args(i))=='NOBB' .or.
c     &      trim(args(i))=='noBB' .or. trim(args(i))=='NoBB') then
c            write(*,*) 'Barzilai-Borwein method switched off.'
c            barbor=.false.
c          elseif (trim(args(i))=='m' .or. trim(args(i))=='M') then
c            write(*,*) 'Using Molpro interface instead of Gaussian.'
c            qcprog=2
c          endif
c        enddo
c        deallocate(args)
c      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Reads the Input File for the Job Information
      read(5,*)
c     get the jobtype
      read(5,*)
      read(5,*) jobtype
c     get the number of atoms
      read(5,*)
      read(5,*) nclu
      if (nclu.gt.mnat) then
        print *,"nclu (",nclu,") > mnat (",mnat,")"
        stop
      endif
      do i=1,nclu
        read(5,*)symb(i),(xx(j,i),j=1,3)
        x(i)=xx(1,i)/autoang
        y(i)=xx(2,i)/autoang
        z(i)=xx(3,i)/autoang
      enddo
c     get the masses
      read(5,*)
      do i=1,nclu
        read(5,*)mm(i)
        mm(i)=mm(i)*amutoau
      enddo
c     get the zero of energy
      read(5,*)
      read(5,*) el_zero
c     get energy grid spacing
      read(5,*)
      read(5,*) es,emax
c     get angular grid spacing
      read(5,*)
      read(5,*) js,jmax
c     get so-coupling
      read(5,*)
      read(5,*) hso12
c     get scale factor for elec part'n fxn
      read(5,*)
      read(5,*) sc_qelec
c     get logical variable for normal mode cuts 
      read(5,*)
      read(5,*) icut


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calls Appropriate Subroutine based on Job Type

      if (jobtype.eq."OPTGM") then
        write(6,*) "-------------------------------"
        write(6,*) "Request to Optimize MSX Geometry"
        write(6,*) "-------------------------------"
        write(6,*) 
        call opt(symb,x,y,z,nclu,mnat,mnsurf,el_zero,barbor) 
      end if
      
      if (jobtype.eq."ROTGM") then
        write(6,*) "------------------------------------------------"
        write(6,*) "Request to Rotate MSX Geometry to Principal Axes"
        write(6,*) "------------------------------------------------"
        write(6,*) 
        call msxfreq(xx,mm,nclu,symb,el_zero,jobtype,icut,
     &               es,emax,js,jmax,hso12,sc_qelec)
      end if

      if (jobtype.eq."HESSR") then
        write(6,*) "-------------------------------------------------"
        write(6,*) "Request to Calculate MSX Freqs By Read'g Hessians"
        write(6,*) "-------------------------------------------------"
        write(6,*) 
        call msxfreq(xx,mm,nclu,symb,el_zero,jobtype,icut,
     &               es,emax,js,jmax,hso12,sc_qelec)
      end if
      
      if (jobtype.eq."HESSC") then
        write(6,*) "---------------------------------------------------"
        write(6,*) "Request to Determine MSX Freqs By Calc'ing Hessians"
        write(6,*) "---------------------------------------------------"
        write(6,*) 
        call msxfreq(xx,mm,nclu,symb,el_zero,jobtype,icut,
     &               es,emax,js,jmax,hso12,sc_qelec)
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Exits the Program with Happy Message
      write(6, *) "Exiting Program..."

      end


