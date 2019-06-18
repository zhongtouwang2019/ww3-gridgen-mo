!  This module is a common block similar in all AFT Model programs and is
!  written in FORTRAN 90.
!                     J G Li   26 Oct 2000
!!
!! Adapted for multiple cell 2D advection tests using UNO schemes.
!!                    J G Li    8 Aug 2007
!! Reformated for global multiple cell advection tests.
!!                    J G Li   16 Nov 2007
!! Modified for extended global SMC grid advection tests.
!!                    J G Li   28 Nov 2008
!! Modified for Arctic 2-D spectral transportation tests.
!!                    J G Li   16 Jul 2009
!! Adapted for UK3km multi-resolution grid spectral transport.
!!                    J G Li    8 Feb 2010
!! Changed for multi-step, multi-resolution G6kSMC grid transport.
!!                    J G Li   23 Feb 2010
!! Adapted for 2-part, 4-step, 8-resolution G6kSMC grid spectral transport.
!!                    J G Li    5 Mar 2010
!! Add great circle turning term in G6kSMC grid spectral transport.
!!                    J G Li   22 Mar 2010
!! Add refraction term and use shallow water wave speed. 
!!                    J G Li   26 Mar 2010
!! Modify refraction term with rotation sub.
!!                    J G Li   18 May 2010
!! Add diffusion term in advection subs.
!!                    J G Li    3 Jun 2010
!! Add Atlantic round patch for comparison with Arctic one.
!!                    J G Li    2 Jun 2011
!! New refraction formulation using cg only.
!!                    J G Li    3 Jun 2011
!! Old refraction plus refraction limiter to gradient direction.
!!                    J G Li   16 Jun 2011
!! Modified for SMC625 global part only spectral transport.
!!                    J G Li   12 Dec 2011
!! Modified to use new cell and face array files.  JGLi28Feb2012
!!
!! Test redued Arctic part on SMC6-25 grid.   JGLi21Mar2012
!!
!! Test refined UK waters on SMC6125 grid.   JGLi08Jan2013
!!
!! Adapted for UK 3km global 25km SMC36125 grid.   JGLi28Feb2014
!!
!! Adapted for global 25km G25SMC grid.   JGLi07Apr2014
!!
!! Add OpenMP directives and Arctic refraction.  JGLi08Jul2015
!!
!! Modifie for UK 3km global 25km SMC36125 grid.   JGLi09Jul2015
!!
!! Modifie for Andy6125 Atlantic grid.    JGLi13Mar2019
!!

                      
      MODULE Constants
          IMPLICIT NONE

! Parameters fixed in the program
    INTEGER,PARAMETER::NCL=900000, NFC=900000, NBDY=256, NDIR=36
!     INTEGER,PARAMETER::NCL=180000, NFC=180000, NBDY=256, NDIR=36,  &
!                    &   MRL=4, NLat=768*2**(MRL-1), NLon=1024*2**(MRL-1), NLat2=NLat/2

!    REAL,PARAMETER:: PoLAT= 0.0, PoLON=-180.0, DLON=0.3515625, DLAT=0.2343750, &
!     REAL,PARAMETER:: PoLAT= 0.0, PoLON=-180.0,                              &
!        &             DLON=(360.0/1024.0)/(2.0**(MRL-1)), DLAT=DLON*2.0/3.0, &
!        &             ZrLAT=-24.257812-DLAT*2.0, ZrLON=-98.261719-DLON*2.0,  &
!        &             Agu36=4.8481E-5, Frqcy=0.0625,   &
!        &             Pie=3.141592654, RAD2D=180.0/Pie, D2RAD=Pie/180.0
     REAL,PARAMETER:: Agu36=4.8481E-5, Frqcy=0.0625,   &
        &             Pie=3.141592654, RAD2D=180.0/Pie, D2RAD=Pie/180.0

!    REAL,PARAMETER:: DT=300.0, DTR=1.0/DT, AKH=9000.0
!    REAL,PARAMETER:: DT=150.0, DTR=1.0/DT, AKH=5000.0
     REAL,PARAMETER:: DT=150.0, DTR=1.0/DT, AKH=100.0

!  Some physical and atmospheric constants
       REAL,PARAMETER:: GRVTY=9.806,CPVAP=1004.5,RDRY=287.05, &
        &    CT0=273.16,CALJO=4.1868,PATM=101325.0,ANGUL=7.2921E-5,  &
        &    EPSLN=0.6220,CLIGHT=2.99792458E8, GeoPie=3.141592654,   &
        &    REARTH=6.371E6

! Array variables to be used for data storage
       REAL:: BXSB, BYSB, BX0, BY0, PoLAT, PoLON, DLAT, DLON, ZrLAT, ZrLON
       REAL::  AMG, CMX, CTT, UMX, DY, DYR, DX0, DThta, SWH0, Alpha
       REAL::  AKHDT2, CGCMX, CRFMX 
       REAL, DIMENSION(-9:NCL):: A, C, D, F, AU, AV, DX, DXR, UC, VC, RCELA
       REAL, DIMENSION(-9:NCL):: HCel, DHDX, DHDY, REFR, CGrp, AngCD, CoGCT, ELaCD
       REAL, DIMENSION( NBDY ):: MBGlo, MBArc
       REAL, DIMENSION( NDIR ):: Theta, Spectr, CSeta, SNeta, SpeGCT
       REAL, DIMENSION(NFC)::   U, V, FU, FV, FX, FY, AngU, AngV
       REAL, DIMENSION(NDIR,NFC)::   UDBn, VDBn
       REAL, DIMENSION(NDIR,NCL)::   UCBn, VCBn, CDBn, CoRfr, DHLMT
       REAL, DIMENSION(:), ALLOCATABLE::  YLat, CSLat, CCLat, DXLat, TnLat 

       INTEGER:: NBLAT, NBLON, NLEVS, MRL, NLon, NLat, NLat2, NSLON, NSLAT 
       INTEGER:: NU, NV, NC, NS, NT, ND, NA, NB, N1, N2, N4, N8, N9
       INTEGER:: NU1, NV1, NU2, NV2, NU4, NV4, NU8, NV8, NU9, NV9
       INTEGER:: NGLo, NGLA, NGLB, NArc, NArA, NArB, NUGL, NUAr, NVGL, NVAr 
       INTEGER:: ICE(4,-9:NCL), KG(NCL), NTS, NWP, MFct, JShft
       INTEGER, DIMENSION(7,NFC)::  ISD
       INTEGER, DIMENSION(8,NFC)::  JSD
       INTEGER, DIMENSION(:), ALLOCATABLE::  NRLCel, NRLUFc, NRLVFc
       INTEGER:: I,II,IJ,IJK,J,JJ,JK,K,KK,L,LL,LM,LMN,M,MM,MN,N,NN
       CHARACTER(LEN=9):: FL9NM='Cn10000.d'

       LOGICAL:: Arctic = .false.
!      LOGICAL:: Arctic = .true.

!      Date and time for timing of program by calling Date_And_Time
       CHARACTER(LEN=10):: CDate, CTime

!  Cell and face array files.
       CHARACTER(LEN=26)::  CelFile='ww3Cels.dat', &
        &                   GISFile='ww3GISide.dat', &
        &                   GJSFile='ww3GJSide.dat', &
        &                   ArcFile='DatGMC/SMC36125BArc.dat', &
        &                   AISFile='DatGMC/S325AISide.dat',   &
        &                   AJSFile='DatGMC/S325AJSide.dat' 


!      CHARACTER(LEN=20):: RUNDATE=' G25SMCAr 07Apr2014 '
       CHARACTER(LEN=20):: RUNDATE='A36125  13Mar2019   '

      END MODULE Constants

!!
!! Adapted for multiple cell 2D advection tests using UNO schemes.
!!                    J G Li   26 Jul 2007
!!
!! Adapted for global multiple cell advection tests using UNO2 scheme.
!!                    J G Li   22 Nov 2007
!!
!! Modified for global SMC grid extended to cover the N Pole.
!!                    J G Li   26 Nov 2008
!!
!! Adapted for 2-part, 3-step, 3-resolution G6kSMC grid spectral transport.
!!                    J G Li    5 Mar 2010
!!
!! Changed for global only SMC625 grid spectral transport.
!!                    J G Li   12 Dec 2011
!!
!! Automatic setting of multi-resolution loops  with MRL and MFct.
!!                    J G Li   28 Feb 2014
!!


       PROGRAM GMC2Dvct 
       USE constants
       IMPLICIT NONE

       INTEGER::  icl, jcl, iuf, juf, ivf, jvf, LvR, Lvm
       REAL:: CNST, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNST8

       NAMELIST /GRID_NML/ NLEVS, NBLAT, NBLON, BXSB, BYSB, BX0, BY0 
       NAMELIST /PROP_NML/ NSLAT, NSLON, PoLAT, PoLON 

!  Read the grid NAMELIST
       OPEN(UNIT=11, FILE='ww3Grid.nml',STATUS='OLD',IOSTAT=nn,ACTION='READ')
       READ (UNIT=11, NML=GRID_NML) 
       READ (UNIT=11, NML=PROP_NML) 
       CLOSE(11)

!  Calculate cell values for highest tier
       NLon = NSLON * 2**(NLEVS-1)
       NLat = NSLAT * 2**(NLEVS-1)
       NLat2 = NLat / 2
       DLON = BXSB / 2.0**(NLEVS-1)
       DLAT = BYSB / 2.0**(NLEVS-1)
       ZrLAT = BY0
       ZrLON = BX0
       MRL = NLEVS

! Set array sizes
       ALLOCATE( YLat(-NLat2:NLat2), CSLat(-NLat2:NLat2),CCLat(-NLat2:NLat2), &
                 DXLat(-NLat2:NLat2),TNLat(-NLat2:NLat2) )
       ALLOCATE( NRLCel(0:NLEVS), NRLUFc(0:NLEVS), NRLVFc(0:NLEVS) )

!  Read Global and Arctic part Multiple-Cell info
       PRINT*, ' Opening '//CelFile//' ...'
       OPEN(UNIT=8, FILE= TRIM(CelFile),STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*,  Celfile//' was not opened! '
          READ (8,*) NGLo, NRLCel(1:MRL) 
       DO J=1,NGLo
          READ (8,*) ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), KG(J)
       END DO
       CLOSE(8)
       PRINT*,  CelFile//' read done ', NGLo, NRLCel(1:MRL) 

!!  Arctic part becomes optional.  JGLi12Dec2011
       IF( Arctic ) THEN
       PRINT*, ' Opening '//ArcFile//' ...' 
       OPEN(UNIT=9, FILE= TRIM(ArcFile),STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*,  ArcFile//' was not opened! '
          READ (9,*) NArc, NArB, NGLB
       DO J=NGLo+1, NGLo+NArc
          READ (9,*) ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), KG(J)
       END DO
       CLOSE(9)
       PRINT*,  ArcFile//' read done  NArc=', NArc

!!  Total Arctic boundary cells
          NB=NArB+NGLB
          PRINT*, ' With Arctic part', NArc, NArB, NGLB

       ELSE
          NArc = 0
          NArB = 0
          NGLB = 0
          PRINT*, ' No Arctic part', NArc, NArB, NGLB
       ENDIF

!  Total cell number will be sum of two parts
       NC = NGLo + NArc

!! Shift Cel(2,*) so that Cel(2,i)=0 on Equator if ZrLat is not on Equator
       JShft = NINT(-ZrLat/DLat + 0.01)
       PRINT*, ' Cell j shifted by ', JShft
       ICE(2,1:NC)=ICE(2,1:NC) - JShft

!!   Set boundary cell counts.  Boundary cells for the global part are at the end
!!   of SMC625Cels.dat and for the Arctic part at the start of SMC625Budy.dat.
!!   Boundary cell will then from NGLo-NGLB+1 to NGLo for lower part and NGLo+1 to NGLo+NArB
!!   NGLA and NArA are the extra numbers to be added for boundary loop 1, NGLB and 1, NArB
       NGLA=NGLo-NGLB
       NArA=NGLo

!  Output a few to check input values
       DO J=1, NC, 10000
          WRITE(6,'(i8,2i6,2i4,i6)') J, ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), KG(J)
       END DO

!! Matching boundary cells for links between global and Arctic parts.
!! Assuming global boundary cells are at the end and Arctic boundary
!! cells are at the begginning of their cell list files.
       IF( Arctic ) THEN

!!   Match global boundary cells with Arctic inner cells
       DO i=1, NGLB
          ii=i+NGLA
!!   Search arctic part its cells to match global part boundary cells
          mm=0
          DO k=NArA+NArB+1, NC-1
             IF(ICE(1,ii) .EQ. ICE(1,k) .AND. ICE(2,ii) .EQ. ICE(2,k)) THEN
                MBGLo(i)=k
                mm=1
             ENDIF
          ENDDO
          IF( mm .EQ. 0 ) PRINT*,' Miss global part boundary cell i=',i
       ENDDO

!!   Match Arctic boundary cells with global inner cells
       DO i=1, NArB
          ii=i+NArA
!!   Search global part to match arctic part boundary cells
          mm=0
          DO k=NGLA-2*NArB, NGLA
             IF(ICE(1,ii) .EQ. ICE(1,k) .AND. ICE(2,ii) .EQ. ICE(2,k)) THEN
                MBArc(i)=k
                mm=1
             ENDIF
          ENDDO
          IF( mm .EQ. 0 ) PRINT*,' Miss Arctic part boundary cell i=',i
       ENDDO
       PRINT*, ' Boundary cells matched for', NGLB, NArB

!!   End of boundary cell matching if (Arctic).
       ENDIF


!    Boundary -9 to 0 cells for cell size 2**n
!    Note the position indice for bounary cell are not used.
       ICE(1,-9:0)=0
       ICE(2,-9:0)=0
       ICE(3,   0)=1
       ICE(4,   0)=1
!!   Restrict boundary cell y-size no more than base cell size 2**(MRL-1).
       mm = 2**(MRL - 1)
       DO i=1,9
          ICE(3,-i)=ICE(3,-i+1)*2
          ICE(4,-i)=MIN(mm, ICE(3,-i))
       ENDDO

!!   Evaluate dx length along latitude in rad
       DX0=DLON*D2RAD
       DO n=-NLat2, NLat2
          YLat(n)=Real( n )*DLAT
          TnLat(n)=TAN(  YLat(n)*D2RAD )
          CSLat(n)=COS(  YLat(n)*D2RAD )
          CCLat(n)=COS( (YLat(n)+0.5*DLAT)*D2RAD )
          DXLat(n)=CCLat(n)*DX0
       ENDDO
     
!!  Define cell area for cell update
       DY=DLAT*D2RAD
       DYR=1.0/DY

       IF( Arctic ) THEN 
!!  North Polar cell NC is a round cell of radius DY*ICE(4,NC), 
!!  equivalent DX=Area/(DY*ICE(4,NC))
!!  The RCELA(NC) represent the net V-flux factor for the polar cell
!!  CSLat will be cancelled as all cell update are divided by the factor.
!!  As polar cell does not need this factor, the factor is included in RCELA
!!  so that it will be cancelled at update when it is divided by the factor.
          DX(NC)=Real(ICE(4,NC))*DY*Pie
          DXR(NC)=1.0/DX(NC)
          RCELA(NC)=DXR(NC)*DX0*CSLat(ICE(2,NC)+ICE(4,NC)/2)/Real( ICE(4,NC) )
          LMN = NC - 1
       ELSE
          LMN = NC
       ENDIF

       WRITE(6, *) "Cell number LMN NC=", LMN, NC

!!   Evaluate other cell dx length in rad and cell area in rad**2
!!   except for the polar (last) cell.
!      DO L=-9,NC-1
       DO L=-9, LMN
          J=ICE(2,L)
          DX(L)=DXLat(J)*Real( ICE(3,L) )
          DXR(L)=1.0/DX(L)
          RCELA(L)=1.0/Real( ICE(3,L)*ICE(4,L) )
       ENDDO

!!  Read sorted ISD JSD variables for global part.
!      OPEN(UNIT=10,FILE='DatGMC/G25SMCGISide.dat',STATUS='OLD',IOSTAT=nn,ACTION='READ')
       PRINT*, ' Opening '//GISFile//' ...'
       OPEN(UNIT=10, FILE= TRIM(GISFile),STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*,  GISFile//' was not opened!'
       READ(10,*) NUGL, NRLUFc(1:MRL)      
       WRITE(6,*) " Read u face numbers NUGL, NRLUFc(1:MRL)"     
       WRITE(6,*)                       NUGL, NRLUFc(1:MRL)      
       DO I=1,NUGL
          READ(10,*)  (ISD(N,I), N=1,7)
       END DO
       CLOSE(10)

       PRINT*, ' Opening '//GJSFile//' ...'
       OPEN(UNIT=11, FILE= TRIM(GJSFile),STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*,  GJSFile//' was not opened!'

       READ(11,*) NVGL, NRLVFc(1:MRL)     
       WRITE(6,*) " Read v face numbers NVGL, NRLVFc(1:MRL) "
       WRITE(6,*)                       NVGL, NRLVFc(1:MRL)  
       DO J=1,NVGL
          READ(11,*)  (JSD(N,J), N=1,8)
       END DO
       CLOSE(11)

!!  Read sorted ISD JSD variables for Arctic part.
       IF( Arctic ) THEN

!      OPEN(UNIT=10,FILE='DatGMC/G25SMCAISide.dat',STATUS='OLD',IOSTAT=nn,ACTION='READ')
       PRINT*, ' Opening '//AISFile//' ...'
       OPEN(UNIT=10, FILE= TRIM(AISFile),STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*,  AISFile//' was not opened!'

       READ(10,*) NUAr
       WRITE(6,*) " Read u face numbers NUAr =", NUAr
       DO I=1,NUAr
          READ(10,*)  (ISD(N,I+NUGL), N=1,7)
       END DO
       CLOSE(10)

!      OPEN(UNIT=11,FILE='DatGMC/G25SMCAJSide.dat',STATUS='UNKNOWN',IOSTAT=nn,ACTION='READ')
       PRINT*, ' Opening '//AJSFile//' ...'
       OPEN(UNIT=11, FILE= TRIM(AJSFile),STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*,  AJSFile//' was not opened!'

       READ(11,*) NVAr
       WRITE(6,*) " Read v face numbers NVAr =", NVAr
       DO J=1,NVAr
          READ(11,*)  (JSD(N,J+NVGL), N=1,8)
       END DO
       CLOSE(11)

!!  Set total face nubmers
       NU=NUGL+NUAr
       NV=NVGL+NVAr

!!  Reset arctic part cell numbers in I/JSD by adding NGLo for positive cells only.
!!  The 0 and negative cells for boundary useage will be shared by the two parts.
       DO I=NUGL+1, NU
          DO M=4,7
             IF(ISD(M,I) > 0) ISD(M,I)=ISD(M,I)+NGLo
          END DO
       END DO

       DO J=NVGL+1, NV
          DO M=4,7
             IF(JSD(M,J) > 0) JSD(M,J)=JSD(M,J)+NGLo
          END DO
       END DO

       WRITE(6,*) " Arctic u v face cell values have been adjusted."

!!  Without Arctic part, set total NU NV equal to global part.
       ELSE
           NUAr = 0
           NVAr = 0
           NU=NUGL
           NV=NVGL
       ENDIF

!!  Shift face array j value as ZrLat is not on Equator.
       PRINT*, ' ISD and JSD j shifted by ', JShft
       ISD(2,1:NU)=ISD(2,1:NU) - JShft
       JSD(2,1:NV)=JSD(2,1:NV) - JShft

!!  Directional bins in rad, shifted by half-width from axises
       DThta=2.0*Pie/FLOAT(NDIR)
       DO K=1, NDIR
          Theta(K)=(FLOAT(K) - 0.5)*DThta
         CSeta(K)=COS(Theta(K))
          SNeta(K)=SIN(Theta(K))
       ENDDO

!!  Convert diffusivity into radian^2 s-1 and multiply by 2*DT
!!  The 2.0 factor is added to cancel the 0.5 factor in the gradient.
       AKHDT2 = 2.0*AKH*DT/(REARTH*REARTH)

!!  Calculate bathymetry gradient DHDX/DY and group speed
        CALL CGDHDXDY

!   Whole array assignment for i=-8,NCL
        C=0.0

!!  Input selection for rotation or deformation test
!100    PRINT*,'Please enter 0 for rotation or 1 for deformation:'
!       READ(5,*) NRoDfm
!       NRoDfm=0

!!  Multiple factor for multi-resolution loops
        MFct=2**(MRL - 1)
        WRITE(6,*) ' Multi-Resolution level and factor=', MRL, MFct

!!  Specify run time steps and writeup interval
        NT=0
        NWP= 12*MFct
        NTS= 60*NWP 
        WRITE(6,*) ' Run/writeup length NTS NWP=', NTS, NWP

!  Initialise AngU/V/CD to be zero.
        AngU=0.0
        AngV=0.0
        AngCD=0.0
        ELaCD=0.0

!  Generate flux face and cell rotation Angle for Arctic cells
        IF( Arctic )  CALL ArctAngd

!  Initialise wave spectra and bin velocities,
!  including GCT and refraction Courant numbers.
        CALL SPECUUVV

!  Read restart time step or set default = 0
!  Restart from NS > 0 overwrite C
!       WRITE(6,*) " Enter restart time step NS, default 0 "
!       READ(5,*)  NS
!       IF (NS .LT. 0) NS=0
        NS = 0
        WRITE(6,*) " Restart time step NS set to be ", NS

!  Save a copy of initial values at writeup time step.
        NT=NS
        write(FL9NM(3:7), FMT='(i5)' )  10000+NT
        OPEN(UNIT=26, FILE=FL9NM, STATUS='NEW',IOSTAT=nn)
        IF(nn /= 0) PRINT*,' File FL9NM was not opened! '
        WRITE(UNIT=6,FMT='(2x,"NT= ",i6,3x,A9)') NT,FL9NM

!!   Filter very small C(n) value so it is greater than E-36
        DO i=1, NC
           IF( Abs(C(i)) .LT. 1.0E-36 ) THEN
               D(i)=SIGN(1.0E-36, C(i))
           ELSE
               D(i)=C(i)
           ENDIF
        ENDDO

!    Central basic cells only, all at end of cell array
        WRITE(UNIT=26, FMT='(2x,2i8)' )  NT, NC
        WRITE(UNIT=26, FMT=7113)  (D(n), n=1,NC)

        CLOSE(26)

!    Define cell and face sub-step counts.  JGLi18Jan2012
         NRLCel(0)=0
         NRLUFc(0)=0
         NRLVFc(0)=0
         DO i=1,MRL-1
            NRLCel(i)=NRLCel(i)+NRLCel(i-1)
            NRLUFc(i)=NRLUFc(i)+NRLUFc(i-1)
            NRLVFc(i)=NRLVFc(i)+NRLVFc(i-1)
         ENDDO
         NRLCel(MRL)=NC
         NRLUFc(MRL)=NU
         NRLVFc(MRL)=NV
         
!     Open files to store writups
        OPEN(UNIT=16,FILE='CMesgs.txt',STATUS='UNKNOWN',IOSTAT=nn, &
        &           ACTION='WRITE')
        IF(nn /= 0) PRINT*,' File CMesgs.txt was not opened! '

!     Header messages and configuration information 
        WRITE(UNIT=16,FMT='(1x/   &
        &  "  Global Multiple-Cell 2-D Spherical Advection Model" /    &
        &  "         SMC Version 2.0   J G  Li  Mar 2010  " /)' )

        CALL DATE_AND_TIME(CDate, CTime)
        WRITE(UNIT=16,FMT="(1x,' Run time date ',A10,2x,A10)") CTime, CDate

        WRITE(UNIT=16,FMT='(1x," Size-1 Units DLON DLAT = ",2f14.10)')  DLON, DLAT
        WRITE(UNIT=16,FMT='(1x," Equatorial PoLat PoLon = ",2f8.2)')  PoLat, PoLon
        WRITE(UNIT=16,FMT='(1x," Standard grid ZrLatLon = ",2f9.5)')  ZrLat, ZrLon
        WRITE(UNIT=16,FMT='(1x," Angular speed Eq Agu36 = ",ES12.3)') Agu36
        WRITE(UNIT=16,FMT='(1x," Horizontal diffusivity = ",f8.1)' )  AKH
        WRITE(UNIT=16,FMT='(1x," Basic time step DT (s) = ",f8.1)' )  DT
        WRITE(UNIT=16,FMT='(1x," Maximum grid speed s-1 = ",ES12.3)') UMX
!       WRITE(UNIT=16,FMT='(1x," Maxm group speed m s-1 = ",f8.3)' )  MAX(CGrp(1:NGLo))
        WRITE(UNIT=16,FMT='(1x," Maximum Courant number = ",f8.3)' )  CMX
        WRITE(UNIT=16,FMT='(1x," Max GCT Courant number = ",f8.3)' )  CGCMX
        WRITE(UNIT=16,FMT='(1x," Max Rfr Courant number = ",f8.3)' )  CRFMX
!       WRITE(UNIT=16,FMT='(1x," Max bathy grad DH/DXDY = ",2f8.3)' ) MAX(ABS(DHDX(1:NGLo))), MAX(ABS(DHDY(1:NGLo)))
        WRITE(UNIT=16,FMT='(1x," Initial integrated SWH = ",f8.3)' )  SWH0
        WRITE(UNIT=16,FMT='(1x," Wave Frequency FrqcyHz = ",f8.4)' )  Frqcy
        WRITE(UNIT=16,FMT='(1x," Total time step no NTS = ",i8)' )  NTS
        WRITE(UNIT=16,FMT='(1x," Restart time step   NS = ",i8)' )  NS
        WRITE(UNIT=16,FMT='(1x," Writeup timestep every = ",i8)' )  NWP
        WRITE(UNIT=16,FMT='(1x," Multi-reso levl factor = ",2i8)')  MRL, MFct
        WRITE(UNIT=16,FMT='(1x," Horizontal cell number = ",2i8)')  NC, NCL
        WRITE(UNIT=16,FMT='(1x," Globl/Arctic cell No.s = ",2i8)')  NGLo, NArc
        WRITE(UNIT=16,FMT='(1x," Globl/Arctic bndy No.s = ",2i8)')  NGLB, NArB
        WRITE(UNIT=16,FMT='(1x," Total number of U-face = ",2i8)')  NU, NFC
        WRITE(UNIT=16,FMT='(1x," Globl/Arctic U-face No = ",2i8)')  NUGL, NUAr
        WRITE(UNIT=16,FMT='(1x," Total number of V-face = ",i8)' )  NV
        WRITE(UNIT=16,FMT='(1x," Globl/Arctic V-face No = ",2i8)')  NVGL, NVAr
        WRITE(UNIT=16,FMT='("Global N9,    N8,    N4,    N2,    N1")')
        WRITE(UNIT=16,FMT='(1x,5i8)')  N9, N8, N4, N2, N1 
        WRITE(UNIT=16,FMT='("Global NU9,   NU8,   NU4,   NU2,   NU1")')
        WRITE(UNIT=16,FMT='(1x,5i8)')  NU9, NU8, NU4, NU2, NU1 
        WRITE(UNIT=16,FMT='("Global NV9,   NV8,   NV4,   NV2,   NV1")')
        WRITE(UNIT=16,FMT='(1x,5i8)')  NV9, NV8, NV4, NV2, NV1 
        WRITE(UNIT=16,FMT='("Sub-step cell count NRLCel(0:MRL)=",6i8)') NRLCel
        WRITE(UNIT=16,FMT='("Sub-step Ufce count NRLUFc(0:MRL)=",6i8)') NRLUFc
        WRITE(UNIT=16,FMT='("Sub-step Vfce count NRLVFc(0:MRL)=",6i8)') NRLVFc

 3912   FORMAT(1x,i4,3F9.1,ES12.3)

        WRITE(16,FMT='(/1x," YLat CSLat at step of  ",i6)' )  10
        WRITE(16,FMT='(8F9.3)')  ( YLat(n), n=-NLat2,NLat2,10)
        WRITE(16,FMT='(8F9.5)')  (CSLat(n), n=-NLat2,NLat2,10)
        WRITE(16,FMT='(/1x," First and last ICE values ")' ) 
        WRITE(16,FMT='(1x,6i8)')  1,(ICE(i, 1),i=1,4)
        WRITE(16,FMT='(1x,6i8)') NC,(ICE(i,NC),i=1,4)
        WRITE(16,FMT='(/1x," First and last ISD values ")' ) 
        WRITE(16,FMT='(1x,8i8)')  1,(ISD(i, 1),i=1,7)
        WRITE(16,FMT='(1x,8i8)') NU,(ISD(i,NU),i=1,7)
        WRITE(16,FMT='(/1x," First and last JSD values ")' ) 
        WRITE(16,FMT='(1x,9i8)')  1,(JSD(i, 1),i=1,8)
        WRITE(16,FMT='(1x,9i8)') NV,(JSD(i,NV),i=1,8)
        WRITE(16,FMT='(1x,8i8)') 

        CALL DATE_AND_TIME(CDate, CTime)
        WRITE(6,"(' Loop start time ',A,'  ',A)")  CTime

!     Start of major time step loop
 TSLoop:  DO  NT=NS, NS+NTS-1, MFct

!!    36 directional bin loop
  DirLop:  DO  ND=1, NDIR

!!    Assign CDBn(ND,1:NC) to C(1:NC)
         C(1:NC)=CDBn(ND,1:NC)

!!    Assign U/VDBn(ND,1:NU/V) to U/V(1:NU/V)
         U(1:NU)=UDBn(ND,1:NU)
         V(1:NV)=VDBn(ND,1:NV)
         UC(1:NC)=UCBn(ND,1:NC)
         VC(1:NC)=VCBn(ND,1:NC)

!! Reset net flux arrays for main step
           F  = 0.0
           AU = 0.0
           AV = 0.0

!! Sub-stepping for different sizes of cells
        DO  NB=1, MFct

!! Loop over 3 refined levels
           DO LvR=1, MRL

!! Sub-level multiple number
              Lvm=2**(LvR-1) 
              CNST2=FLOAT(Lvm)

!! Only calculate the level when MOD(NB, Lvm) = 0
              IF( MOD(NB, Lvm) .EQ. 0 ) THEN

!! Assign sub-time step counts
              icl=NRLCel(LvR-1)+1
              iuf=NRLUFc(LvR-1)+1
              ivf=NRLVFc(LvR-1)+1
              jcl=NRLCel(LvR)
              juf=NRLUFc(LvR)
              jvf=NRLVFc(LvR)

!! Check index ranges at first main time step
           IF( NT .EQ. NS .AND. ND .EQ. 1 ) THEN
              WRITE(6, *) "NB, LvR, Lvm =", NB, LvR, Lvm
           ENDIF

!  Call subroutines to calculate advection
!  Call GMCxUNO2 to calculate MFx value
           CALL GMCxUNO2(iuf, juf, Lvm)

!  Store conservative flux in F advective one in A
!  Add diffusion flux as well, note FX with an negative sign
           DO i=iuf, juf
              M=ISD(5,i)
              N=ISD(6,i)
              F(M) = F(M) - FU(i)*U(i)*CNST2 + FX(i)
              F(N) = F(N) + FU(i)*U(i)*CNST2 - FX(i)
              AU(M) = AU(M) - FU(i)*UC(M)*CNST2 + FX(i)
              AU(N) = AU(N) + FU(i)*UC(N)*CNST2 - FX(i)
           ENDDO

!  Store conservative update in D and advective update in C
!  The side length in MF value has to be cancelled with cell length
!  Also divided by another cell length as U UC is in basic unit
           DO n=icl, jcl
              D(n)=C(n) + F(n)*RCELA(n)
              C(n)=C(n) + AU(n)*RCELA(n)
              F(n)=0.0
              AU(n)=0.0
           ENDDO
!! Note the N Polar cell does not have any U-flux input.

!  Call GMCyUNO2 to calculate MFy value
           CALL GMCyUNO2(ivf, jvf, Lvm)

!  Store conservative flux in F
!  Add diffusion flux FY, note FY with an negative sign
           DO j=ivf, jvf
              M=JSD(5,j)
              N=JSD(6,j)
              AV(M) = AV(M) - FV(j)*V(j)*CNST2 + FY(j)
              AV(N) = AV(N) + FV(j)*V(j)*CNST2 - FY(j)
           ENDDO

!  Store conservative update of D in C 
!  The v side length in MF value has to be cancelled with cell length
!! One cosine factor is also needed to be divided for GMC grid
           DO n=icl,jcl
              CNST=RCELA(n)/CSLat( ICE(2,n)+ICE(4,n)/2 )
              C(n)=D(n) + AV(n)*CNST
              AV(n)=0.0
           ENDDO

!! End of refined level IF( MOD(NB, Lvm) .EQ. 0 ) block
              ENDIF
!! End of refined level loop LvR 
           ENDDO

!! End of sub-step loops NB
        ENDDO

!!    Store C(1:NC) back to CDBn(ND,1:NC)
         CDBn(ND,1:NC)=C(1:NC)

!!    End of directional bin loop
      ENDDO  DirLop

!!    Great circle turning (GCT) for lower part cells at 4 times of substeps,
!!    excluding lower part bounary cells. 
!      CALL GMCGtCrTn(1, NGLA, 4)
!!    Extended to include Arctic part if any
       CALL GMCGtCrTn(1, NC, MFct)

!!    Update boundary cells after proper rotation if Arctic part is
!!    included. 
       IF( Arctic ) THEN

!!    Arctic cells for global boundary cells
       DO i=1,NGLB
          ii=i+NGLA
          kk=MBGLo(i)

!!   Rotate the Arctic spectra by -AnglD before assigning to the lower part
!!   Note that it is equivalent to rotated the directional bins by AnglD.
          Spectr=CDBn(1:NDIR,kk)
          Alpha=  AngCD(kk)

          CALL Specturn( NDir, 1, Alpha, Spectr )

          CDBn(1:NDIR,ii)=Spectr

       ENDDO

!!    Global cells for Arctic boundary cells
       DO i=1,NArB
          ii=i+NArA
          kk=MBArc(i)

!!   Rotate the lower part spectra by AnglD before assigning to the Arctic part
!!   Or turn the directional bins by -AnglD.   21 Jul 2009
          Spectr=CDBn(1:NDIR,kk)
!!   Note only the Arctic cells are assigned the AngCD value
          Alpha= -AngCD(ii)

          CALL Specturn( NDir, 1, Alpha, Spectr )

          CDBn(1:NDIR,ii)=Spectr

       ENDDO

!!   End of updating boundary cells IF( Arctic ). 
       ENDIF


!  Output tracer concentration if at selected writeup time steps
      IF( (NT+MFct .LT. 10*NWP .AND. MOD(NT+MFct,NWP/2) .eq. 0) .OR.    &
     &    (MOD(NT+MFct,NWP) .eq. 0) .OR. (NT+MFct .eq. NWP/4) )  THEN
        write(FL9NM(3:7), FMT='(i5)' )  10000+NT+MFct
        OPEN(UNIT=26, FILE=FL9NM, STATUS='NEW',IOSTAT=nn)
        IF(nn /= 0) PRINT*,' File FL9NM was not opened! '
        WRITE(UNIT=6,FMT='(2x,"NT= ",i6,3x,A9)') NT+MFct,FL9NM

           ii=0
        DO i=1, NC
           CTT=0.0
           DO k=1, NDIR
              CTT = CTT + CDBn(k,i)
           ENDDO
           C(i)=CTT*DThta
           D(i) = SIGN( SQRT( ABS(C(i)) ), C(i) )

!!   Filter very small C(n) value so it is greater than E-36
           IF( Abs(D(i)) .LT. 1.0E-36 ) THEN
               D(i)=SIGN(1.0E-36, C(i))
           ENDIF

        ENDDO

!    All cells are saved 
        WRITE(UNIT=26, FMT='(2x,2i8)' )  NT+MFct, NC
        WRITE(UNIT=26, FMT=7113)  (D(n),  n=1, NC)

        CLOSE(26)
      ENDIF
 7113 FORMAT( 1x, 7ES11.3 )

!!  End of time step loop
      ENDDO  TSLoop

       CALL DATE_AND_TIME(CDate, CTime)
       WRITE(6,"(' Loop ended time ',A,'  ',A)")  CTime

 9999  PRINT*, ' GMC2Dvct completed '

       CALL DATE_AND_TIME(CDate, CTime)
       WRITE(UNIT= 6,FMT="(1x,' End time date ',A10,2x,A10)") CTime, CDate
       WRITE(UNIT=16,FMT="(1x,' End time date ',A10,2x,A10)") CTime, CDate

       END PROGRAM GMC2Dvct 
!  End of main program


! Subroutine that calculate mid-flux values for x dimension 
       SUBROUTINE GMCxUNO2(NUA, NUB, MTP)
         USE Constants
         IMPLICIT NONE
         INTEGER, INTENT(IN):: NUA, NUB, MTP
         REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNST8

!    Two layer of boundary cells are added to each boundary cell face
!    with all boundary cells refer to C(-2:0)=0.0.

!    Notice an extra side length L is multiplied to mid-flux to give correct
!    proportion of flux into the cells.  This length will be removed by the
!    cell length when the tracer concentration is updated.

      DO i=NUA, NUB

!    Diffusion k*dt*2/dx multiplied by MTP for multiple steps
         CNST0=AKHDT2*MTP/DXLat( ISD(2,i) )

!    Select Upstream, Central and Downstream cells
         K=ISD(4,i)
         L=ISD(5,i)
         M=ISD(6,i)
         N=ISD(7,i)

!    Face bounding cell lengths and central gradient
         CNST2=REAL( ICE(3,L) )
         CNST3=Real( ICE(3,M) )
         CNST5=(C(M)-C(L))/(CNST3+CNST2)

!    Courant number in local size-1 cell unit
         CNST6=U(i)*MTP

!    Multi-resolution SMC grid requires flux multiplied by face factor.
         CNST8 = FLOAT( ISD(3,i) )

!    For positive velocity case
         IF(CNST6 >= 0.0)  THEN

!    Upstream cell length and gradient, depending on UFLX sign.
           CNST1=REAL( ICE(3,K) )
           CNST4=(C(L)-C(K))/(CNST1+CNST2)

!    Use minimum gradient all region
           CNST=Sign(1.0, CNST5)*min( Abs(CNST4), Abs(CNST5) )

!    Mid-flux value inside central cell
           FU(i)=(C(L) + CNST*(CNST2-CNST6))*CNST8
           
!    For negative velocity case
         ELSE

!    Upstream cell length and gradient, depending on UFLX sign.
           CNST1=REAL( ICE(3,N) )
           CNST4=(C(N)-C(M))/(CNST3+CNST1)

!    Use minimum gradient outside monotonic region
           CNST=Sign(1.0, CNST5)*min( Abs(CNST4), Abs(CNST5) )

!    Mid-flux value inside central cell M
           FU(i)=(C(M) - CNST*(CNST3+CNST6))*CNST8

         ENDIF

!    Diffusion flux by face gradient x DT x face_width
         FX(i)=CNST0*CNST5*CNST8 

      END DO

! 999  PRINT*, ' Sub GMCxUNO2 ended.'

      RETURN
      END SUBROUTINE GMCxUNO2


! Subroutine that calculate mid-flux values for x dimension 
       SUBROUTINE GMCyUNO2(NVA, NVB, MTP)
         USE Constants
         IMPLICIT NONE
         INTEGER, INTENT(IN):: NVA, NVB, MTP
         REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNST8

!    Two layer of boundary cells are added to each boundary cell face
!    with all boundary cells refer to C(-2:0)=0.0.

!    Diffusion k*dt*2/DY multiplied by MTP for multiple steps
         CNST0=AKHDT2*MTP*DYR

!    Notice an extra side length L is multiplied to mid-flux to give correct
!    proportion of flux into the cells.  This length will be removed by the
!    cell length when the tracer concentration is updated.

      DO j=NVA, NVB

!    Select Upstream, Central and Downstream cells
           K=JSD(4,j)
           L=JSD(5,j)
           M=JSD(6,j)
           N=JSD(7,j)

!    Face bounding cell lengths and gradient
           CNST2=Real( ICE(4,L) )
           CNST3=Real( ICE(4,M) )
           CNST5=(C(M)-C(L))/(CNST3+CNST2)

!    Courant number in basic cell unit
         CNST6=V(j)*MTP

!    Face size integer and cosine factor
         CNST8=CSLat( JSD(2,j) )*Real( JSD(3,j) )

!    For positive velocity case
         IF(CNST6 >= 0.0)  THEN

!    Upstream cell size and irregular grid gradient, depending on VFLY.
           CNST1=Real( ICE(4,K) )
           CNST4=(C(L)-C(K))/(CNST1+CNST2)

!    Use minimum gradient outside monotonic region
           CNST=Sign(1.0, CNST5)*min( Abs(CNST4), Abs(CNST5) )

!    Mid-flux value multiplied by face width and cosine factor
           FV(j)=( C(L) + CNST*(CNST2-CNST6) )*CNST8

!    For negative velocity case
         ELSE

!    Upstream cell size and irregular grid gradient, depending on VFLY.
           CNST1=REAL( ICE(4,N) )
           CNST4=(C(N)-C(M))/(CNST3+CNST1)

!    Use minimum gradient outside monotonic region
           CNST=Sign(1.0, CNST5)*min( Abs(CNST4), Abs(CNST5) )

!    Mid-flux value multiplied by face width and cosine factor
           FV(j)=( C(M) - CNST*(CNST3+CNST6) )*CNST8 

         ENDIF

!    Diffusion flux by face gradient x DT x face_width x cos(lat)
         FY(j)=CNST0*CNST5*CNST8

      END DO

! 999  PRINT*, ' Sub GMCyUNO2 ended.'

      RETURN
      END SUBROUTINE GMCyUNO2


! Subroutine that calculate mid-flux values for x dimension 
       SUBROUTINE CGDHDXDY
         USE Constants
         IMPLICIT NONE
         REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6

!!   Assign water depth to HCel from KG integer values.
!!   set all depth zero for negative cells.
       CGRP(-9:0)=0.0
       DHDX(-9:0)=0.0
       DHDY(-9:0)=0.0
       HCel(-9:0)=0.0
       HCel(1:NC)=FLOAT( KG(1:NC) )
       DO n=1, NC
           IF (HCel(n) .LT. 10.0) HCel(n) = 10.0
       ENDDO

!!   Calculate group speed and refaction factor for all cells
       CNST=2.0*GeoPie*Frqcy
       CNST0=1.0E-8

       DO n=1, NC
            CNST3=HCel(n)
            CNST4=CNST*CNST/GRVTY
!!  Iteration to calculate kH value for this cell
            CNST1=CNST4/TANH(CNST3)
            CNST2=CNST4/TANH(CNST1*CNST3)
            DO WHILE ( ABS(CNST2 - CNST1) .GT. CNST0 )
               CNST1=CNST2
               CNST2=CNST4/TANH(CNST1*CNST3)
            ENDDO
            
            CNST1=CNST2*CNST3 
            CNST2=1.0/COSH(CNST1)
!!  Group speed
            CGrp(n)=(GRVTY*0.5/CNST)*(TANH(CNST1)+CNST2*CNST2*CNST1)

!!  Refraction rate factor, without the dt/d(theta) factor
            REFR(n)=CNST/SINH(2.0*CNST1)
!!  Refraction rate factor, without the dt/d(theta) factor, New formulation
!           REFR(n)=2.0*CNST*( 1.0-CNST3*CNST4 )/( 2.0*CNST1+SINH(2.0*CNST1) )

       ENDDO

!!   Calculate DH/DX using ISD flux array to find neighbouring cells.
!!   0.5 divided by the basic cell x-length at Equator in meter
       CNST6=0.5/(DX0*REARTH)

       AU=0.0
!!   Use the face arrays to calculate the bathymetry x-gradients.
       DO i=1, NU

!    Select Central and Downstream cells
           M=ISD(5,i)
           N=ISD(6,i)

!    Cell length of UCD cells
           CNST2=Real( ICE(3,M) )
           CNST3=Real( ICE(3,N) )

!    Side gradients over basic cell length for central cells include face length factor
           FU(i)=CNST6*ISD(3,i)*(HCel(N)-HCel(M))/(CNST3+CNST2)

!    Store side gradient in two neighbouring cells
           AU(M) = AU(M) + FU(i)
           AU(N) = AU(N) + FU(i)

       END DO

!  Assign averaged side-gradient to DHDX, plus latitude factor
!  Note averaging over 2 times of cell y-width factor. 
!  Excluding polar cell and restrict gradient < 0.1
       DO n=1,NC-1
            CNST3=CCLat( ICE(2,n) )*ICE(4,n)*2.0
            DHDX(n)=AU(n)/CNST3
            IF(DHDX(n) .GT. 0.1) DHDX(n)=0.1
       ENDDO

!! Set polar cell gradient to be zero
            DHDX(NC)=0.0


!!   Calculate DH/DY using JSD flux array to find neighbouring cells.
!!   0.5 divided by the basic cell y-length at Equator in meter
       CNST6=0.5/(DY*REARTH)

       AV=0.0
!!   Use the face arrays to calculate the bathymetry y-gradients.
       DO j=1, NV

!    Select Central and Downstream cells
           M=JSD(5,j)
           N=JSD(6,j)

!    Cell length of UCD cells
           CNST2=Real( ICE(4,M) )
           CNST3=Real( ICE(4,N) )

!    Side gradients over basic cell length for central cells include face length factor
           FV(j)=CNST6*JSD(3,j)*(HCel(N)-HCel(M))/(CNST3+CNST2)

!    Store side gradient in two neighbouring cells
           AV(M) = AV(M) + FV(j)
           AV(N) = AV(N) + FV(j)

       END DO

!  Assign averaged side-gradient to DHDY.
!  Note averaging over 2 times of cell x-width factor. 
!  Excluding the polar cell
       DO n=1,NC-1
            CNST3=ICE(3,n)*2.0
            DHDY(n)=AV(n)/CNST3
            IF(DHDY(n) .GT. 0.1) DHDY(n)=0.1
       ENDDO

!! Set polar cell gradient to be zero
            DHDY(NC)=0.0

!!  Output DHDX DHDY variables for examination
!     WRITE(6,*) " Storing bathy gradient array 1, ", NGLo

!     OPEN(UNIT=10,FILE='G6kmDHDX.d',STATUS='UNKNOWN',IOSTAT=nn)
!     IF(nn /= 0) PRINT*,' File Pros was not opened! '
!     WRITE(10,FMT='(1x,i8)') NGLo
!     WRITE(10,FMT='((10F8.4))')  (DHDX(N), N=1,NGLo)
!     CLOSE(10)

!     OPEN(UNIT=11,FILE='G6kmDHDY.d',STATUS='UNKNOWN',IOSTAT=nn)
!     IF(nn /= 0) PRINT*,' File Pros was not opened! '
!     WRITE(11,FMT='(1x,i8)') NGLo
!     WRITE(11,FMT='((10F8.4))')  (DHDY(N), N=1,NGLo)
!     CLOSE(11)

!     OPEN(UNIT=12,FILE='G6kmCgrp.d',STATUS='UNKNOWN',IOSTAT=nn)
!     IF(nn /= 0) PRINT*,' File Pros was not opened! '
!     WRITE(12,FMT='(1x,i8)') NGLo
!     WRITE(12,FMT='((10F8.3))')  (CGrp(N), N=1,NGLo)
!     CLOSE(12)

!     OPEN(UNIT=13,FILE='G6kmRefr.d',STATUS='UNKNOWN',IOSTAT=nn)
!     IF(nn /= 0) PRINT*,' File Pros was not opened! '
!     WRITE(13,FMT='(1x,i8)') NGLo
!     WRITE(13,FMT='((10F8.4))')  (Refr(N), N=1,NGLo)
!     CLOSE(13)

! 999  PRINT*, ' Sub CGDHDXDY ended.'

       RETURN
       END SUBROUTINE CGDHDXDY


! Subroutine that calculate great circle turning (GCT)
       SUBROUTINE GMCGtCrTn(NCA, NCB, MTP)
         USE Constants
         IMPLICIT NONE
         INTEGER, INTENT(IN):: NCA, NCB, MTP
         REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6

!    Notice an extra side length L is multiplied to mid-flux to give correct
!    proportion of flux into the cells.  This length will be removed by the
!    cell length when the tracer concentration is updated.

      DO n=NCA, NCB

!!   GCT Courant number in MTP times of basic time step without cos(theta)
!!   and it is assumed to be less than 1 in magnitude.
!        CNST0=CoGCT(n)*MTP
!!   Combined into CoRfr

!!   Restrict Courant number to be less than 1.0
!        IF(ABS(CNST0) > 1.0) THEN
!           WRITE(6,*) "GCT C > 1 at cell n=", n, CNST0
!           CNST0=Sign( 1.0, CNST0 )
!        ENDIF

!!   Asign cell spectrum to temporary variable Spcetr
         Spectr=CDBn(1:NDIR,n)
         SpeGCT=0.0 

!!   Loop through NDIR directional bins for each cell spectrum
         DO j=1, NDIR

!    Final GCT Courant number for this dirctional bin
!    Add refraction Courant nubmer as well
!        CNST6=CNST0*CSeta(j) + CoRfr(j,n)*MTP
         CNST6=CoRfr(j,n)*MTP
!!   CoRfr now contains GCT term and refraction. 

!    For positive turning case
         IF(CNST6 > 0.0)  THEN

!    Work out integer number of bins to be skipped.
!    If K is great than NDIR, full circle turning is removed.
           K=MOD( INT(CNST6), NDIR )
 
!    Select the upstream and downstream bins to rotate in, wrap at end
           L=j+K
           M=j+K+1
           IF( L .GT. NDIR ) L = L - NDIR
           IF( M .GT. NDIR ) M = M - NDIR

!!   Divide the j bin energy by fraction of CNST6 and store in SpeGCT
           CNST1=CNST6 - FLOAT( INT(CNST6) )
           CNST2=1.0 - CNST1
           SpeGCT(L)=SpeGCT(L)+Spectr(j)*CNST2
           SpeGCT(M)=SpeGCT(M)+Spectr(j)*CNST1

!    For negative or no turning case
         ELSE 

!    Work out integer number of bins to be skipped.
!    If K is great than NDIR, full circle turning is removed.
           K=MOD( INT(-CNST6), NDIR )

!    Select the upstream and downstream bins to rotate in, wrap at end
           L=j-K
           M=j-K-1
           IF( L .LT. 1 ) L = L + NDIR
           IF( M .LT. 1 ) M = M + NDIR

!!   Divid the bin energy by fraction of CNST6 and store in SpeGCT
           CNST1=-CNST6 - FLOAT( INT(-CNST6) )
           CNST2=1.0 - CNST1
           SpeGCT(L)=SpeGCT(L)+Spectr(j)*CNST2
           SpeGCT(M)=SpeGCT(M)+Spectr(j)*CNST1

         ENDIF

!!   End of directional loop j
         END DO

!!   Store GCT spectrum
         CDBn(1:NDIR,n) = SpeGCT

!!   End of cell loop n
      END DO

! 999  PRINT*, ' Sub GMCGtCrTn ended.'

      RETURN
      END SUBROUTINE GMCGtCrTn


! Subroutine that initialise C U V UC VC for rotating or deform 
!   flow field, using the same grid and time step.

      SUBROUTINE SPECUUVV
        USE Constants
        IMPLICIT NONE
        REAL:: CNST, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNST7
        REAL:: Spec1D(NDIR), Spec2D(NDIR)

!!  All velocities are in unit of basic cell length and time step or 
!!  divided by the grid velocity BX/DT or BY/DT

!! Initialise C U V for spherical rotation test in Arctic

!! Whole array assignment for DDBn(NDir,NCL)
        CDBn=0.0

!! All spectra are identical to a cosine distribution with main direction
!! at 45 deg NE and total wave energy equal to 25.0 so SQRT(E)=5.0.
!! or at -45 deg SE
        CNST=50.0/Pie
        Spec1D = 0.0
        Spec2D = 0.0
        CNST6 = 0.0
        DO k=1, NDIR
           CNST1=COS(Theta(k) + Pie/4.0)
           CNST2=COS(Theta(k) - Pie/4.0)
           IF(CNST1 .GT. 0.0) THEN
              Spec1D(k)=CNST*CNST1*CNST1
           ENDIF
           IF(CNST2 .GT. 0.0) THEN
              Spec2D(k)=CNST*CNST2*CNST2
           ENDIF
           CNST6 = CNST6 + Spec1D(k)
        ENDDO
        SWH0=SQRT(CNST6*DThta)
        
!! Initialise a strip below 45S-60S or j = -768 to -1024 in Southern Ocean 
!! and 45N to 60N or j = 768 to 1024 in upper part of the G6-25SMC domain.
!! Also put a non-zero spectral zone in the Arctic above 86.25N or 16 rows.
!! Note ICE(2,*) is on SW cell corner.
!! Equator is on cell face and ICE is equal to 0 at Equator.
      ij=NINT(45.0/DLat)
      jk=NINT(60.0/DLat)
      mn=NINT(86.25/DLat)

!! N Atlantic round patch centred at (30W 5N)
      ijk=NINT(330.0/(4.0*DLon))*4
      lmn=NINT(  5.0/(4.0*DLat))*4
!! N Atlantic round patch about 3.75 deg latitude in radius, minus
!! half-size-4 cell for cell centre. Note 11 longitudial cells are 
!! used for a radius of 3.868 deg logitude. A ratio of 1.001 is used 
!! to count for the round centre at a face corner.
      CNST3 = 1.001 
      CNST2 = DLat/(3.750-2.0*DLat)
      CNST1 = DLon/(3.868-2.0*DLon)
!! Cell centre is used so round patch radius has to take half-cell away.

      DO i=1,NC

         kk=ICE(2,i)
         mm=ICE(2,i)+ICE(4,i)/2

!!  Southern ocean uses Spec2D
         IF( -jk .LT. kk  .AND.  kk .LT. -ij ) THEN
             CDBn(1:NDIR,i)=Spec2D
             C(i)=SWH0
         ENDIF
!!  Northern Atlantic and Pacific use Spec1D
         IF( ij .LT. kk  .AND.  kk .LT. jk ) THEN
             CDBn(1:NDIR,i)=Spec1D
             C(i)=SWH0
         ENDIF

!! A round patch in N Atlantic to match Arctic one.
!! Use distance to cell centre for symmetry.  JGLi26Mar2012
         CNST5 = (REAL( ICE(1,i)+ICE(3,i)/2 - ijk ))*CNST1*CSLat(mm)
         CNST6 = (REAL( mm                  - lmn ))*CNST2
         CNST4 = CNST5*CNST5 + CNST6*CNST6
         IF( CNST4 .LT. CNST3 ) THEN
             CDBn(1:NDIR,i)=Spec2D
             C(i)=SWH0
         ENDIF

!!  Arctic uses Spec2D as well
         IF( mn .LT. kk ) THEN
             CDBn(1:NDIR,i)=Spec2D
             C(i)=SWH0
         ENDIF

      ENDDO

!  Set angular speed to be 36 hr per cycle
       AMG= Agu36

!!  Factor to convert group speed Cg into Courant number with cell legth and time step
      CNST=DT/REarth

!!  The same angular speed of 36 hr per cycel is used for 24 bin directins
!!  Note dirctions for the Arctic part are deducted by AngU/V in radian.
!!  Face size is for the proportion of face length, not gradient distance. 
      DO L=1, NU
!!  Use average Cg on cell face
         CNST3=0.5*(CGrp( ISD(5,L) ) + CGrp( ISD(6,L) ))*CNST/DXLat( ISD(2,L) )
!!  24 bins are done in one vector assignment 
         UDBn(1:NDIR,L)=CNST3*COS( Theta(1:NDIR)-AngU(L) )
      ENDDO
 
      CNST2=CNST/DY
      DO L=1, NV
         CNST3=0.5*(CGrp( JSD(5,L) ) + CGrp( JSD(6,L) ))*CNST2
         VDBn(1:NDIR,L)=CNST3*SIN( Theta(1:NDIR)-AngV(L) )
      ENDDO
      WRITE(6,*) " Converting U V to Courant number done"

!!  Cell centre U V for advective flux update.
!!  Set boundary cell velocity to be zero
      UC(-9:0)=0.0
      VC(-9:0)=0.0

!!  Find maximum Courant number with cell centre speed
!!  while converting cell centre UC VC to Courant number
!!  in size-1 unit, not really the cell size for other sized cells.
!!  Note that Polar cell DX(NC) is used to store its area
!!  rather than the desired DX(NC), which has no definition.
!!  So UC(NC) is not properly specified here as Courant number.
!!  But this value is never used because no U-flux is associated 
!!  with the Polar cell.  VC(NC) is fine here as DY is cancelled
!!  at the cell update line to give the proper polar cell area.
!!  Note size-1 dx and dy are included in UC and VC, respectively.
!!  Note longitudinal merging factor is divided here as
!!  UC is only divided by the single siz-1 dx.  y-size is
!!  added because subtime steps is proportional to it.

!     CNST=AMG*DT
      CNST=DT/REarth
      CNST1=0.0
      CNST2=0.0
      CNST3=0.0
      WRITE(6,*) " Loop to set max grid speed"
      DO i=1, NC
         UC(i)=CGrp(i)*CNST/DXLAT( ICE(2,i) )
         VC(i)=CGrp(i)*CNST*DYR
         CNST1=Max( CNST1, Abs(UC(i)*ICE(4,i)/ICE(3,i)) )
         CNST2=Max( CNST2, Abs(VC(i)) )
!!  24 bins are done in one vector assignment 
         UCBn(1:NDIR,i)=UC(i)*COS( Theta(1:NDIR)-AngCD(i)*D2RAD )
         VCBn(1:NDIR,i)=VC(i)*SIN( Theta(1:NDIR)-AngCD(i)*D2RAD )
      ENDDO
      WRITE(6,*) " Loop ended"
      CMX=Max(CNST1, CNST2)
!!  Maximum grid speed, i.e. the ratio of speed to grid length
!!  Not necessrily the maximum speed as grid length varies.
      UMX=CMX/DT
      WRITE(6,*) " Max grid speed set"

!!  Great circle tuning Courant number without cos(theta)
      CNST=DT/(REarth*DThta)
      CNST3=0.0

!!  Initialise GCT coefficient to be zero for all
      CoGCT = 0.0

!!  Only specify cells in the lower part (without the Arctic part )
      DO i=1, NGLo
!!  Add GCT term Courant number part 1 (excluding cos(theta))
         CoGCT(i)= -Cgrp(i)*CNST*TnLat( ICE(2,i)+ICE(4,i)/2 )
         CNST3=Max( CNST3, Abs(CoGCT(i)) )
      ENDDO

!!  Add Arctic GCT term.  JGLi09Jul2015
      IF( Arctic ) THEN

      DO i=NGLo+1, NC
!!  Add GCT term Courant number part 1 (excluding cos(theta))
         CoGCT(i)= -Cgrp(i)*CNST*TAN( ELaCD(i)*D2RAD ) 
         CNST3=Max( CNST3, Abs(CoGCT(i)) )
      ENDDO

      ENDIF
      WRITE(6,*) " GCT set"

!!  Maximum GCT Courant number is at main time step (MFct*dt)
      CGCMX=MFct*CNST3

!!  Initialise refraction coefficient to be zero for all
      CoRfr = 0.0
      DHLMT = 0.0

!!  Final refraction Courant number for lower part only
      CNST=DT/DThta
      CNST3=0.0
      DHLMT=0.0

!!  Include Arctic part from NGLo+1 to NC
!     DO i=1, NGLo
      DO i=1, NC 
      DO k=1, NDIR
!!  For Arctic part theta has to be rotated by AngCD(i).  JGLi09Jul2015
         IF( i .GT. NGLo ) THEN
             CNST7=Theta(k) - AngCD(i)*D2RAD
             CNST6=CNST*REFR(i)*( DHDX(i)*SIN(CNST7) - DHDY(i)*COS(CNST7) )
         ELSE
             CNST6=CNST*REFR(i)*( DHDX(i)*SNeta(k) - DHDY(i)*CSeta(k) )
         ENDIF

!!  New refraction limiter to depth gradient direction.  JGLi16Jun2011
!Li   Work out magnitude of depth gradient
         CNST4 = 1.0001*SQRT(DHDX(i)*DHDX(i) + DHDY(i)*DHDY(i))
!Li   Directional depedent depth gradient limiter.  JGLi16Jun2011
         IF ( CNST4 .GT. 1.0E-5 ) THEN
!Li   Refraction is done only when depth gradient is non-zero
            IF( i .GT. NGLo ) THEN
                CNST7=Theta(k) - AngCD(i)*D2RAD
                CNST5 = ACOS(-(DHDX(i)*COS(CNST7) + DHDY(i)*SIN(CNST7))/CNST4 )
            ELSE
                CNST5 = ACOS(-(DHDX(i)*CSeta(k)+DHDY(i)*SNeta(k))/CNST4 )
            ENDIF
            DHLMT(k,i)=0.75*MIN(CNST5,ABS(Pie-CNST5))/DThta

         ENDIF

         CNST1 = Sign( MIN(ABS(CNST6), DHLMT(k,i)), CNST6 )
         CNST3 = Max( CNST3, Abs(CNST1) )

!!  Combining GCT term and refraction term. Note the Arctic part also
!!  the CSeta(k) as GCT uses the rotated theta direction. 
         CoRfr(k,i)= CNST1 + CoGCT(i)*CSeta(k)

      ENDDO
      ENDDO

!!  Maximum refraction Courant number is at main time step (4*dt)
      CRFMX=MFct*CNST3
      WRITE(6,*) '  Max refraction courant number set!'

! 999  PRINT*, ' Sub SPECUUVV ended.'

      RETURN

      END SUBROUTINE SPECUUVV


!  This subroutine turn the wave spectrum by an fixed angle anti-clockwise
!  so that it may be used in the rotated or stanadard system.
!  First created:   26 Aug 2005   Jian-Guo Li
!  Last modified:   20 Jul 2009   Jian-Guo Li
!
! Subroutine Interface:

      Subroutine Specturn( NDirc, NFreq, Alphad, Spectr )
 
! Description:
!   Rotates wave spectrum anticlockwise by angle alphad
!
! Subroutine arguments
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NFreq, NDirc         ! No. frequ and direc bins
       REAL,    INTENT(IN) :: Alphad               ! Turning angle in degree
       REAL, INTENT(INOUT) :: Spectr(NDirc,NFreq)  ! Wave spectrum in and out

! Local variables
       INTEGER :: ii, jj, kk, nsft
       REAL    :: Ddirc, frac, CNST
       REAL, Dimension(NFreq)      ::  Wrkfrq, Tmpfrq
       REAL, Dimension(NDirc,NFreq)::  Wrkspc

! Check input bin numbers
       IF( (NFreq .LT. 0) .OR. (NDirc .LT. 0) )  THEN
          PRINT*, " Invalid bin number NF or ND", NFreq, NDirc
          RETURN
       ELSE
          Ddirc=360.0/FLOAT(NDirc)
       ENDIF

! Work out shift bin number and fraction

      CNST=Alphad/Ddirc
      nsft=INT( CNST )
      frac= CNST - FLOAT( nsft )
!     PRINT*, ' nsft and frac =', nsft, frac

! Shift nsft bins if >=1
        IF( ABS(nsft) .GE. 1 )  THEN
      DO ii=1, NDirc

! Wave spectral direction bin number is assumed to increase clockwise from North
! So shift nsft bins anticlockwise results in local bin number increases by nsft
         jj=ii + nsft
 
! As nsft may be either positive or negative depends on alphad, wrapping may
! happen in either ends of the bin number train
         IF( jj > NDirc )  jj=jj - NDirc
         IF( jj < 1     )  jj=jj + NDirc

! Copy the selected bin to the loop bin number
         Wrkspc(ii,:)=Spectr(jj,:)
 
      Enddo

! If nsft=0, no need to shift, simply copy
        ELSE
        Wrkspc = Spectr
        ENDIF

! Pass fraction of wave energy in frac direction
! Positive or anticlock case, larger bin upstream
        IF( frac > 0.0 ) THEN
      Tmpfrq=Wrkspc(1,:)*frac
      DO kk=NDirc, 1, -1
         Wrkfrq=Wrkspc(kk,:)*frac 
         Spectr(kk,:)=Wrkspc(kk,:) - Wrkfrq + Tmpfrq 
         Tmpfrq=Wrkfrq
      ENDDO
        ELSE
! Negative or clockwise case, smaller bin upstream
      Tmpfrq=Wrkspc(NDirc,:)*frac
      DO kk=1, NDirc
         Wrkfrq=Wrkspc(kk,:)*frac
         Spectr(kk,:)=Wrkspc(kk,:) + Wrkfrq - Tmpfrq
         Tmpfrq=Wrkfrq
      ENDDO
        ENDIF

! Specturn completed

       Return 
       End Subroutine Specturn
!

! Subroutine that generates the Arctic reference direction angle
      SUBROUTINE ArctAngd
        USE Constants
        IMPLICIT NONE
        REAL:: CNST, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6
        REAL, ALLOCATABLE, DIMENSION(:)::  XLon, WLat, ELon, ELat, AnglD
        REAL ::  DfPolat, DfPolon, Omega, DfSpeed

!!    Note only the Arctic part needs the rotation angles.
!     Work out u-face central position XLon, WLat in standard grid
      CNST1=DLon*0.5
      CNST2=DLat*0.5
      DfPolat=Polat
      DfPolon=Polon

      ALLOCATE( XLon(NUAr), WLat(NUAr), ELon(NUAr), ELat(NUAr), AnglD(NUAr) )
      WRITE(6,*) " Calculating U component ..."

!!  Initialisation done before call this sub for Arctic part.
!        AngU=0.0
!        AngV=0.0
!        AngCD=0.0

      DO L=1, NUAr
         i=L+NUGL 
!!  U-face latitude with half dlat increase from SW corner
!!  Longitude is measured from ZrLon where cell i=0. 
!!  Cel and face j values have been shifted so that j=0 is on the !Equator.
         XLon(L)= Float( ISD(1,i) )*DLon + ZrLon 
         WLat(L)= Float( ISD(2,i) )*DLat + Float( ISD(3,i) )*CNST2 
      END DO

!!  Convert standard lat/lon into rotated lat/lon for deformation wind
      CALL LLTOEQANGLE( WLat, XLon, ELat, ELon,     &
     &                 AnglD, DfPolat, DfPolon, NUAr)

      DO L=1, NUAr
         i=L+NUGL 
!!  Convert the AnglD into rad and store in AngU(L). True U value for all 
!!  directional bins will be generated from the value later.
         AngU(i)=AnglD(L)*D2RAD 
      END DO

!!  Output AngU for checking
!     WRITE(6,        *)  "(AngU(L+NUGL), L=1, NUAr)"
!     WRITE(6,'(8ES12.3)') (AngU(L+NUGL), L=1, NUAr)

      DEALLOCATE( XLon, WLat, ELon, ELat, AnglD )

        ALLOCATE( XLon(NVAr), WLat(NVAr), ELon(NVAr), ELat(NVAr), AnglD(NVAr) )

!     Work out v-face central position XLon, WLat in standard grid
      DO L=1, NVAr
         j=L+NVGL
!!  V-face latitude with half_dlon*JSD(3) increase from SW corner
!!  Longitude is measured from ZrLon where cell i=0. 
!!  Cel and face j values have been shifted so that j=0 is on the !Equator.
         XLon(L)= Float( JSD(1,j) )*DLon+CNST1*Float( JSD(3,j) )+ZrLon 
         WLat(L)= Float( JSD(2,j) )*DLat
      END DO

!!  Convert standard lat/lon into rotated lat/lon for deformation wind
      CALL LLTOEQANGLE( WLat, XLon, ELat, ELon,     &
     &                 AnglD, DfPolat, DfPolon, NVAr)

      DO L=1, NVAr
         j=L+NVGL
!!  Convert the AnglD into rad and store in AngV(L). True V value for all 
!!  directional bins will be generated from the value later.
         AngV(j)= AnglD(L)*D2RAD
      END DO

!!  Output AngV for checking
!     WRITE(6,        *)  "(AngV(L+NVGL), L=1, NVAr)"
!     WRITE(6,'(8ES12.3)') (AngV(L+NVGL), L=1, NVAr)

      DEALLOCATE( XLon, WLat, ELon, ELat, AnglD )


!! Specific cell centre velocity compenents
      ALLOCATE( XLon(NArc), WLat(NArc), ELon(NArc), ELat(NArc), AnglD(NArc) )
      WRITE(6,*) " Calculating UC, VC component ..."

      CNST1=DLon*0.5
      CNST2=DLat*0.5
!! All cells include the polar cell
!! Note the wlat is not at 90N for the polar cell as direction will be undefined.
!! Here Wlat is half dlat from the polar cell edge and half dlat from the NP.
      DO L=1, NArc-1
         i=L+NGLo

!!  Cell centre latitude equal to west side centre latitude.
!!  Cell centre longitude with half cell width increase from West side centre
!!  Although the polar cell is of angular radius dlat (not dlat/2) the 
!!  transformation location is still used dlat/2 from its SW corner. The error
!!  will be negeligible as only the AnglD is used.
         XLon(L)= Float( ICE(1,i) )*DLon + CNST1*Float( ICE(3,i) )+ZrLon 
         WLat(L)= Float( ICE(2,i) )*DLat + CNST2*Float( ICE(4,i) )

      END DO

!! North Polar cell centre coincide with NP
         XLon(NArc)=0.0
         WLat(NArc)=90.0
!! AnglD will be undefined with NP location as no local east at NP.

!!  Convert standard lat/lon into rotated lat/lon for deformation wind
      CALL LLTOEQANGLE( WLat, XLon, ELat, ELon,     &
     &                 AnglD, DfPolat, DfPolon, NArc )

      DO L=1, NArc
         i=L+NGLo
!!  Keep the AnglD in Deg and store in AngCD(L).  Spectral rotation for
!!  boundary cell update will use this angle later.
         AngCD(i)=  AnglD(L) 
!!  Save rotated latitude for refraction in Arctic part.  JGLi08Jul2015
         ELaCD(i)=  ELat(L) 
      END DO

!!  Output AngCD for checking
!     WRITE(6,        *)  "(AngCD(L+NGLo), L=1, NArc)"
!     WRITE(6,'(8ES12.3)') (AngCD(L+NGLo), L=1, NArc)

 999  PRINT*, ' Sub ArctAngd ended.'

      RETURN

      END SUBROUTINE ArctAngd


!Li
!Li  Merged UM source code for rotated grid, consiting the following
!Li  original subroutines in UM 6.1
!Li    LLTOEQ1A  WCOEFF1A  and  LBCROTWINDS1
!Li  The last subroutine is modified to process only one level winds
!Li  cpp directives are removed and required header C_Pi.h inserted.
!Li	    Jian-Guo Li     26 May 2005
!Li
!Li  The WCOEFF1A subroutine is merged into LLTOEQ to reduce repetition
!Li  of the same calculations. Subroutine interface changed to 
!Li  LLTOEQANGLE
!Li	    Jian-GUo Li     23 Aug 2005
!Li
!LL  Subroutine LLTOEQANGLE--------------------------------------------    
!LL                                                                        
!LL  Purpose:  Calculates latitude and longitude on equatorial             
!LL            latitude-longitude (eq) grid used in regional               
!LL            models from input arrays of latitude and                    
!LL            longitude on standard grid. Both input and output           
!LL            latitudes and longitudes are in degrees.                    
!Li	       Also calculate rotation angle in degree to tranform
!Li            standard wind velocity into equatorial wind.
!Li	       Valid for 0<PHI_POLE<90 or new pole in N. hemisphere.
!LL                                                                        
!* Arguments:--------------------------------------------------------    
      SUBROUTINE LLTOEQANGLE( PHI, LAMBDA, PHI_EQ, LAMBDA_EQ,     &  
     &                 ANGLED, PHI_POLE, LAMBDA_POLE, POINTS )       

      IMPLICIT NONE 

      INTEGER:: POINTS    !IN  Number of points to be processed             

      REAL :: PHI_POLE,  & !IN  Latitude of equatorial lat-lon pole
     &        LAMBDA_POLE  !INOUT  Longitude of equatorial lat-lon pole

      REAL, DIMENSION(POINTS) ::         &
     &        PHI,       & !IN  Latitude
     &        LAMBDA,    & !IN  Longitude
     &        ANGLED,    & !OUT turning angle in deg for standard wind
     &        LAMBDA_EQ, & !OUT Longitude in equatorial lat-lon coords
     &        PHI_EQ       !OUT Latitude in equatorial lat-lon coords

! Define local varables:-----------------------------------------------
      REAL :: A_LAMBDA, A_PHI, E_LAMBDA, E_PHI, SIN_PHI_POLE, COS_PHI_POLE,  &
     &        TERM1, TERM2, ARG, LAMBDA_ZERO, LAMBDA_POLE_KEEP
      INTEGER   :: I   
      REAL, PARAMETER :: SMALL=1.0E-6

! Constants from comdecks:---------------------------------------------

      Real, Parameter :: Pi = 3.14159265358979323846  , &
     &                   Pi_Over_180 = Pi/180.0       , &
     &                   Recip_Pi_Over_180 = 180.0/Pi        

!*----------------------------------------------------------------------   

! 1. Initialise local constants
! Scale lambda pole to range -180 to 180 degs
      LAMBDA_POLE_KEEP=LAMBDA_POLE
      IF (LAMBDA_POLE.GT. 180.0) then
          LAMBDA_POLE=LAMBDA_POLE-360.0
      ENDIF

! Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.0
! Sine and cosine of latitude of eq pole
      IF (PHI_POLE >= 0.0) THEN
        SIN_PHI_POLE =  SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE =  COS(PI_OVER_180*PHI_POLE)
      ELSE
        SIN_PHI_POLE = -SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE = -COS(PI_OVER_180*PHI_POLE)
      ENDIF

! 2. Transform from standard to equatorial latitude-longitude

      DO 200 I= 1, POINTS

! Scale longitude to range -180 to +180 degs

      A_LAMBDA=LAMBDA(I)-LAMBDA_ZERO
      IF(A_LAMBDA.GT. 180.0) A_LAMBDA=A_LAMBDA-360.
      IF(A_LAMBDA.LE.-180.0) A_LAMBDA=A_LAMBDA+360.

! Convert latitude & longitude to radians

      A_LAMBDA=PI_OVER_180*A_LAMBDA
      A_PHI=PI_OVER_180*PHI(I)

! Compute eq latitude using equation (4.4)

      ARG=-COS_PHI_POLE*COS(A_PHI)*COS(A_LAMBDA)   &
     &    +SIN_PHI_POLE*SIN(A_PHI)
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      E_PHI=ASIN(ARG)
      PHI_EQ(I)=RECIP_PI_OVER_180*E_PHI

! Compute eq longitude using equation (4.6)

      TERM1 = SIN_PHI_POLE*COS(A_PHI)*COS(A_LAMBDA)   &
     &       +COS_PHI_POLE*SIN(A_PHI)
      TERM2 = COS(E_PHI)
      IF(TERM2 .LT. SMALL) THEN
        E_LAMBDA=0.0
      ELSE
        ARG=TERM1/TERM2
        ARG=MIN(ARG, 1.0)
        ARG=MAX(ARG,-1.0)
        E_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
        E_LAMBDA=SIGN(E_LAMBDA,A_LAMBDA)
      ENDIF

! Scale longitude to range 0 to 360 degs

      IF(E_LAMBDA.GE.360.0) E_LAMBDA=E_LAMBDA-360.0
      IF(E_LAMBDA.LT.  0.0) E_LAMBDA=E_LAMBDA+360.0
      LAMBDA_EQ(I)=E_LAMBDA

!Li  Calculate turning angle for standard wind velocity

      E_LAMBDA=PI_OVER_180*LAMBDA_EQ(I)

! Formulae used are from eqs (4.19) and (4.21)

      TERM2=SIN(E_LAMBDA)
      ARG= SIN(A_LAMBDA)*TERM2*SIN_PHI_POLE      &
     &    +COS(A_LAMBDA)*COS(E_LAMBDA)
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      TERM1=RECIP_PI_OVER_180*ACOS(ARG)
      ANGLED(I)=SIGN(TERM1,TERM2)
!Li

 200  CONTINUE

! Reset Lambda pole to the setting on entry to subroutine
      LAMBDA_POLE=LAMBDA_POLE_KEEP

      RETURN
      END SUBROUTINE LLTOEQANGLE
!Li

