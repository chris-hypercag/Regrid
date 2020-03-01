      PROGRAM REGRID
!
!  Generate new mov files for modified domain 
!  and MPI zones.
!
      IMPLICIT NONE
      INTEGER(4) :: INCAL, JNCAL, KNCAL, IPL, JPL, KPL, &
                    MY_ID, ITEMP_ID, I, J, NDEC_ID, JDEC_ID, & 
                    IL, JL, IL1_G, JL1_G, K, KL, &         
                    NADV, II, JJ, ILN
      INTEGER(4) :: INODE, JNODE, IBC0, IBC1, JBC0, JBC1, &
                    IETR, JETR, JDIFF, &
                    INCOM, JNCOM, M, NX, NY
      REAL(8) :: xx1, CCFL, TIME, XX, GBETA
      REAL(8) :: GI1, GI2, GI3, TT
      REAL(8) :: BC1, Z1, TZ, TZ1, ERR, MAXERR
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: RHOS,RHOU,RHOV,RHOW,E, &
                                                X, Y, Z
      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: HSHK, DST
      REAL(8), ALLOCATABLE, DIMENSION(:)   :: RATEJ,    RATEI,   &
                                              RATEJNEW, RATEINEW, &
                                              RATEITEMP, RATEJTEMP
      CHARACTER(64) :: FILEN
      INTEGER(4) :: IPLNEW,JPLNEW,INCALNEW,JNCALNEW,INODENEW,JNODENEW
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: XNEW, YNEW, ZNEW, &
                              RHOSNEW, RHOUNEW,RHOVNEW,RHOWNEW,ENEW
      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: HSHKNEW, DSTNEW, RHOTEMP
                                          

      WRITE(*,*) 'CODE STARTED'

      WRITE(*,*) 'READ: DATA1'
      OPEN(9,FILE='DATA1')
      READ(9,*) &
          xx1,     xx1, &
          xx1,     xx1,     xx1,     xx1, &
          INCAL,   JNCAL,   KNCAL, &
     
          xx1,     xx1,     xx1,     xx1, &
          xx1,     xx1,     xx1, &
          xx1,     xx1, &
          xx1,     xx1,     xx1,     xx1, &
          xx1,     xx1, &
          xx1,     xx1,     ILN, &
          xx1,     xx1,     xx1,     GBETA, &
          xx1,     xx1,     xx1, &
          xx1,     xx1,     xx1,     xx1, &
          xx1,     xx1,     xx1,     xx1, &
          xx1,     xx1,     xx1,     xx1, &
          xx1,     xx1 
      CLOSE(9)

      OPEN(10,FILE='DATA1.WEDGECONE')
      READ(10,*) xx1,  xx1,  xx1,  xx1,  xx1, &
                 GI1,  GI2, &
                 xx1,  xx1
      CLOSE(10)

      WRITE(*,*) 'READ: DATA1.MPI'
      OPEN(11,FILE='DATA1.MPI')
      READ(11,*) &
          INODE, JNODE, IBC0, IBC1 , JBC0, JBC1, &
          IETR, JETR, JDIFF, &
          INCOM, JNCOM
      CLOSE(11)

      WRITE(*,*) 'READ: DATA3.MPI'
      OPEN(13,FILE='DATA3.MPI')
      READ(13,*) INODENEW, JNODENEW, &
                 INCALNEW, JNCALNEW
      CLOSE(13)
 
      IPL = INCAL*INODE + IBC0 + IBC1
      JPL = JNCAL*JNODE + JBC1
      KPL = KNCAL

      WRITE(*,*) 'Allocation dimension:', IPL, JPL, KPL

      ALLOCATE(RHOS(IPL,JPL,KPL),RHOU(IPL,JPL,KPL),RHOV(IPL,JPL,KPL), &
               RHOW(IPL,JPL,KPL),   E(IPL,JPL,KPL), &
               HSHK(IPL,KPL),     DST(IPL,KPL), &
               X(IPL,JPL,KPL), Y(IPL,JPL,KPL), Z(IPL,JPL,KPL), &
               RATEI(IPL), RATEJ(JPL) ) 

      WRITE(*,*) 'Number of MPI zones:', (INODE*JNODE)

      ! Obtain frustum corner (XXB,YYB)
      OPEN(12,FILE="previous-endzone-grids-bottom.dat")
      READ(12,*) IL, JL, KL
      II = IL-ILN
      DO K=1,KL
        DO J=1,JL
          DO I=1,IL
            IF (I.EQ.II.AND.J.EQ.1) THEN
              READ(12,*) X(1,1,K), Y(1,1,K), Z(1,1,K)
            ELSE 
              READ(12,*) xx1, xx1, xx1
            END IF
          END DO
        END DO
      END DO
      CLOSE(12)

      ! Set up J-grid stretching
      XX = 1.0D0 / DBLE(JPL-1)
      BC1=(GBETA+1.0D0)/(GBETA-1.0D0)
      DO J=1,JPL
        Z1=1.0D0-(DBLE(J-1))/(DBLE(JPL-1))
        TZ=BC1**Z1
        TZ1=GBETA*(TZ-1.0D0)/(TZ+1.0D0)
        RATEJ(J)=1.0D0-TZ1
      END DO

      ! Set up I-grid stretching
      GI3 = 1.0D0-GI1-GI2
      XX  = 1.0D0/DBLE(IPL-1)
      DO I=1,IPL
        TT = DBLE(I-1)*XX
        RATEI(I) = (GI1+(GI2+GI3*TT*TT)*TT*TT)*TT
      END DO

      !
      ! Start reading-in data
      !
      DO MY_ID = 0,(INODE*JNODE-1)
        DO I = 1,INODE
          DO J = 1,JNODE
            ITEMP_ID = (I-1)+INODE*(J-1)
            IF (MY_ID.EQ.ITEMP_ID) THEN
              NDEC_ID = I-1 ! Column numbering
              JDEC_ID = J-1 ! Depth numbering
              GO TO 1010
            END IF
          END DO
        END DO

1010    CONTINUE

!        WRITE(*,*) 'WORKING ON:', MY_ID

        ! Set global index for end points of each MPI zone
        IF (NDEC_ID.eq.0) THEN
          ! Inlet
          IL1_G = 1
        ELSE
          ! Interior columns and exit
          IL1_G = IBC0 + NDEC_ID*INCAL - IETR + 1
        END IF

        IF (JDEC_ID.EQ.0) THEN
          ! Surface 
          JL1_G = 1
        ELSE
          ! Interior rows and shock 
          JL1_G = JDEC_ID*JNCAL - JETR + JDIFF + 1
        END IF

!        WRITE(*,*) 'Global Index positions',IL1_G, JL1_G
!        WRITE(*,*) 'MPI indicies', IL, JL, 4

        II = IL1_G-1
        JJ = JL1_G-1

        WRITE(FILEN,'("mov-",I3.3,".dat")') MY_ID
        OPEN(11,FILE=TRIM(FILEN),FORM='unformatted')
        READ(11) NADV, CCFL, TIME, IL, JL, KL
        READ(11) (((RHOS(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        READ(11) (((RHOU(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        READ(11) (((RHOV(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        READ(11) (((RHOW(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        READ(11) (((   E(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        READ(11) (( HSHK(II+I,K),I=1,IL),K=1,KL)
        READ(11) (( DST(II+I,K),I=1,IL),K=1,KL)
        CLOSE(11)
      END DO

      ! Construct New Grids
      IPLNEW = INCALNEW*INODENEW + IBC0 + IBC1
      JPLNEW = JNCALNEW*JNODENEW + JBC1

      WRITE(*,*) 'New dimension:', IPLNEW, JPLNEW, KPL
      WRITE(*,*) 'New MPI zones:', INODENEW*JNODENEW

      M=10
      WRITE(*,*) 'Interpolation Order, M=',M
      ALLOCATE(XNEW(IPLNEW,JPLNEW,KPL),YNEW(IPLNEW,JPLNEW,KPL), &
               ZNEW(IPLNEW,JPLNEW,KPL), &
               HSHKNEW(IPLNEW,KPL), DSTNEW(IPLNEW,KPL), &
               RATEINEW(IPLNEW), RATEJNEW(JPLNEW), &
               RHOSNEW(IPLNEW,JPLNEW,KPL), &
               RHOUNEW(IPLNEW,JPLNEW,KPL), &
               RHOVNEW(IPLNEW,JPLNEW,KPL), &
               RHOWNEW(IPLNEW,JPLNEW,KPL), &
                  ENEW(IPLNEW,JPLNEW,KPL), &
               RATEITEMP(M),RATEJTEMP(M),RHOTEMP(M,M) )

      ! Set up New J-grid stretching
      XX = 1.0D0 / DBLE(JPLNEW-1)
      BC1=(GBETA+1.0D0)/(GBETA-1.0D0)
      DO J=1,JPLNEW
        Z1=1.0D0-(DBLE(J-1))/(DBLE(JPLNEW-1))
        TZ=BC1**Z1
        TZ1=GBETA*(TZ-1.0D0)/(TZ+1.0D0)
        RATEJNEW(J)=1.0D0-TZ1
      END DO

      ! Set up I-grid stretching
      GI3 = 1.0D0-GI1-GI2
      XX  = 1.0D0/DBLE(IPLNEW-1)
      DO I=1,IPLNEW
        TT = DBLE(I-1)*XX
        RATEINEW(I) = (GI1+(GI2+GI3*TT*TT)*TT*TT)*TT
      END DO

      ! Begin Interpolation
      MAXERR = 0.0D0
      DO I = 1,IPLNEW
        CALL HUNT(RATEI,IPL,RATEINEW(I),II)
        NX=MIN(MAX(II-(M-1)/2,1),IPL+1-M)
        RATEITEMP(1:M) = RATEI(NX:NX+M-1)

        CALL POLINT(RATEITEMP,HSHK(NX:NX+M-1,1),M, &
                    RATEINEW(I),HSHKNEW(I,1),ERR)
        MAXERR = MAX(ERR,MAXERR)

        CALL POLINT(RATEITEMP,DST(NX:NX+M-1,1),M, &
                    RATEINEW(I),DSTNEW(I,1),ERR)
        MAXERR = MAX(ERR,MAXERR)

        HSHKNEW(I,:)=HSHKNEW(I,1)
        DSTNEW(I,:)=DSTNEW(I,1)

        DO J = 1,JPLNEW
          CALL HUNT(RATEJ,JPL,RATEJNEW(J),JJ)
          NY=MIN(MAX(JJ-(M-1)/2,1),JPL+1-M)
          RATEJTEMP(1:M) = RATEJ(NY:NY+M-1)
          DO K = 1,KPL

            RHOTEMP(1:M,1:M) = RHOS(NX:NX+M-1,NY:NY+M-1,K)
            CALL POLIN2(RATEITEMP,RATEJTEMP,RHOTEMP,M,M, &
                RATEINEW(I),RATEJNEW(J),RHOSNEW(I,J,K),ERR)
            MAXERR = MAX(ERR,MAXERR)

            RHOTEMP(1:M,1:M) = RHOU(NX:NX+M-1,NY:NY+M-1,K)
            CALL POLIN2(RATEITEMP,RATEJTEMP,RHOTEMP,M,M, &
                RATEINEW(I),RATEJNEW(J),RHOUNEW(I,J,K),ERR)
            MAXERR = MAX(ERR,MAXERR)

            RHOTEMP(1:M,1:M) = RHOV(NX:NX+M-1,NY:NY+M-1,K)
            CALL POLIN2(RATEITEMP,RATEJTEMP,RHOTEMP,M,M, &
                RATEINEW(I),RATEJNEW(J),RHOVNEW(I,J,K),ERR)
            MAXERR = MAX(ERR,MAXERR)

            RHOTEMP(1:M,1:M) = RHOW(NX:NX+M-1,NY:NY+M-1,K)
            CALL POLIN2(RATEITEMP,RATEJTEMP,RHOTEMP,M,M, &
                RATEINEW(I),RATEJNEW(J),RHOWNEW(I,J,K),ERR)
            MAXERR = MAX(ERR,MAXERR)

            RHOTEMP(1:M,1:M) = E(NX:NX+M-1,NY:NY+M-1,K)
            CALL POLIN2(RATEITEMP,RATEJTEMP,RHOTEMP,M,M, &
                RATEINEW(I),RATEJNEW(J),ENEW(I,J,K),ERR)
            MAXERR = MAX(ERR,MAXERR)

          END DO
        END DO
      END DO
      WRITE(*,*) "POLINT MAXERR=", MAXERR

      ! New mov data output
      DO MY_ID = 0,(INODENEW*JNODENEW-1)
        DO I = 1,INODENEW
          DO J = 1,JNODENEW
            ITEMP_ID = (I-1)+INODENEW*(J-1)
            IF (MY_ID.EQ.ITEMP_ID) THEN
              NDEC_ID = I-1 ! Column numbering
              JDEC_ID = J-1 ! Depth numbering
              GO TO 2020
            END IF
          END DO
        END DO

2020    CONTINUE

!        WRITE(*,*) 'WORKING ON:', MY_ID

        IF (INODENEW.EQ.1) THEN
!          IL1 = 1
          IL = IPLNEW
        ELSE IF (NDEC_ID.EQ.0) THEN
          ! Index limits along inlet
!          IL1 = 1
          IL = IBC0 + INCALNEW + IETR
        ELSE IF (NDEC_ID.EQ.INODENEW-1) THEN
          ! Index limits along exit
!          IL1 = 1
          IL = IETR + INCALNEW + IBC1
        ELSE
          ! Index limits for interior columns
!          IL1 = 1
          IL = IETR + INCALNEW + IETR
        END IF

        IF (JNODENEW.EQ.1) THEN
          ! Index limits for a single row zone
!          JL1 = 1
          JL  = JPLNEW
        ELSE IF (JDEC_ID.EQ.0) THEN
          ! Index limits for surface row
!          JL1 = 1
          JL  = JETR + JNCALNEW + JDIFF
        ELSE IF (JDEC_ID.EQ.JNODENEW-1) THEN
          ! Index limits for shock row
!          JL1 = 1
          JL  = JETR + JNCALNEW - JDIFF + JBC1
        ELSE
          ! Index limits for interior rows
!          JL1 = 1
          JL  = JETR + JNCALNEW + JETR
        END IF

        ! Set global index for end points of each MPI zone
        IF (NDEC_ID.eq.0) THEN
          ! Inlet
          IL1_G = 1
!          IL_G  = IL1_G + IL - 1
        ELSE
          ! Interior columns and exit
          IL1_G = IBC0 + NDEC_ID*INCALNEW - IETR + 1
!          IL_G  = IL1_G + IL - 1
        END IF

        IF (JDEC_ID.EQ.0) THEN
          ! Surface 
          JL1_G = 1
        ELSE
          ! Interior rows and shock 
          JL1_G = JDEC_ID*JNCALNEW - JETR + JDIFF + 1
        END IF

!        WRITE(*,*) 'Global Index positions',IL1_G, JL1_G, IL, JL

        ! Write mov files
        II = IL1_G-1
        JJ = JL1_G-1
        KL = KPL
        WRITE(FILEN,'("movnew/mov-",I3.3,".dat")') MY_ID
        OPEN(11,FILE=TRIM(FILEN),FORM='unformatted')
        WRITE(11) NADV, CCFL, TIME, IL, JL, KL
        WRITE(11) (((RHOSNEW(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        WRITE(11) (((RHOUNEW(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        WRITE(11) (((RHOVNEW(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        WRITE(11) (((RHOWNEW(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        WRITE(11) (((   ENEW(II+I,JJ+J,K),I=1,IL),J=1,JL),K=1,KL)
        WRITE(11) (( HSHKNEW(II+I,K),I=1,IL),K=1,KL)
        WRITE(11) (( DSTNEW(II+I,K),I=1,IL),K=1,KL)
        CLOSE(11)

      END DO

      WRITE(*,*) "Remesh Complete"

      END PROGRAM REGRID



      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
!
!   General polynomial interpolation of a function. The code is from 
!   the book: Numerical Recipes.
!
      INTEGER(4) N, NMAX
      REAL(8) DY, X, Y, XA(N), YA(N)
      PARAMETER (NMAX=10)
      INTEGER(4) I, M, NS
      REAL(8) DEN, DIF, DIFT, HO, HP, W, C(NMAX), D(NMAX)
      NS=1
      DIF=DABS(X-XA(1))
      DO I=1,N
         DIFT=DABS(X-XA(I))
         IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
         END IF
         C(I)=YA(I)
         D(I)=YA(I)
      END DO
      Y=YA(NS)
      NS=NS-1
      DO M=1,N-1
         DO I=1,N-M
            HO=XA(I)-X
            HP=XA(I+M)-X
            W=C(I+1)-D(I)
            DEN=HO-HP
!            IF (DEN.EQ.0.0D0) WRITE(*,*) 'ERROR: SUB POLINT'
            DEN=W/DEN
            D(I)=HP*DEN
            C(I)=HO*DEN
         END DO
         IF (2*NS.LT.N-M)THEN
            DY=C(NS+1)
         ELSE
            DY=D(NS)
            NS=NS-1
         END IF
         Y=Y+DY
      END DO
      RETURN
      END SUBROUTINE POLINT



      SUBROUTINE POLIN2(X1A,X2A,YA,M,N,X1,X2,Y,DY1)
      INTEGER(4) M,N,NMAX,MMAX
      REAL*8 DY,DY1,X1,X2,Y,X1A(M),X2A(N),YA(M,N)
      PARAMETER (NMAX=20,MMAX=20)
!     USES POLINT
      INTEGER(4) J,K
      REAL*8 YMTMP(MMAX),YNTMP(NMAX)
      DY1 = 0.0D0
      DO J=1,M
        DO K=1,N
          YNTMP(K)=YA(J,K)
        END DO
        CALL POLINT(X2A,YNTMP,N,X2,YMTMP(J),DY)
        DY1 = MAX(DY,DY1)
      END DO
      CALL POLINT(X1A,YMTMP,M,X1,Y,DY)
      DY1 = MAX(DY,DY1)
      RETURN
      END SUBROUTINE POLIN2



      SUBROUTINE HUNT(xx,n,x,jlo)  
      INTEGER(4) jlo,n  
      REAL(8) x,xx(n)  
      INTEGER(4) inc,jhi,jm  
      LOGICAL ascnd  
      ascnd=xx(n).GE.xx(1)  
      IF (jlo.LE.0.OR.jlo.GT.n) THEN
        jlo=0  
        jhi=n+1  
        GOTO 3  
      END IF
      inc=1  
      IF (x.GE.xx(jlo).EQV.ascnd) THEN
1       jhi=jlo+inc  
        IF (jhi.GT.n) THEN
          jhi=n+1  
        ELSE IF (x.GE.xx(jhi).EQV.ascnd) THEN
          jlo=jhi  
          inc=inc+inc  
          GOTO 1  
        END IF
      ELSE  
        jhi=jlo  
2       jlo=jhi-inc  
        IF (jlo.LT.1) THEN
          jlo=0  
        ELSE IF (x.LT.xx(jlo).EQV.ascnd) THEN
          jhi=jlo  
          inc=inc+inc  
          GOTO 2  
        END IF   
      END IF
3     IF (jhi-jlo.EQ.1) THEN
        IF (DABS(x-xx(n)).LT.1.0D-10) jlo=n-1
        IF (DABS(x-xx(1)).LT.1.0D-10) jlo=1
!        IF (x.EQ.xx(n)) jlo=n-1  
!        IF (x.EQ.xx(1)) jlo=1  
        RETURN  
      END IF  
      jm=(jhi+jlo)/2  
      IF(x.GE.xx(jm).EQV.ascnd)THEN
        jlo=jm  
      ELSE  
        jhi=jm  
      END IF   
      GOTO 3  
      END SUBROUTINE HUNT
