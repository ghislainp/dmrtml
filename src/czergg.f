C	   SUBROUTINE CZERGG(ZEROS,RHO,RHOF,N,ALPHA,M,NMAX,B,B2,D)
C      *******************************************************
C BUT: CALCULER LES NOEUDS ET LES POIDS DE LA FORMULE GENERALISEE DE GAUSS
C ------------------------------------------------------------
C TRANSLATED FROM PASCAL PROCEDURE CZEROS OF C.BERNARDI
C BY O.KOUTCHMY      AUGUST 1990
C ------------------------------------------------------------
C PARAMETRES D'ENTREE:
C -------------------
C N      -NOMBRE DE NOEUDS INTERIEURS
C ALPHA  -TYPE DE POLYNOME DE JACOBI:
C         0 POUR LES POLYNOMES DE LEGENDRE
C M      -TYPE DE LA FORMULE DE GAUSS GENERALISEE:
C         1 POUR LA FORMULE DE GAUSS-LOBATTO
C NMAX   -NOMBRE DE NOEUDS INTERIEURS + 1
C B,B2,D -TABLEAUX DE TRAVAIL
C
C PARAMETRES DE SORTIE:
C --------------------
C ZEROS -NOEUDS DE LA FORMULE DE GAUSS GENERALISEE
C RHO   -POIDS ASSOCIES AUX NOEUDS INTERIEURS
C RHOF  -POIDS ASSOCIES AUX VALEURS DE LA FONCTION ET DE SES DERIVEES 1,...,M-1
C        AUX POINTS -1 ET +1
CC      *******************************************************
C       SUBROUTINE CZERGG(ZEROS,RHO,RHOF,N,ALPHA,M,NMAX,B,B2,D)
CC      *******************************************************
C       IMPLICIT REAL*8 (A-H,O-Z)
C       DIMENSION ZEROS(NMAX),RHO(NMAX),RHOF(0:M)
C       DIMENSION B(NMAX),B2(NMAX),D(NMAX)
CC
C       ALPHAM=ALPHA+M
C       IF(ALPHA .GT. -1.) THEN
C       CALL CACOEF(N,ALPHAM,B,B2,NMAX)
C       CALL CALZER(N,ALPHAM,ZEROS,B2,D,NMAX)
C       CALL PGAUSS(N,ALPHAM,B,ZEROS,RHO,NMAX)
C       IF(M .GE. 1) THEN
C       CALL PFRONT(N,M,ALPHA,RHOF,Z)
C       DO 10 I=1,N
C       RHO(I)=RHO(I)/DEXP(M*DLOG(1.D0-ZEROS(I)**2))
C10     CONTINUE
C       END IF
CC
C       END IF
C       END

C     function ajoutée par Ghislain pour le cas gauss simple 
C     avec polynome de Legendre
C      *******************************************************
       SUBROUTINE CZERG(ZEROS,RHO,N,NMAX,B,B2,D)
C      *******************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION ZEROS(NMAX),RHO(NMAX)
       DIMENSION B(NMAX),B2(NMAX),D(NMAX)
C
       ALPHAM=0.0
       CALL CACOEF(N,ALPHAM,B,B2,NMAX)
       CALL CALZER(N,ALPHAM,ZEROS,B2,D,NMAX)
       CALL PGAUSS(N,ALPHAM,B,ZEROS,RHO,NMAX)
C
       END



C          ************************************
C          SUBROUTINE ITER(N,S,TEST,D,B2,NMAX)
C          ************************************
C N IS THE NUMBER OF INTERIOR NODES
C NMAX=N+1
C ITER IS CALLED BY CALZER
C DATA: N,S,D,B2,NMAX
C RESULTS : TEST,D,B2
C -------------------------------------------------
C TRANSLATED FROM C.BERNARDI PASCAL PROCEDURE ITERE
C BY OLGA KOUTCHMY
C JULY 1990
C          ************************************
           SUBROUTINE ITER(N,S,TEST,D,B2,NMAX)
C          ************************************
	   IMPLICIT REAL*8 (A-H,O-Z)
	   DIMENSION D(NMAX),B2(NMAX)
C
	   Q=D(1)-S
	   E=B2(1) / Q
	   QQ= Q + E
	   B2(N) = 0.D0
	   P = 1.D0
	   DO 10 I = 2, N
	   Q = D(I) - S - E
	   X = Q / QQ
	   EE = E * X
	   D(I-1) = QQ + EE
	   E = B2(I) / Q
	   QQ = Q - EE + E
	   B2(I-1) = EE* QQ
	   P = P * X + 1.D0
10         CONTINUE
	   D(N) = QQ
	   S = QQ / P
	   TEST = Q - EE
	   END
C          ************************************
C	   SUBROUTINE CACOEF(N,ALPHA,B,B2,NMAX)
C          ************************************
C TRANSLATED FROM C.BERNARDI PASCAL PROCEDURE CALCULCOEF
C BY OLGA KOUTCHMY
C AUGUST 1990
C ---------------------------------------------------------
C NMAX=N+1
C          ************************************
	   SUBROUTINE CACOEF(N,ALPHA,B,B2,NMAX)
C          ************************************
	   IMPLICIT REAL*8 (A-H,O-Z)
	   DIMENSION B(NMAX),B2(NMAX)
C
	   IF( ALPHA .EQ. -0.5 ) THEN
	   B2(1) = 0.5D0
	   ELSE
	   B2(1) = ( 1.D0 + 2.D0 * ALPHA ) /
     S             ( 4.D0 * ( ALPHA + 1.D0 ) **2 - 1.D0 )
	   END IF
	   B(1) = DSQRT( B2(1) )
	   DO 10 I = 2, N - 1
	   B2(I) = DBLE(I) * ( DBLE(I) + 2.D0 * ALPHA ) /
     S             ( 4.D0 * ( DBLE(I) + ALPHA ) ** 2 - 1.D0 )
	   B(I) = DSQRT( B2(I) )
10         CONTINUE
	   END
C          **************************************
C          SUBROUTINE CALZER(N,ALPHA,Z,B2,D,NMAX)
C          **************************************
C TRANSLATED FROM C.BERNARDI PASCAL PROCEDURE CALCULZEROS
C BY OLGA KOUTCHMY
C AUGUST 1990
C ----------------------------------------------------------
C NMAX=N+1
C          **************************************
           SUBROUTINE CALZER(N,ALPHA,Z,B2,D,NMAX)
C          **************************************
	   IMPLICIT REAL*8 (A-H,O-Z)
	   DIMENSION Z(NMAX),B2(NMAX),D(NMAX)
C
	   EPS = 1.D-15
	   DO 10 I = 1, N
	   D(I) = 0.D0
10         CONTINUE
	   S = -1.D0
	   T = 0.D0
	   DO 20 M = N, 2, -1
15         CONTINUE
	   T = T + S
	   CALL ITER( M, S, TEST, D, B2, NMAX)
	   IF( TEST .GE. EPS ) GOTO 15
	   Z( N + 1 - M ) = D(M) + T
20         CONTINUE
	   Z(N) = D(1) + T
	   END
C          ***************************************
C          SUBROUTINE PGAUSS(N,ALPHA,B,Z,RHO,NMAX)
C          ***************************************
C
C TRANSLATED FROM C.BERNARDI PASCAL PROCEDURE CALCULPOIDSGAUSS
C BY O.KOUTCHMY    AUGUST 1990
C REMARK. MO IS RENAMED RMO, LA IS RENAMED RLA
C
C          ***************************************
           SUBROUTINE PGAUSS(N,ALPHA,B,Z,RHO,NMAX)
C          ***************************************
	   IMPLICIT REAL*8 (A-H,O-Z)
	   DIMENSION B(NMAX),Z(NMAX),RHO(NMAX)
C
	   PI=3.141592653589793239D0
	   ALPHA1=ALPHA+1.D0
	   ALPHA2=ALPHA+1.5D0
	   RMO=DSQRT(PI)*GAMMAR(ALPHA1)/GAMMAR(ALPHA2)
	   DO 10 I=1,N
	   RLA=Z(I)
	   X=1.D0
	   Y=RLA/B(1)
	   S=1.D0+Y*Y
	   DO 5 J=3,N
	   T=(RLA*Y-B(J-2)*X)/B(J-1)
	   X=Y
	   Y=T
	   S=S+Y*Y
5          CONTINUE 
           RHO(I)=RMO/S
10         CONTINUE 
           END
C          ***********************************
C          SUBROUTINE PFRONT(N,M,ALPHA,RHOF,Z)
C          ***********************************
C
CC TRANSLATED FROM C.BERNARDI PASCAL PROCEDURE CALCULPOIDSFRONTIERE
CC BY O.KOUTCHMY    AUGUST 1990
CC WARNING ! RECURSIVE CALL                       
CC
CC          ***********************************
C           SUBROUTINE PFRONT(N,M,ALPHA,RHOF,Z)
CC          ***********************************
C	   IMPLICIT REAL*8 (A-H,O-Z)
C	   DIMENSION RHOF(0:M)
CC
C	   IF( M .EQ. 1 ) THEN
CC
C	   RHOF(0)=DEXP((2.D0*ALPHA+1.D0)*DLOG(2.D0))*(ALPHA+1.D0)*
C     S             GAMMAR(ALPHA+1.D0)**2
C     S            *GAMMAR(N+1.D0)/GAMMAR(N+2.D0*ALPHA+3.D0)
C           Z=RHOF(0)
CC
C	   ELSE               
CC
C	   CALL PFRONT(N,M-1,ALPHA+1.D0,RHOF,Z)
C	   RHOF(M)=0.D0
C	   SOMME=0.D0 
C	   DO 5 K=M-1,1,-1
C	   RHOF(K)=-0.5D0*(K+1.D0)*RHOF(K+1)-0.5D0*RHOF(K-1)/DBLE(K)
CC	   IF(ODD(N)) THEN
C           IF(MOD(N,2) .EQ. 1) THEN
C	   SOMME=(SOMME+RHOF(K)+DBLE(K+1)*RHOF(K+1))*DBLE(N+1-K)*
C     S       (DBLE(N)+2.D0*ALPHA+DBLE(2*M+K))/(2.D0*(ALPHA+DBLE(M+K)))
C	   ELSE
C	   SOMME=(SOMME+RHOF(K))*DBLE(N+1-K)*
C     S       (DBLE(N)+2.D0*ALPHA+DBLE(2*M+K))/(2.D0*(ALPHA+DBLE(M+K)))
C	   END IF
C5          CONTINUE
CC 
CC          IF(ODD(N)) THEN
C           IF(MOD(N,2) .EQ. 1) THEN
C	   Z=Z*DBLE(N+2*M-3)*(DBLE(N)+2.D0*ALPHA+4.D0)/
C     S       (4.D0*(ALPHA+1.D0)*DBLE(M-1))
C	   ELSE
C	   Z=Z*DBLE(N+2*M-2)*(DBLE(N)+2.D0*ALPHA+3.D0)/
C     S       (4.D0*(ALPHA+1.D0)*DBLE(M-1))
C	   END IF
CC 
C           RHOF(0)=Z-SOMME
CC 
CC	   IF(ODD(N)) THEN
C           IF(MOD(N,2) .EQ. 1) THEN
C	   RHOF(0)=RHOF(0)-RHOF(1)
C	   END IF
CC 
C	   END IF
C	   END
C      **************************
C       REAL*8 FUNCTION GAMMAR (X )
C      **************************
C
C TRANSLATED FROM PASCAL FUNCTION GAMMA OF C.BERNARDI
C BY O.KOUTCHMY      AUGUST 1990
C WARNING ! RECURSIVE CALL  
C
C      **************************
       REAL*8 FUNCTION GAMMAR (X )
C      **************************
       IMPLICIT REAL*8 (A-H,O-Z)
C
       C2  = +0.5772156649015328606D0
       C3  = -0.655878071520253881D0
       C4  = -0.0420026350340952D0
       C5  = +0.1665386113822915D0
       C6  = -0.0421977345555443D0
       C7  = -0.0096219715278770D0
       C8  = +0.0072189432466630D0
       C9  = -0.0011651675918591D0
       C10 = -0.0002152416741149D0
       C11 = +0.0001280502823882D0
       C12 = -0.0000201348547807D0
       C13 = -0.0000012504934821D0
       C14 = +0.0000011330272320D0
       C15 = -0.0000002056338417D0
       C16 = +0.0000000061160950D0
       C17 = +0.0000000050020075D0
       C18 = -0.0000000011812746D0
       C19 = +0.0000000001043427D0
C
       GAMMAR=1.0D0
       DO WHILE ( (DABS(X) .GT. 0.5).AND.(X .NE. 1.0D0) )
C    
          X = X - 1.D0
          GAMMAR = X * GAMMAR
       END DO

       IF( X .EQ. 1.) THEN
          ! rien a faire
       ELSE

       Y = C17 + X * (C18 + X * C19)
       Y = C9 + X * (C10 + X * (C11 + X * 
     S      (C12 + X * (C13 + X * (C14 + X *
     S      (C15 + X * (C16 + X * Y)))))))
       Y = X * (1.D0 + X * (C2 + X * (C3 + X * 
     S      (C4 + X * (C5 + X * (C6 + X *
     S      (C7 + X * (C8 + X * Y))))))))
       GAMMAR = GAMMAR * 1.D0/ Y
       END IF
       END
