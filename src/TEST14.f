
      SUBROUTINE TEST14
C-----------------------------------------------------------------
C     THE TOTAL OXYGEN PARTICLE FLUX AS A FUNCTION OF NORMALIZED
C     COLLISIONALITY
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR,L
      REAL M,T,DEN,DS,CFF4,XI,TAU, EPS, SIGMA, NORM, RESUL
      REAL RHO, RN, E, Q, BN,  OMTH, ZSP
      LOGICAL NEOGEO, NEOFRC

      PARAMETER(NAR = 5)
      PARAMETER(NZM = 30)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF4(NAR,NZM,4),XI(NAR,NZM),
     +          TAU(NAR,NAR),SIGMA(4),RESUL(40,6)


C     SET THE PARAMETERS FOR THE CIRCULAR GEOMETRY
      RHO = 0.05
      E = 0.2
      Q = 2
      RN = 1.65
      BN = 2.5
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
C     USE THE CIRCULAR GEOMETRY APPROXIMATION
      ISEL = 2
C     SET THE ACCURACY
      EPS = 1E-5
C     ALL REGIMES IN THE CALCULATION OF THE VISCOSITY
      NREG = 0
c     THE NUMBER OF LEGENDRE HARMONICS
      NLEG = 3
C     USE ENERGY SCATTERING IN THE CALC. OF VISCOSITY
      NENERGY = 1
C     ION-ELECTRON COLLISIONS
      NCOF = 1
C     IN THE FIRST CALL THE MATRICES HAVE TO BE CALCULATED
      NEOFRC = .FALSE.
C     RECALCULATE THE GEOMETRY DEPENDENT PARAMETERS EVERY TIME
      NEOGEO = .TRUE.

C     CALCULATE THE TOTAL FLUX 
      IC = 3

C     THE NUMBER OF SPECIES IS 2
      NS = 2
C     IONS AND ELECTRONS HAVE ONLY ONE CHARGE 
      NC(1) = 1
      NC(2) = 1
C     THE MASS OF THE HYDROGEN AND OXYGEN
      M(1) = 1.6727E-27
      M(2) = 16*1.6727E-27
C     THE CHARGE OF THE HYDROGEN AND OXYGEN
      ZSP(1,1) = 1
      ZSP(2,1) = 4
C     THE DENSITY OF THE SPECIES IN 10^19 M^-3
      DEN(1,1) = 1.
      DEN(2,1) = 0.03 
C     THE TEMPERATURE IN KEV
C     SET THE COUPLING
      SIGMA(1) = 0
      SIGMA(2) = 0
      SIGMA(3) = 0 
      SIGMA(4) = 1

C     LOOP nu star
      DO 10000 L = 1, 40

      T(1) = 1.5*EXP(-(L-1)/20. * 1.5 * LOG(10.))
      T(2) = T(1)


      DO 304 I = 1, NS
        DO 304 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 304  CONTINUE
      DS(1,1,1) = 1.


      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,CFF4)


      CALL COLXI(NAR,NZM,NS,NC,ZSP,DEN,T,M,TAU,XI)

      OMTH = SQRT(2*1.6E-16*T(1)*M(1))*DEN(1,1)*1E19/TAU(1,1)
     +       /(Q*RN)

      NORM = 1.6E-22*BN**2*ZSP(2,1)/(2*(Q)**2*T(2)*TAU(2,2))

C     USE INVERSE ASPECT RATIO AS X AXIS
      RESUL(L,1) = 1/(SQRT(2.)*OMTH)
      RESUL(L,2) = CFF4(2,1,1)*NORM*RESUL(L,1)
      RESUL(L,1) = RESUL(L,1) / SQRT(E**3)


10000 CONTINUE

      OPEN(11, FILE = 'TEST14c.DAT')
C     THE COLLUMS HAVE THE FOLLOWING MEANING
C     COLLUM 1. : THE NORMALIZED COLLISIONALITY
C     COLLUM 2. : THE OXYGEN PARTICLE FLUX
      DO 30000 L = 1, 40
        WRITE(11,30001)(RESUL(L,I),I = 1,2)
30000 CONTINUE
30001 FORMAT(9(1X,1PE13.5))

      RETURN
      END
