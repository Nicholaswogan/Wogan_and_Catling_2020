      SUBROUTINE H2aq
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      CHARACTER*8 ISPEC,REACTYPE,PLANET,CHEMJ
      CHARACTER*20 fmtstr,fmtdatastr
      CHARACTER*60 fmtheadstr

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/FBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'  !this has flux in it - compare against below
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/WBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ZBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SATBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SULBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RRATS.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/MBLOK.inc'

      COMMON/LifeTime/TAUO2,TAUCH4,TAUSO2

      DIMENSION FUP(NQ1),FLOW(NQ1),CON(NQ1),FLUXO(NQ1,NZ)
     2  ,ZF(NZ),TAUAER(NP+1)
      dimension TAUHCABS(kw),TAUHCEXT(kw), TAUHCEFF(kw)
      dimension vdep(nq), flux_H(nz),USOLORIG(NQ,NZ)
      dimension PSO2MIF(nz),PROSO2MIF(nz),PROSO(nz)


C
C-PK  This subroutine calculates the concentration of dissolved H2
C-PK  based on free energy considerations applied to the reaction
C-PK      CO2 + 4H2 --> CH4 + 2H2O (H2-based methanogenesis).
C-PK  It follows Kral et al.,OLEB,1998 and Kasting et al.,OLEB,2001.
C-PK  It uses the following Gibbs free energy relation:
C-PK      DG = DG0 + RTlnQ ("precursor" to Nernst equation)
C-PK  Used T = 25C for now, but have incorporated temp dependence in DG0
      T = 298.
C-PK  Use DG = 30 since that is energy needed to make 1 ATP
      DELTA_G = -30.
      DELTA_G0 = -253. + 0.41*298
      R = 0.008314
      Q = EXP((DELTA_G - DELTA_G0)/(R*298))
C-PK      PRINT *,'Q=',Q

      ALPHA_CO2 = 3.4E-2
      CO2_aq = ALPHA_CO2*USOL(LCO2,1)

      ALPHA_CH4 = 1.4E-3
      SCALE = 6.02E20
      DIFF_CH4 = 1.8E-5
      Z_FILM = 4.E-3
      EQUIL_CH4 = ALPHA_CH4*USOL(LCH4,1)
      VPIS_CH4 = DIFF_CH4/Z_FILM
C-PK  Calculate [CH4]aq assuming piston velocity formulation for CH4 flux:
C-PK     PHI_CH4 = VPIS_CH4 * SCALE * ([CH4]aq - ALPHA_CH4*pCH4)
      CH4_aq = FLOW(LCH4)/(VPIS_CH4*SCALE) + EQUIL_CH4

      ALPHA_H2 = 7.8E-4
      EQUIL_H2 = ALPHA_H2*USOL(LH2,1)
C-PK  [H2]aq calculated from Q = [CH4]aq*aCO2/([CO2]aq*aCH4 * ([H2]aq*aH2)^4)
      H2_aq = ALPHA_H2*(CH4_aq*ALPHA_CO2/(CO2_aq*ALPHA_CH4*Q))**0.25

      PRINT 100
 100  FORMAT(1X,'DISSOLVED CH4 AND H2 (ECO MODEL)'/)
      CH4_H2_RATIO = USOL(LCH4,1)/USOL(LH2,1)
C-PK  See main program (deposition velocities section) for derivation of vmax(H2)
      VMAX_H2 = 2.4E-4
      VECO_H2 = VMAX_H2*(EQUIL_H2 - H2_aq)/EQUIL_H2
      VPHOTO_H2 = FLOW(LH2)/(USOL(LH2,1)*DEN(1))
      VRATIO = VPHOTO_H2/VECO_H2
      PRINT 110, VPIS_CH4,CH4_aq,EQUIL_CH4,H2_aq,EQUIL_H2,CH4_H2_RATIO,
     2           VMAX_H2,VECO_H2,VRATIO
 110  FORMAT(5X,'VPIS(CH4) =',1PE10.3,3X,'[CH4]aq =',1PE11.4,
     2  3X,'alphaCH4*pCH4 =',1PE11.4,3X,'[H2]aq =',1PE11.4,
     3  3X,'alphaH2*pH2 =',1PE11.4,//5X,'f(CH4)/f(H2) =',1PE11.4,
     4  3X,'VMAX(H2) =',1PE11.4,3X,'VECO(H2) =',1PE11.4,
     5  3X,'VPHOTO(H2):VECO(H2) =',1PE11.4/)

      RETURN
      END
