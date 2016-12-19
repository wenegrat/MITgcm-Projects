C $Header: /u/gcmpack/MITgcm/pkg/diagnostics/diagnostics_fill_state.F,v 1.30 2006/02/16 23:03:16 molod Exp $
C $Name:  $

#include "DIAG_OPTIONS.h"

CBOP
C     !ROUTINE: DIAGNOSTICS_FILL_STATE
C     !INTERFACE:
      SUBROUTINE ADV_DIAGNOSTICS_FILL(fldType, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DIAGNOSTICS_FILL_STATE
C     | o Fill-in main code, state-variables diagnostics
C     *==========================================================*
C     \ev
C     ASSUMPTIONS: f-plane, Linear EOS function of T only, 
C      
C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "NH_VARS.h"
#include "SURFACE.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid     :: my Thread Id number
      INTEGER fldType
      INTEGER myThid 
C     aField     :: field to advect
#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     == Local variables ==
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
      _RL advOut  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL aFld    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL uFld    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vFld    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL wFld    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL dFdx    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL dFdy    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL dFdz    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL adFdx   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL adFdy   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL adFdz   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
  
      INTEGER i,j,K,bi,bj,jj,ii
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      

      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO K=1,Nr
C     Make local copies of flow field and zero variables
              DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
                advOut(i,j,K,bi,bj) = 0. _d 0
                dFdx(i,j,K)   = 0. _d 0
                dFdy(i,j,K)   = 0. _d 0
                dFdz(i,j,K)   = 0. _d 0
                uFld(i,j,K) = uVel(i,j,K,bi,bj)
                vFld(i,j,K) = vVel(i,j,K,bi,bj)
                wFld(i,j,K) = wVel(i,j,K,bi,bj)
              ENDDO
              ENDDO
            ENDDO


C     Select field and Interpolate to rho points
      IF (fldType.EQ.0) THEN
                aFld = gravity*tAlpha*theta
      ENDIF
C     U 
      IF (fldType.EQ.1) THEN 
              DO k=1,Nr
              DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx-1
                aFld(i,j,k,bi,bj) = 0.5*(uVel(i+1,j,k,bi,bj)
     &                                   +uVel(i,j,k,bi,bj))
              ENDDO
              ENDDO
              ENDDO
      ENDIF 
C     V Grid
      IF (fldType.EQ.2) THEN 
              DO k=1,Nr
              DO j=1-OLy,sNy+OLy-1
              DO i=1-OLx,sNx+OLx
                aFld(i,j,k,bi,bj) = 0.5*(vVel(i,j+1,k,bi,bj)
     &                                   +vVel(i,j,k,bi,bj))
              ENDDO
              ENDDO
              ENDDO
      ENDIF 


C  Calculate X Derivatives
            DO j=1-OLy,sNy+OLy
            DO k=1,Nr
             DO i=1-OLx,sNx+OLx-1
               dFdx(i+1,j,k) = recip_dxC(i+1,j,bi,bj)*
     &                            (aFld(i+1,j,k,bi,bj)-
     &                             aFld(i,j,k,bi,bj))
             ENDDO
             DO i=1-OLx,sNx+OLx-1
               adFdx(i,j,k,bi,bj) = 0.5*(uFld(i,j,k)*dFdx(i,j,k) +
     &                               uFld(i+1,j,k)*dFdx(i+1,j,k))
             ENDDO
            ENDDO
            ENDDO

C  Calculate Y Derivatives
            DO i=1-OLx,sNx+OLx
            DO k=1,Nr
             DO j=1-OLy,sNy+OLy-1
               dFdy(i,j+1,k) = recip_dyC(i,j+1,bi,bj)*
     &                            (aFld(i,j+1,k,bi,bj)-
     &                             aFld(i,j,k,bi,bj))
             ENDDO
             DO j=1-OLy,sNy+OLy-1
             adFdy(i,j,k,bi,bj) = 0.5*(vFld(i,j,k)*dFdy(i,j,k) +
     &                               vFld(i,j+1,k)*dFdy(i,j+1,k))
             ENDDO
            ENDDO
            ENDDO

C  Calculate Vertical Derivatives
            DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             DO k=1,Nr-1
               dFdz(i,j,k+1) = recip_drC(k)*
     &                            (aFld(i,j,k,bi,bj)-
     &                             aFld(i,j,k+1,bi,bj))
             ENDDO
             DO k=1,Nr-1
               adFdz(i,j,k,bi,bj) = 0.5*(wFld(i,j,k)*dFdz(i,j,k) +
     &                               wFld(i,j,k+1)*dFdz(i,j,k+1))
             ENDDO
            ENDDO
            ENDDO


C Enddo: Bi, Bj
        ENDDO
      ENDDO
C  Calculate advection, defined as RHS term.
    
      advOut = adFdx + adFdy + adFdz
      IF (fldType.EQ.0) THEN
      CALL DIAGNOSTICS_FILL(advOut,'BADVTERM',0,Nr,0,1,1,myThid)
      ENDIF
      IF (fldType.EQ.1) THEN
      CALL DIAGNOSTICS_FILL(advOut,'UADVTERM',0,Nr,0,1,1,myThid)
      ENDIF
      IF (fldType.EQ.2) THEN
      CALL DIAGNOSTICS_FILL(advOut,'VADVTERM',0,Nr,0,1,1,myThid)
      ENDIF 
             

#endif /* ALLOW_DIAGNOSTICS */
      
      RETURN 
      END
