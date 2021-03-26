#include "cppdefs.h"
      MODULE sst_nudge_mod
#ifdef COIN_SST 
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine implement a combined optimal interpolation and nudging !
!  scheme for assimilating OSISAF sea ice concentration into ROMS.     !
!                                                                      !
!  reference: 
!  Wang, K., J. Debernard, A. Sperrevik, P. Isachsen, T. Lavergne:     !
!    A combined optimal interpolation and nudging scheme to assimilate !
!    OSISAF sea ice concentration into ROMS, Annals of Glaciology,     !
!    DOI: https://doi.org/10.3189/2013AoG62A138                        !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC sst_nudge

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sst_nudge (ng, tile)
!***********************************************************************

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
      USE mod_forces

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 84)
# endif
!
      CALL sst_nudge_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
     &                     OCEAN(ng) % t,                               &
     &                     FORCES(ng) % sst,                            &
     &                     FORCES(ng) % sste)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 84)
# endif
      RETURN
      END SUBROUTINE sst_nudge
!
!***********************************************************************
      subroutine sst_nudge_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs, t,                               &
     &                           sst, sste)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_parallel
      USE mod_iounits

      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile

# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: sst(LBi:,LBj:)
      real(r8), intent(inout) :: sste(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: sst(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: sste(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
! Local variable definitions
!
      integer :: i, j, itrc
      logical :: nudge_now = .FALSE.

      real(r8) :: rda, mod_err, mod_err2, gain, wt
      real(r8), parameter :: frln = -0.0543_r8

# include "set_bounds.h"

! Deside if we will do nudging now or not 
      nudge_now = assim_SST(ng)

      IF (.not. nudge_now) RETURN

      rda = dt(ng) * sec2day / SSTnudass(ng)
      IF (Master) THEN
         WRITE(stdout,*) 'assim_SST = ', assim_SST(ng)
         WRITE(stdout,*) 'p_nudging  = ', p_nudging(ng)
         WRITE(stdout,*) 'SSTnudass = ', SSTnudass(ng)
         WRITE(stdout,*) 'dt = ', dt(ng)
         WRITE(stdout,*) 'rda = ', rda
      ENDIF

      DO j=JstrR,JendR
        DO i=IstrR,IendR

           IF ((sst(i,j) >= -1.9_r8) .AND. (sst(i,j) <= 50.0_r8)) THEN
              mod_err  = t(i,j,N(ng),nrhs,itemp) - sst(i,j)
              mod_err2 = mod_err**2 + sste(i,j)**2

              gain   = mod_err2 / (mod_err2 + sste(i,j)**2 + 1.0e-6_r8)
              IF (p_nudging(ng)) gain = 0.9999_r8
              wt = 1.0_r8 - (1.0_r8 - gain)**rda

              t(i,j,N(ng),nrhs,itemp) = t(i,j,N(ng),nrhs,itemp) &
     &                                - wt * mod_err
!              sfrz = max(0.0, t(i,j,N(ng),nrhs,itemp) / frln)
!              t(i,j,N(ng),nrhs,isalt) = max(sfrz,t(i,j,N(ng),nrhs,isalt))
           ENDIF

        ENDDO
      ENDDO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
           CALL exchange_r3d_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nrhs,itrc))
        END DO
      ENDIF

#  ifdef DISTRIBUTE   
      DO itrc=1,NT(ng)
        CALL mp_exchange3d (ng, tile, iNLM, 1,                          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    t(:,:,:,nrhs,itrc))
      END DO
#  endif

      END SUBROUTINE sst_nudge_tile

#endif
      END MODULE sst_nudge_mod

