! ***********************************************************************
!
!   Copyright (C) 2011  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

!===================================================================================
!                      WIND SCHEME ID LEGEND
!-----------------------------------------------------------------------------------
! Low-Mass Winds:
!  10.0: R75   (Reimers 1975)
!  11.0: B95   (Bloecker 1995)
!
! Thin Winds (Hot Stars):
!  20.0: GM23  (Gormaz-Matamala+ 2023)
!  21.0: B23   (Bjorklund+ 2023)
!  22.0: V01   (Vink+ 2001)
!  23.0: V17   (Vink 2017)
!  24.0: K24   (Krticka+ 2024)
!  25.0: VS21  (Vink & Sander 2021)
!  26.0: P25   (Pauli+ 2025)
!  27.0: V17   (Vink 2017)
!
! Dust/Cool Winds:
!  30.0: dJ88  (de Jager+ 1988)
!  31.0: vL05  (van Loon+ 2005)
!  32.0: B23   (Beasor+ 2023)
!  33.0: Y23   (Yang+ 2023)
!  34.0: A24   (Antoniadis+ 2024)
!  35.0: D24   (Decin+ 2024)
!
! Thick/WR Winds:
!  40.0: V11   (Vink+ 2011)
!  41.0: B20   (Bestenlehner 2020)
!  42.0: GH08  (Grafener & Hamann 2008)
!  43.0: Y06   (Yoon+ 2006)
!  44.0: NL00  (Nugis & Lamers 2000)
!  45.0: SV20  (Sander & Vink 2020)
!
! Special Cases:
!  90.0: LBV   (Belczynski+ 2010 LBV)
!===================================================================================

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use binary_def
      use utils_lib


      implicit none

      integer :: already_thick = 0                                              ! 0/1/2:    not yet thick/thick with eta/thick with gamma
      real(dp) :: Mdot_switch,L_switch,Mhom_switch,gamma_edd_switch
      real(dp) :: eta,eta_trans,gamma_edd,gamma_edd_old                         ! gamma_edd_old is to check the previous timestep

      real(dp) :: wind_scheme,wind_scheme_interp           ! To know which winds model I am using at each timestep
      ! real(dp) :: wind_scheme,wind_scheme_interp
      !                                                                           ! 0/0.5/1/1.25/1.5/2/2.5/3/4/4.5/5/6: GM23/Bj23/V01/V17/dJ88/A24/K24/B20/V11/SV20/NL00
      real(dp) :: max_years_dt_old
      ! these routines are called by the standard run_star check_model
      contains



      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns


         s% other_wind => my_other_wind
      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         max_years_dt_old = s% max_years_for_timestep

      end subroutine extras_startup


      subroutine my_other_wind(id, L, M, R, Tsurf, X, Y, Z, w, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: L, M, R, Tsurf, X, Y, Z ! surface values (cgs)
        real(dp), intent(out) :: w !wind in units of Msun/year (value is >= 0)
        real(dp) :: gmrstar,gmlogg,lteff,xlmdot,logZ_div_Zsun,hehratio,vterm,Mhom,f,Zsolar
        real(dp) :: Z_div_Z_solar,Teff_jump,alfa,log_gamma_edd,gamma_trans,logL_div_Lsun
        real(dp) :: F1,F2,F3,F4,F5,F6,F7,F8,F9,vesc_eff,vesc,vinf_fac,const_k
        real(dp) :: w1, w2, w_2ndDust
        real(dp) :: divisor,beta,alfa_mid
        logical :: thick_met

        logical :: gamma_condition, eta_condition, xsurf_condition, HPoor_WR_condition

        ! real(dp) :: Xs,Ys,h1,he4
        integer, intent(out) :: ierr
        type(star_info), pointer :: s

        ierr = 0

        call star_ptr(id, s, ierr)

        gmlogg=log10(s%grav(1))
        gmrstar=R/Rsun
        lteff=log10(s% Teff/1000)

        Zsolar = s% x_ctrl(6)

        Z_div_Z_solar = s% kap_rq% Zbase/Zsolar
        logZ_div_Zsun=log10(Z_div_Z_solar)
        logL_div_Lsun=log10(L/Lsun)

        w = 0

        ! -----------------------------------------------

        if ( s% x_character_ctrl(1) /= 'eta' .and. s% x_character_ctrl(5) == 'V11' ) then
          call mesa_error(__FILE__,__LINE__,'V11 winds can only be computed with the eta condition')
        end if

        if ( s% x_logical_ctrl(2) .and. s% x_character_ctrl(5) == 'V11' ) then
          call mesa_error(__FILE__,__LINE__,'V11 winds cannot work with a cool WR wind condition')
        end if

        ! -----------------------------------------------

        if (ierr /= 0) return

        s% max_years_for_timestep = max_years_dt_old                            ! Sometimes I change the timesteps in the code to avoid big shifts at the onset of different mass loss recipes
                                                                                !  This is to set everything back

        if (s% stop_near_zams) then
          write(*,*) 'pre-ZAMS -- no mass loss'
          w=0
          return
        else if (s% kap_rq% Zbase < 1d-40) then                                 ! PopIII
          write(*,*) 'popIII -- no mass loss'
          w=0
          return

        elseif (s% star_mass<8.5) then
          if (((Tsurf < 10**3.7 .and. logL_div_Lsun > 1.0) .or. &               ! GENEC calibrations for low-mass stars + my linear Teff fit
            (log10(Tsurf) < 0.385*logL_div_Lsun+3.705)) .and. &                 ! RGB
            s% center_he4 > s% RGB_to_AGB_wind_switch) then
            call eval_Reimers_wind(w)
            write(*,*) "Here are Reimers RGB cool winds: log(Mdot [Msun/yr]) =", log10(ABS(w))

          elseif (Tsurf < 10**3.7 .and. &                                       ! AGB
            s% center_he4 < s% RGB_to_AGB_wind_switch) then
            call eval_Blocker_wind(w)
            write(*,*) "Here are Blocker AGB cool winds: log(Mdot [Msun/yr]) =", log10(ABS(w))

          elseif ( s% w_div_w_crit_avg_surf > 1.01*s% surf_omega_div_omega_crit_limit ) then
            w=1d-11                                                             ! Mechanical mass loss; just enough to take away angular momentum
            write(*,*) 'Small MS star -- Mechanical mass loss: log(Mdot [Msun/yr]) =', log10(ABS(w))
            s% max_years_for_timestep = 1d3                                     ! I constrain this

          else
            w=0
            write(*,*) 'Small MS star -- no mass loss'

          end if

          return


        end if



        if ( logZ_div_Zsun >= log10(0.2) ) then
          gamma_trans = 0.5
        else
          gamma_trans = 0.5 - 0.301*logZ_div_Zsun-0.045*logZ_div_Zsun**2
        end if

        log_gamma_edd = -4.813d0+log10(1+X)+log10(L/Lsun)-log10(M/Msun)
        gamma_edd=10**log_gamma_edd
        ! gamma_edd = L/s% prev_Ledd

        if ( (12 < M/Msun .and. M/Msun < 250) .and. log10(L/Lsun)<6.5 ) then

          F1 = 4.026
          F2 = 4.277
          F3 = -1.0
          F4 = 25.48
          F5 = 36.93
          F6 = -2.792
          F7 = -3.226
          F8 = -5.317
          F9 = 1.648

        elseif ( M/Msun <= 12 .and. log10(L/Lsun)<6.5 ) then

          F1 = 2.582
          F2 = 0.829
          F3 = -1.0
          F4 = 9.375
          F5 = 0.333
          F6 = 0.543
          F7 = -1.376
          F8 = -0.049
          F9 = 0.036

        elseif ( M/Msun <= 4000 .or. log10(L/Lsun)>=6.5) then
          F1 = 10.05
          F2 = 8.204
          F3 = -1.0
          F4 = 151.7
          F5 = 254.5
          F6 = -11.46
          F7 = -13.16
          F8 = -31.68
          F9 = 2.408

        end if

        const_k = (clight*Msun*1d5)/(Lsun*secyer)

        f = F4 + F5*X + F6*X**2 + (F7+F8*X)*log10(L/Lsun)
        Mhom=10**((F1+F2*X+F3*sqrt(f))/(1+F9*X))


        vesc = sqrt(2d0*standard_cgrav*M/R)/1d5
        vterm = 2.6 * sqrt(2d0*standard_cgrav*(Mhom*Msun)*(1-gamma_edd)/R)/1d5*Z_div_Z_solar**0.20d0

        eta_trans = 0.75/(1+(vesc**2)/(vterm**2))
        ! eta = (ABS(s% mstar_dot /Msun)*secyer * vterm)/(L/(clight))
        eta = const_k*(ABS(s% mstar_dot/Msun*secyer)*vterm)/(L/Lsun)


        write(*,*)
        write(*,*) "Edd_factor_e:", gamma_edd, "vs        Gamma_trans:", gamma_trans
        write(*,*) "previous log(Mdot):", log10(ABS(s% mstar_dot/Msun*secyer))
        write(*,*) "vterm:", vterm, " km/s    vs        vesc:", vesc, " km/s"
        write(*,*) "eta factor:", eta , "vs        eta trans:", eta_trans
        write(*,*) "log(g):", gmlogg
        write(*,*) "Mhom:", Mhom, " Msun"
        write(*,*)

        ! --------------------------------------------- Check for thick winds ---------------------------------------------------
        thick_met = .false.

        if ( s% x_character_ctrl(1) == 'gamma' ) then
          if ( gamma_edd >= gamma_trans .or. already_thick==1 ) then
            thick_met = .true.
          end if
        elseif (s% x_character_ctrl(1) == 'eta') then
          if (eta>=eta_trans .or. already_thick==1) then
            if ( already_thick == 0 ) then
              gamma_edd_switch = gamma_edd
              Mdot_switch = ABS(s% mstar_dot/Msun*secyer)
              L_switch = L
              Mhom_switch = Mhom
            end if
            thick_met = .true.
          end if
        elseif ( s% x_character_ctrl(1) == 'Xsurf'  ) then
          if (X<=0.4) then
            thick_met = .true.
          end if
        elseif ( s% x_character_ctrl(1) == 'gamma_eta') then
          if ( gamma_edd >= gamma_trans .or. eta>=eta_trans .or. already_thick==1 ) then
            thick_met = .true.
          end if
        elseif ( s% x_character_ctrl(1) == 'any') then
          if ( gamma_edd >= gamma_trans .or. eta>=eta_trans .or. X<=0.4 .or. already_thick==1 ) then
            thick_met = .true.
          end if
        end if

        if ( (s% x_character_ctrl(10) == "Xsurf" .and. X <=s%  x_ctrl(2)) .or. &
        (s% x_character_ctrl(10) == "Teff" .and. Tsurf >= s% x_ctrl(2)) ) then
          HPoor_WR_condition = .true.
        else
          HPoor_WR_condition = .false.
        end if

        ! -----------------------------------------------------------------------------------------------------------------------

        if(s% x_logical_ctrl(5)) then ! Belczynski+2010 LBV2 winds (eq. 8) with factor 1
           if (s% center_h1 < 1.0d-3) then  ! postMS
             if (L/Lsun > 6.0d5 .and. 1.0d-5 * R/Rsun * (L/Lsun)**0.5d0 > 1.0d0) then ! Humphreys-Davidson limit
               wind_scheme = 90.0
               w  = 1.0d-4
               s% max_years_for_timestep = 1d2
               write(*,*) "Here are LBV winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
               return ! This is a way of saying to not check anything below if we reached LBV winds in this phase
             endif
           endif
        endif


      if ( Tsurf < s% x_ctrl(4) ) then                                          ! Dust-driven winds

          s% max_years_for_timestep = 1d2                                       ! To have more resolution during this phase

          call eval_dust_winds(w)

      elseif (.not. thick_met) then                                             ! Thin winds part
          call eval_thin_winds(w)

      elseif(thick_met) then
          already_thick = 1

          if ((Tsurf/1000 < 30 .and. s% x_logical_ctrl(2) .and. &               ! Cool WR winds, if requested
          s% x_character_ctrl(5) /= "V11") .or. (.not. HPoor_WR_condition .and. gamma_edd < 0.4)) then                                 ! If it is V11, we go through Vink/Sabhahit procedure
            call eval_thin_winds(w)
          else
            call eval_thick_winds(w)
          end if
      end if

        !  ---------------------- Thick winds interpolation -------------------------

        ! Here I check which thick winds condition is used and if we are near to that switch
        !  In case of multiple conditions, the philosophy is also: first come, first served

        gamma_condition = ("gamma" == s%x_character_ctrl(1) .or. &
        "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1)) .and. &
                  gamma_edd > gamma_trans - s%x_ctrl(3) .and. &
                  gamma_edd < gamma_trans + s%x_ctrl(3) .and. &
                  .not. ((eta > eta_trans + s%x_ctrl(3) .and. ("eta" == s%x_character_ctrl(1) .or. &
                  "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))) .or. &
                  (X < 0.4 - s%x_ctrl(3) .and. ("Xsurf" == s%x_character_ctrl(1) .or. &
                  "any" == s%x_character_ctrl(1))))

        eta_condition = ("eta" == s%x_character_ctrl(1) .or. &
        "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1)) .and. &
                eta > eta_trans - s%x_ctrl(3) .and. &
                eta < eta_trans + s%x_ctrl(3) .and. &
                .not. ((gamma_edd > gamma_trans + s%x_ctrl(3) .and. ("gamma" == s%x_character_ctrl(1) .or. &
                "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))) .or. &
                (X < 0.4 - s%x_ctrl(1) .and. ("Xsurf" == s%x_character_ctrl(1) .or. &
                "any" == s%x_character_ctrl(1))))

        xsurf_condition = ("Xsurf" == s%x_character_ctrl(1) .or. &
        "any" == s%x_character_ctrl(1)) .and. &
                X > 0.4 - s%x_ctrl(3) .and. &
                X < 0.4 + s%x_ctrl(3) .and. &
                .not. ((eta > eta_trans + s%x_ctrl(3) .and. ("eta" == s%x_character_ctrl(1) .or. &
                "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))) .or. &
                (gamma_edd > gamma_trans + s%x_ctrl(3) .and. ("gamma" == s%x_character_ctrl(1) .or. &
                "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))))

        ! These conditions above are a way to say "If I am within this range AND
        !   I am not well deep in thick winds due to another transition condition,
        !   then interpolate"

        if (s% x_logical_ctrl(4) .and. Tsurf > s% hot_wind_full_on_T .and. &
        (gamma_condition .or. eta_condition .or. xsurf_condition) .and. s% x_character_ctrl(5) /= "V11") then

            if ( Tsurf/1000 < 30 .and. ((gmlogg>3.0d0 .and. s% x_character_ctrl(2) == s% x_character_ctrl(4)) .or. &
            (gmlogg<=3.0d0 .and. s% x_character_ctrl(3) == s% x_character_ctrl(4))) ) then
              write(*,*) "Cool WR winds == Thin winds ; no need for transition"

            else
              write(*,*) "Near optically-thick winds threshold interpolation"
              wind_scheme_interp = wind_scheme

              if (Tsurf/1000 < 30 .and. s% x_logical_ctrl(2)) then
                call eval_thin_winds(w1)
              else
                call eval_thick_winds(w1)
              end if

            call eval_thin_winds(w2)

            if(s% x_ctrl(3) == 0) then
               w = 0.5d0*(w1 + w2)
              else
                divisor = 2*s%x_ctrl(3)
              if (gamma_condition) then
                beta = min( (0.5+s%x_ctrl(3) - gamma_edd) / divisor, 1d0)
              else if (eta_condition) then
                beta = min( (eta_trans+s%x_ctrl(3) - eta) / divisor, 1d0)
              else
                beta = min( (X+s%x_ctrl(3) - 0.4) / divisor, 1d0)
              end if
              alfa_mid = 1d0 - beta
              w = alfa_mid*w1 + beta*w2
            end if
            write(*,*) "Interpolated mass loss: log(Mdot [Msun/yr]) =", log10(ABS(w))
            wind_scheme=wind_scheme_interp                                      !This is to not mess up with the tags
          end if
        end if


        !  ---------------------- low Teff interpolation -------------------------


      if (s% cool_wind_full_on_T<=Tsurf .and. Tsurf<=s% hot_wind_full_on_T .and. &
      s% x_logical_ctrl(6)) then

        write(*,*) "Near dust-driven winds threshold interpolation"
        wind_scheme_interp = wind_scheme

        call eval_dust_winds(w1)

        if (.not. thick_met) then                                               ! No need to combine the two interpolations because
                                                                                !  if I am close to dust-winds, I get at most cool WR
            call eval_thin_winds(w2)                                            !  winds. So even if I am near to the thick winds
                                                                                !  transition, mass loss rates are not about to super-change

       else
          if (Tsurf/1000 < 30 .and. s% x_logical_ctrl(2) .and. &
          s% x_character_ctrl(5) /= "V11") then
            call eval_thin_winds(w2)
          else
            call eval_thick_winds(w2)
          end if
          if(s% hot_wind_full_on_T == s% cool_wind_full_on_T)then
             w = 0.5d0*(w1 + w2)
          else
            divisor = s% hot_wind_full_on_T - s% cool_wind_full_on_T
            beta = min( (s% hot_wind_full_on_T - Tsurf) / divisor, 1d0)
            alfa_mid = 1d0 - beta
            w = alfa_mid*w2 + beta*w1
          end if

          write(*,*) "Interpolated mass loss: log(Mdot [Msun/yr]) =", log10(ABS(w))
          wind_scheme=wind_scheme_interp                                        !This is to not mess up with the tags
        end if
      end if

      gamma_edd_old = gamma_edd

      !  -----------------------------------------------------------------------

    contains

      subroutine eval_Reimers_wind(w)
         real(dp), intent(out) :: w
         include 'formats'

         if ( s% center_h1 > 0.001 ) then                                       ! GENEC calibrations
           s% Reimers_scaling_factor = 0.85d0
         elseif (L/s% prev_Ledd > 1) then
           s% Reimers_scaling_factor = 3.0d0
         else
           s% Reimers_scaling_factor = 1.0d0
         end if

         wind_scheme = 10.0

         w = 4d-13*(L*R/M)/(Lsun*Rsun/Msun)
         w = w*s% Reimers_scaling_factor

      end subroutine eval_Reimers_wind

      subroutine eval_Blocker_wind(w)
        ! Bloecker, T. 1995, A&A, 297, 727
         real(dp), intent(out) :: w
         include 'formats'
         s% Blocker_scaling_factor = 0.1
         w = 4.83d-9*(M/Msun)**(-2.1)*(L/Lsun)**2.7*4d-13*(L/Lsun)*(R/Rsun)/(M/Msun)
         w = w*s% Blocker_scaling_factor

         wind_scheme = 11.0

      end subroutine eval_Blocker_wind

        subroutine eval_Krticka24_wind(w)
          real(dp), intent(inout) :: w
          real(dp) :: hehratio

          wind_scheme = 24.0

          xlmdot=-13.82d0+0.358*logZ_div_Zsun+(1.52d0-0.11*logZ_div_Zsun)*(log10(L/Lsun)-6.d0) &
         + 13.82d0*log10((1.0+0.73*logZ_div_Zsun)*exp(-((Tsurf/1000.0d0-14.16d0)/3.58d0)**2.d0) &
         + 3.84d0*exp(-((Tsurf/1000.0d0-37.9d0)/56.5d0)**2.d0))

          w=10**xlmdot

        end subroutine   eval_Krticka24_wind

        subroutine eval_GormazMatamala23_wind(w)
          real(dp), intent(inout) :: w
          real(dp) :: hehratio

          wind_scheme = 20.0

          xlmdot=-40.314+15.438*lteff+45.838/gmlogg-8.284*lteff/gmlogg+1.0564*gmrstar
          xlmdot=xlmdot-lteff*gmrstar/2.36-1.1967*gmrstar/gmlogg+11.6*logZ_div_Zsun
          xlmdot=xlmdot-4.223*lteff*logZ_div_Zsun-16.377*logZ_div_Zsun/gmlogg+(gmrstar*logZ_div_Zsun)/81.735
          !hehratio=0.25*(Y/X)    !Alex said it doesn't change much
          hehratio=0.085
          xlmdot=xlmdot+0.0475-0.559*hehratio

          w=10**xlmdot

        end subroutine   eval_GormazMatamala23_wind


        subroutine eval_Vink01_wind(w)
          real(dp), intent(inout) :: w
          real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

          wind_scheme = 22.0

          ! alfa = 1 for hot side, = 0 for cool side
          if (Tsurf > 27500d0) then
            alfa = 1
          else if (Tsurf < 22500d0) then
            alfa = 0
          else
            ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
            Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*logZ_div_Zsun))
            dT = 2000d0!100d0
            if (Tsurf > Teff_jump + dT) then
              alfa = 1
              ! wind_scheme = 1.1
            else if (Tsurf < Teff_jump - dT) then
              alfa = 0
              ! wind_scheme = 1.3
            else
              alfa = 0.5d0*(Tsurf - (Teff_jump - dT)) / dT
            end if
            ! wind_scheme = 1.2
          end if

          if (alfa > 0) then ! eval hot side wind (eqn 24)
            vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
            vinf_div_vesc = vinf_div_vesc*Z_div_Z_solar**0.13d0 ! corrected for Z; previously **0.13d0
            logMdot = &
            - 6.697d0 &
            + 2.194d0*log10(L/Lsun/1d5) &
            - 1.313d0*log10(M/Msun/30d0) &
            - 1.226d0*log10(vinf_div_vesc/2d0) &
            + 0.933d0*log10(Tsurf/4d4) &
            - 10.92d0*(log10(Tsurf/4d4))**2 &
            + 0.85d0*logZ_div_Zsun
            w1 = 10**(logMdot)
          else
            w1 = 0
          end if


          if (alfa < 1) then ! eval cool side wind (eqn 25)
            vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
            vinf_div_vesc = vinf_div_vesc*Z_div_Z_solar**0.13d0 ! corrected for Z; previously **0.13d0
            logMdot = &
            - 6.688d0 &
            + 2.210d0*log10(L/Lsun/1d5) &
            - 1.339d0*log10(M/Msun/30d0) &
            - 1.601d0*log10(vinf_div_vesc/2d0) &
            + 1.07d0*log10(Tsurf/2d4) &
            + 0.85d0*logZ_div_Zsun
            w2 = 10**(logMdot)
          else
            w2 = 0
          end if

          w = alfa*w1 + (1 - alfa)*w2

        end subroutine eval_Vink01_wind

        subroutine eval_VinkSander21_wind(w)
          real(dp), intent(inout) :: w
          real(dp) :: logMdot, vinf_div_vesc

          wind_scheme = 25.0

          if (Tsurf > 2.5d4) then ! eval hot side wind (eqn 24)
            vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
            vinf_div_vesc = vinf_div_vesc*Z_div_Z_solar**0.13d0 ! corrected for Z
            logMdot = &
            - 6.697d0 &
            + 2.194d0*log10(L/Lsun/1d5) &
            - 1.313d0*log10(M/Msun/30d0) &
            - 1.226d0*log10(vinf_div_vesc/2d0) &
            + 0.933d0*log10(Tsurf/4d4) &
            - 10.92d0*(log10(Tsurf/4d4))**2 &
            + 0.85d0*logZ_div_Zsun

          elseif (Tsurf > 2.0d4) then
            vinf_div_vesc = 0.7d0
            vinf_div_vesc = vinf_div_vesc*Z_div_Z_solar**0.13d0

            logMdot = &
            - 5.99d0 &
            + 2.21d0*log10(L/Lsun/1d5) &
            - 1.339d0*log10(M/Msun/30d0) &
            - 1.601d0*log10(vinf_div_vesc/2d0) &
            + 1.07d0*log10(Tsurf/2d4) &
            + 0.85d0*logZ_div_Zsun


          else  ! eval cool side wind
            vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
            vinf_div_vesc = vinf_div_vesc*Z_div_Z_solar**0.13d0 ! corrected for Z
            logMdot = &
            - 6.688d0 &
            + 2.210d0*log10(L/Lsun/1d5) &
            - 1.339d0*log10(M/Msun/30d0) &
            - 1.601d0*log10(vinf_div_vesc/2d0) &
            + 1.07d0*log10(Tsurf/2d4) &
            + 0.85d0*logZ_div_Zsun
          end if

          w = 10**(logMdot)

        end subroutine eval_VinkSander21_wind

        subroutine eval_Pauli25_wind(w)
           real(dp), intent(inout) :: w
           real(dp) :: xlmdot, Meff

           wind_scheme = 26.0

           ! eq 5 from Pauli et al, 2025, A&A, 697 (2025) A114
           xlmdot = -3.92 + 4.27*log_gamma_edd + 0.86*logZ_div_Zsun
           w = exp10(xlmdot)

        end subroutine eval_Pauli25_wind

        subroutine eval_Bjorklund23_wind(w)
           real(dp), intent(inout) :: w
           real(dp) :: xlmdot, Meff

           wind_scheme = 21.0

           ! "This recipe is valid within the ranges 4.5 <= log L/LSun <= 6.0,
           !    15 <= M/LSun <= 80, 15 000K <= Teff <= 50 000 K, and 0.2 <= Z/ZSun <= 1.0"

           ! electron opacity is constant 0.34 in their models (Eq. 6)
           Meff = M*(1d0 - 0.34d0*L/(pi4*clight*s% cgrav(1)*M))  ! effective mass

           ! eq 7 from Björklund et al, 2023, A&A, Volume 676, id.A109, 14 pp.
           xlmdot = - 5.52d0 &
                  + 2.39d0 * log10(L/(1d6*Lsun)) &
                  - 1.48d0 * log10(M/(4.5d1*Msun)) &
                  + 2.12d0 * log10(Tsurf/4.5d4) &
                  + (0.75d0 - 1.87d0 * log10(Tsurf/4.5d4)) * logZ_div_Zsun
           w = exp10(xlmdot)

        end subroutine eval_Bjorklund23_wind

        subroutine eval_Vink17_wind(w)
           real(dp), intent(inout) :: w
           real(dp) :: xlmdot

           wind_scheme = 27.0

           xlmdot = - 13.3d0 &
                  + 1.36d0 * log10(L/Lsun) &
                  + 0.61 * logZ_div_Zsun
           w = exp10(xlmdot)

           write(*,*) "Here are V17 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))

        end subroutine eval_Vink17_wind

        subroutine eval_thin_winds(w)
          real(dp), intent(inout) :: w

          if ( X < 0.4 .and. .not. thick_met) then
            if ( s% x_character_ctrl(9) == "V17" ) then
              call eval_Vink17_wind(w)

            elseif ( s% x_character_ctrl(9) == "P25" ) then
              call eval_Pauli25_wind(w)
            end if


          elseif (gmlogg>3.0d0 .and. .not. thick_met) then
            if ( s%x_character_ctrl(2) =='V01' ) then
              call eval_Vink01_wind(w)
              write(*,*) "Here are V01 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(2) =='VS21') then
              call eval_VinkSander21_wind(w)
              write(*,*) "Here are VS21 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(2) =='GM23') then
              call eval_GormazMatamala23_wind(w)
              write(*,*) "Here are GM23 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(2) =='Bj23') then
              call eval_Bjorklund23_wind(w)
              write(*,*) "Here are Bj23 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(2) =='K24') then
              call eval_Krticka24_wind(w)
              write(*,*) "Here are K24 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(2) =='P25') then
              call eval_Pauli25_wind(w)
              write(*,*) "Here are P25 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            end if

          else if (gmlogg<=3.0d0 .and. .not. thick_met) then
            if ( s%x_character_ctrl(3) =='V01' ) then
              call eval_Vink01_wind(w)
              write(*,*) "Here are V01 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(3) =='VS21') then
              call eval_VinkSander21_wind(w)
              write(*,*) "Here are VS21 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(3) =='Bj23') then
              call eval_Bjorklund23_wind(w)
              write(*,*) "Here are Bj23 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(3) =='K24') then
              call eval_Krticka24_wind(w)
              write(*,*) "Here are K24 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(3) =='P25') then
              call eval_Pauli25_wind(w)
              write(*,*) "Here are P25 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            end if

          else
            if ( s%x_character_ctrl(4) =='V01' ) then
              call eval_Vink01_wind(w)
              write(*,*) "Here are V01 (thick) winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(4) =='VS21') then
              call eval_VinkSander21_wind(w)
              write(*,*) "Here are VS21 (thick) winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(4) =='Bj23') then
              call eval_Bjorklund23_wind(w)
              write(*,*) "Here are Bj23 (thick) winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(4) =='K24') then
              call eval_Krticka24_wind(w)
              write(*,*) "Here are K24 (thick) winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            elseif (s%x_character_ctrl(4) =='P25') then
              call eval_Pauli25_wind(w)
              write(*,*) "Here are P25 (thick) winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            end if

          end if

        end subroutine eval_thin_winds


        subroutine eval_Vink11_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: xlmdot

            wind_scheme = 40.0

            if (gamma_edd >= gamma_edd_switch .and. gamma_edd >= gamma_edd_old ) then
              wind_scheme = 4.5

              xlmdot = log10(ABS(Mdot_switch)) + 4.77d0*log10(L/L_switch) - 3.99d0*log10(M/(Mhom_switch*Msun))
              ! write(*,*) "Mdot ", ABS(Mdot_switch), "Lswitch ", log10(L/L_switch), "Mhom_switch ", log10(M/(Mhom_switch*Msun))
              w = 10**(xlmdot)
              write(*,*) "Here are V11 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
            else
              call eval_thin_winds(w)

            end if

         end subroutine eval_Vink11_wind

         subroutine eval_Bestenlehner20_wind(w)
          real(dp), intent(inout) :: w
          real(dp) :: xlmdot

          wind_scheme = 41.0


         !*** Bestenlehner (2020) prescription for hot stars
         !*** with fitting parameters from Brands et al. (2022)
          xlmdot = -5.19d0 + 2.69d0 * log10(gamma_edd) - 3.19d0 * log10(1-gamma_edd)
          xlmdot = xlmdot + (logZ_div_Zsun+0.3d0)*(0.4+15.75d0/M)

          write(*,*) "Here are B20 winds: log(Mdot [Msun/yr]) =", xlmdot

          w = 10**(xlmdot)


        end subroutine eval_Bestenlehner20_wind

        subroutine eval_GrafenerHamann08_wind(w)
         real(dp), intent(inout) :: w
         real(dp) :: xlmdot
         ! Grafener, G. & Hamann, W.-R. 2008, A&A 482, 945

         wind_scheme = 42.0

         xlmdot = 10.046 + 1.727*log10(gamma_edd-0.326) - 3.5*log10(Tsurf) + 0.42*log10(L/Lsun) - 0.45*X

         write(*,*) "Here are GH08 winds: log(Mdot [Msun/yr]) =", xlmdot

         w = 10**(xlmdot)

       end subroutine eval_GrafenerHamann08_wind

        subroutine eval_Yoon06_wind(w)
         real(dp), intent(inout) :: w
         real(dp) :: xlmdot

         wind_scheme = 43.0

         if ( log10(L/Lsun) <= 4.5 ) then
           xlmdot = -36.8 + 6.8*log10(L/Lsun)-2.85*X + 0.85*logZ_div_Zsun
         else
           xlmdot = -12.95 + 1.5*log10(L/Lsun) - 2.85*X + 0.85*logZ_div_Zsun
         end if

         write(*,*) "Here are Y06 winds: log(Mdot [Msun/yr]) =", xlmdot

         w = 10**(xlmdot)

       end subroutine eval_Yoon06_wind

        subroutine eval_NugisLamers_wind(w)
         real(dp), intent(inout) :: w
         real(dp) :: xlmdot

         wind_scheme = 44.0

          ! If I want an universal NL00 for both WN and WC/WO I take this one below

          !xlmdot = log10(1d-11 * (L/Lsun)**1.29d0 * Y**1.7d0 * sqrt(Z))  ! Default MESA setup for everything

          ! Addition of the calibrations from Eldridge & Vink (2006), taken from the study of late-type WN and WC
          !   mass loss predictions of Vink & de Koter (2005). Also this is the GENEC model
           if (.not. HPoor_WR_condition .or. Z<=0.03d0) then
             ! WN
             xlmdot=-13.60d0+1.63d0*logL_div_Lsun+2.22d0*log10(Y)+0.85d0*logZ_div_Zsun
           else
             ! WC + WO
             if (s% kap_rq% Zbase > Zsolar ) then          ! Zinit>Zsolar
               xlmdot=-8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.40d0*logZ_div_Zsun
             else if (s% kap_rq% Zbase  >  0.002d0) then
               xlmdot=-8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.66d0*logZ_div_Zsun
             else
               if (s% kap_rq% Zbase  <  0.00000001d0) then
                    xlmdot = -8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.66d0*log10(0.002d0/Zsolar)+ &
                            0.35d0*log10(Z/0.002d0)

              else

               xlmdot = -8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.66d0*log10(0.002d0/Zsolar)+ &
                           0.35d0*log10(s% kap_rq% Zbase/0.002d0)

             end if
           end if
         end if

         write(*,*) "Here are NL00 winds: log(Mdot [Msun/yr]) =", xlmdot


       w = 10**(xlmdot)


      end subroutine eval_NugisLamers_wind

      subroutine eval_SanderVink20_wind(w)
        real(dp), intent(inout) :: w
        real(dp) :: mdg_a,mdg_cbd,mdg_geddb,mdg_logMdotOff,logMdot_breakdown,logMdot_pureWR
        real(dp) :: xlmdot

        wind_scheme = 45.0

        !*** Sander & Vink formula for WR winds as a function of Gamma_e
        mdg_a = 2.932
        mdg_geddb = -0.324*logZ_div_Zsun+0.244
        mdg_cbd = -0.44*logZ_div_Zsun+9.15
        mdg_logMdotOff = 0.23*logZ_div_Zsun-2.61
        logMdot_pureWR = mdg_a*(log10(-log10(1.0-gamma_edd))) + mdg_logMdotOff
        logMdot_breakdown = log10(2.0) * (mdg_geddb/gamma_edd)**(mdg_cbd)
        xlmdot = logMdot_pureWR - logMdot_breakdown
        !    xlmdot = xlmdot - 6.0d0*log10(teff/141000.0d0) ! Sander et al. 2023

        write(*,*) "Here are SV20 winds: log(Mdot [Msun/yr]) =", xlmdot

        w = 10**(xlmdot)
        wind_scheme = 5.0


      end subroutine eval_SanderVink20_wind


         subroutine eval_thick_winds(w)
           real(dp), intent(inout) :: w
           real(dp) :: logMdot

           include 'formats'
           ! gamma_edd = exp10(-4.813d0)*(1+xsurf)*(L/Lsun)*(Msun/M)

           if (.not. HPoor_WR_condition .or. .not. s% x_logical_ctrl(3)) then         ! H-rich WR winds

             if ( s%x_character_ctrl(5)=='B20' ) then

               call eval_Bestenlehner20_wind(w)

             elseif ( s%x_character_ctrl(5)=='V11' ) then

               call eval_Vink11_wind(w)

             elseif ( s%x_character_ctrl(5)=='GH08' ) then

               call eval_GrafenerHamann08_wind(w)

             elseif ( s%x_character_ctrl(5)=='Y06' ) then

               call eval_Yoon06_wind(w)

             elseif ( s%x_character_ctrl(5)=='NL00' ) then

               call eval_NugisLamers_wind(w)

             end if


           else

             if (  s%x_character_ctrl(6)=='SV20' ) then

               call eval_SanderVink20_wind(w)

             elseif( s%x_character_ctrl(6)=='Y06') then
               call eval_Yoon06_wind(w)

            elseif( s%x_character_ctrl(6)=='NL00') then
              call eval_NugisLamers_wind(w)

             end if


           end if

           ! w = w * scaling_factor ! No need for the scaling factor now

         end subroutine eval_thick_winds

        subroutine eval_Antoniadis24_wind(w)
           ! Antoniadis, K., et al.: A&A, 686, A88 (2024)
           real(dp), intent(out) :: w
           real(dp) :: log10w
           include 'formats'
           wind_scheme = 34.0

           if ( log(L/Lsun) < 4.4 ) then
             log10w = 0.26*log10(L/Lsun) - 14.19*log10(Tsurf/4000)-9.17
           else
             log10w = 2.49536114*log10(L/Lsun) -33.70573385*log10(Tsurf/4000)-19.10219818
           end if

           w = exp10(log10w)

           write(*,*) "Here are A24 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))

        end subroutine eval_Antoniadis24_wind

        subroutine eval_van_Loon05_wind(w)
           real(dp), intent(out) :: w
           real(dp) :: log10w
           include 'formats'
           wind_scheme = 31.0

           log10w = 1.05 * log10(L/Lsun/1d4) - 6.3 * log10(Tsurf/3500) - 5.56
           w = exp10(log10w)
           write(*,*) "Here are vL05 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
        end subroutine eval_van_Loon05_wind

        subroutine eval_Beasor23_wind(w)
           real(dp), intent(out) :: w
           real(dp) :: log10w
           include 'formats'
           wind_scheme = 32.0

           log10w = -0.15 * s% initial_mass + 3.6 * log10(L/Lsun) - 21.5
           w = exp10(log10w)
           write(*,*) "Here are B23 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
        end subroutine eval_Beasor23_wind

        subroutine eval_Yang23_wind(w)
           real(dp), intent(out) :: w
           real(dp) :: log10w
           include 'formats'
           wind_scheme = 33.0

           log10w = 0.45 * log10(L/Lsun)**3 - 5.26 * (log10(L/Lsun))**2 + 20.93 * log10(L/Lsun) - 34.56
           w = exp10(log10w)
           write(*,*) "Here are Y23 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
        end subroutine eval_Yang23_wind

        subroutine eval_Decin24_wind(w)
           real(dp), intent(out) :: w
           real(dp) :: log10w
           include 'formats'
           wind_scheme = 35.0

           log10w = 1.71 -1.63 * s% initial_mass/10 + 3.47 * log10(L/Lsun/1d5) - 5
           w = exp10(log10w)
           write(*,*) "Here are D24 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
        end subroutine eval_Decin24_wind

        subroutine eval_de_Jager88_wind(w)
           ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
           real(dp), intent(out) :: w
           real(dp) :: log10w
           include 'formats'
           wind_scheme = 30.0

           log10w = 1.769d0*log10(L/Lsun) - 1.676d0*log10(Tsurf) - 8.158d0
           w = exp10(log10w)
           write(*,*) "Here are dJ88 winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
        end subroutine eval_de_Jager88_wind

        subroutine eval_dust_winds(w)
          real(dp), intent(inout) :: w
          integer :: which_dust                                     ! Which dust wind to use (in case of multple initiations)
          real(dp) :: logMdot

          include 'formats'

          if (s% x_ctrl(4) /= s% x_ctrl(5) .and. Tsurf < s% x_ctrl(5) .and. &   ! Second switch to dust-driven winds. For now these recipes are the only two implemented
             log10(L/Lsun) <= 5.8) then
            which_dust = 7
          else
            which_dust = 8
          end if


          if ( s%x_character_ctrl(which_dust)=='dJ88' ) then

            call eval_de_Jager88_wind(w)

          elseif ( s%x_character_ctrl(which_dust)=='vL05' ) then

            call eval_van_Loon05_wind(w)

          elseif ( s%x_character_ctrl(which_dust)=='B23' ) then

            call eval_Beasor23_wind(w)

          elseif ( s%x_character_ctrl(which_dust)=='Y23' ) then

            call eval_Yang23_wind(w)

          elseif ( s%x_character_ctrl(which_dust)=='A24' ) then

            call eval_Antoniadis24_wind(w)

          elseif ( s%x_character_ctrl(which_dust)=='D24' ) then

            call eval_Decin24_wind(w)


          end if


          ! w = w * scaling_factor ! No need for the scaling factor now

        end subroutine eval_dust_winds

      end subroutine my_other_wind

      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: vrot_max,tot_L_notSi
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if ( s% initial_mass < 2 .and. s% center_h1 < 1d-3) then
           s% delta_lgT_cntr_limit = 0.05d0                                     ! Stars that enter the RGB need a bit more flexibility for the central temperature
         end if

         extras_check_model = keep_going

        if (s% center_h1 < 1d-3 .and. s% center_he4 < 1d-3 .and. &
        s% center_c12 < s% x_ctrl(1) .and. s% x_logical_ctrl(1)) then
         write (*,*) "End of core-C burning. We are happy with that -- Terminate"
         extras_check_model = terminate
        end if
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 5
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         real(dp) :: Msum,Rsum
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = "wind_scheme"
         names(2) = "eta"
         names(3) = "gamma_edd"
         names(4) = "eta_trans"
         names(5) = "Radius_99perc"

         Msum=0
         Rsum=0

         do i = s% nz, 1, -1
            Msum = Msum + s% m(i)/Msun

            if (Msum>=s% star_mass*0.99) then
              exit
            end if

            Rsum = Rsum + s% r(i)/Rsun
         end do

         vals(1) = wind_scheme
         vals(2) = eta
         vals(3) = gamma_edd
         vals(4) = eta_trans
         vals(5) = Rsum
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = "rot_momentum"
         do k=1,nz
           vals(k,1) = s% j_rot(k)*s% m(k)
         end do
      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         if (s% center_c12<=1d-4 .and. s% center_h1<=1d-4) then


           ! ---------- From Schenider 2024 -----------------------------

           s% tol_residual_norm1 = 1d-5
           s% tol_max_residual1 = 1d-2
           s% iter_for_resid_tol2 = 5
           !tol_residual_norm2 = 1d-3
           !tol_max_residual2 = 1d-2
           s% iter_for_resid_tol3 = 12
           s% min_timestep_limit = 1d-12 ! (seconds)

          s% mesh_delta_coeff = 2.5
          s% mesh_delta_coeff_for_highT = 0.4
          s% max_dq = 1d-3 ! should set a minimum of 2000 zones for each model
          s% delta_lgRho_cntr_limit = 1.2d-2
          s% delta_lgT_cntr_limit = 2.0d-3
          s% varcontrol_target = 1d-4 !7d-4
          s% dX_nuc_drop_limit = 5.0d-2 !1.0d-3
          s% dX_nuc_drop_limit_at_high_T = 2d-3 !5d-3
          s% logT_max_for_standard_mesh_delta_coeff = 9.0
          s% logT_min_for_highT_mesh_delta_coeff = 9.1
          !s% dH_div_H_limit = 0.1d0

          ! high center T limit to avoid negative mass fractions
          s% sig_min_factor_for_high_Tcenter = 0.015
           ! inactive when >= 1d0
             ! if Tcenter >= Tcenter_min_for_sig_min_factor_full_on,
             ! then okay to reduce sig by as much as this factor
             ! as needed to prevent causing negative abundances
          s% Tcenter_min_for_sig_min_factor_full_on = 3.2d9
             ! if Tcenter >= this, factor = sig_min_factor_for_neg_abundances
             ! this should be > Tcenter_max_for_sig_min_factor_full_off.
          s% Tcenter_max_for_sig_min_factor_full_off = 2.8d9
             ! if Tcenter <= this, factor = 1, so has no effect
             ! this should be < Tcenter_min_for_sig_min_factor_full_on.
          ! for T > full_off and < full_on, factor changes linearly with Tcenter

          ! ---------- End Schenider 2024  -----------------------------

           s% delta_HR_hard_limit = -1
           s% adjust_J_q_hard_limit = -1
           s% min_xa_hard_limit = -1

           s% hydro_mtx_max_allowed_logT = 99d0
           s% hydro_mtx_max_allowed_logRho = 99d0

           s% hydro_mtx_min_allowed_logT = -1d5!-1d0
           s% hydro_mtx_min_allowed_logRho = -1d2

           s% maxT_for_gold_tolerances = 5d8

           s% adjust_J_q_limit = 0.5
           s% sum_xa_hard_limit = 5d-4
           s% logT_max_for_sum_xa_hard_limit = -1!10d0

           ! s% time_delta_coeff = 0.5
          !  s% mesh_delta_coeff = 0.5
           !at this point, allow the code to evaluate radial velocities
             !this allows the model to core-collapse.
             call star_set_v_flag(s% id, .true., ierr)
             if (ierr /= 0) then
                write(*,*) "Failure in changing v flag"
                extras_finish_step = terminate
                return ! failure in profile
             end if
             !activate velocities ONLY for the internal, hot meshes
             !I don't care about pulsations, only for the collapse.
             s% velocity_logT_lower_bound = 8.8
             s% max_dt_yrs_for_velocity_logT_lower_bound = 0.1
         end if

      end function extras_finish_step



      end module run_star_extras
