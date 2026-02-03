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
!  20.0: NdJ90 (Nieuwenhuijzen & de Jager 1990) --> Also for cool supergiants!
!  21.0: V01   (Vink+ 2001)
!  21.5: VS21  (Vink & Sander 2021)
!  22.0: V17   (Vink 2017)
!  23.0: Bj23  (Bjorklund+ 2023)
!  24.0: GM23  (Gormaz-Matamala+ 2023)
!  25.0: K24   (Krticka+ 2024)
!  25.5: K25   (Krticka+ 2025) --> Only for Z < 0.2 Zsun!
!  26.0: P25   (Pauli+ 2025) --> Also for WR and He stars!
!  27.0: Sa25  (Sabhahit+ 2025) --> Only for very massive stars!
!
! Dust/Cool Winds:
!  30.0: dJ88  (de Jager+ 1988)
!  31.0: vL05  (van Loon+ 2005)
!  32.0: Be23  (Beasor+ 2023)
!  33.0: Ya23  (Yang+ 2023)
!  34.0: A24   (Antoniadis+ 2024)
!  35.0: D24   (Decin+ 2024)
!
! Thick/WR Winds:
!  40.0: L89   (Langer 1989  + Vink, de Koter 2005)
!  40.5: Ha98  (Hamann & Koesterke 1998 + Vink, de Koter 2005)
!  41.0: NL00  (Nugis & Lamers 2000)
!  41.5: Y06   (Yoon+ 2006)
!  42.0: GH08  (Grafener & Hamann 2008)
!  42.5: V11   (Vink+ 2011)
!  43.0: TSK16 (Tramper, Sana & de Koter 2016)
!  43.5: Y17   (Yoon 2017)
!  44.0: Sh19  (Shenar+ 2019 plus 2020 erratum)
!  44.5: S19   (Sander+ 2019 + Vink 2015 "True WR origin"; calibrations from Iorio+ 2021)
!  45.0: B20   (Bestenlehner 2020)
!  45.5: SV20  (Sander & Vink 2020)
!
! Special Cases:
!  90.0: Vb98   (Vanbeveren+ 1998 + Vink, de Koter 2005 for WR/OB; three different formulae for OB, WR, and cool supergiants)
!  90.5: H00    (Hurley+ 2000; LBV)
!  91.0: Bk10   (Belczynski+ 2010; LBV)
!  91.5: Ch24   (Cheng+ 2024; LBV)
!  92.0: P26    (Pauli+ 2026; LBV)
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

      logical :: already_thick = .true.
      real(dp) :: Mdot_switch,L_switch,M_switch,gamma_edd_switch
      real(dp) :: eta,eta_trans,gamma_edd,gamma_edd_old                         ! gamma_edd_old is to check the previous timestep
      real(dp) :: wind_scheme,wind_scheme_interp                                ! To know which winds model I am using at each timestep

      ! this is used to soften too large changes in wind mass loss rates
      real(dp) :: old_wind = 0

      real(dp) :: max_years_dt_old
      ! these routines are called by the standard run_star check_model
      contains

      function round_3(val) result(res)
      real(dp), intent(in) :: val
      real(dp) :: res
      res = nint(val * 1000.0d0) / 1000.0d0
      end function round_3

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

         if (s% kap_rq% Zbase/0.0142d0 > 0.5) then
           s% use_superad_reduction = .true.
         else
           s% use_superad_reduction = .true.
         end if

         s% overshoot_f(1) = f_ov_fcn_of_mass(s% initial_mass)
         s% overshoot_f0(1) = s% overshoot_f(1)/100

         s% overshoot_f(2) = f_ov_fcn_of_mass(s% initial_mass)/10
         s% overshoot_f0(2) = s% overshoot_f0(2)/100

      end subroutine extras_startup

     function f_ov_fcn_of_mass(m) result(f_ov)
       real(dp), intent(in) :: m
       real(dp) :: f_ov, frac
       real(dp), parameter :: f1 = 1.6d-2, f2=4.15d-2
       if(m < 4.0d0) then
         frac = 0.0d0
       else if(m > 8.0d0) then
         frac = 1.0d0
       else
         frac = 0.5d0 * (1.0d0 - cos(0.25d0 * (m - 4.0d0) *  3.1415927))
       endif

       f_ov = f1 + (f2-f1)*frac

       if (m>=20) then
         f_ov = 5.0d-2
       end if

     end function f_ov_fcn_of_mass


     subroutine my_other_wind(id, L, M, R, Tsurf, X, Y, Z, w, ierr)
       integer, intent(in) :: id
       real(dp), intent(in) :: L, M, R, Tsurf, X, Y, Z ! surface values (cgs)
       real(dp), intent(out) :: w !wind in units of Msun/year (value is >= 0)
       real(dp) :: gmrstar,gmlogg,lteff,logMdot,logZ_div_Zsun,hehratio,vterm,Zsolar
       real(dp) :: Z_div_Z_solar,Teff_jump,alfa,log_gamma_edd,gamma_trans,logL_div_Lsun
       real(dp) :: vesc_eff,vesc,vinf_fac,const_k,Vinf_e
       real(dp) :: w1, w2
       real(dp) :: divisor,beta,alfa_mid
       character(len=10) :: vinf_type
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

       Zsolar = s% x_ctrl(9)

       Z_div_Z_solar = s% kap_rq% Zbase/Zsolar
       logZ_div_Zsun=log10(Z_div_Z_solar)
       logL_div_Lsun=log10(L/Lsun)

       w = 0

       ! -----------------------------------------------

       if ( s% x_character_ctrl(1) /= 'eta' .and. s% x_character_ctrl(5) == 'V11' ) then
         call mesa_error(__FILE__,__LINE__,'V11 winds can only be computed with the eta condition')
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

       gamma_trans = s% x_ctrl(2)

       ! if ( logZ_div_Zsun >= log10(0.2) ) then
       !   gamma_trans = 0.5
       ! else
       !   gamma_trans = 0.5 - 0.301*logZ_div_Zsun-0.045*logZ_div_Zsun**2
       ! end if

       log_gamma_edd = -4.813d0+log10(1+X)+log10(L/Lsun)-log10(M/Msun)
       gamma_edd=10**log_gamma_edd
       gamma_edd = round_3(gamma_edd)
       ! gamma_edd = L/s% prev_Ledd

       vesc = sqrt(2d0*standard_cgrav*M/R)/1d5

       vinf_type = "Hawcroft24"

       select case (trim(vinf_type))
       case ("Lamers95")                                                        ! Usual fits for bistability jump
            if ( Tsurf > 25000 ) then
              Vinf_e = 2.6
            else
              Vinf_e = 1.3
            end if
            vterm = Vinf_e * sqrt(2d0*standard_cgrav*(M)*(1-gamma_edd)/R)/1d5*Z_div_Z_solar**0.20d0
          case ("Crowther06")                                                   ! Fits from Atlas I from Crowther+ (2006). If you want Lamers, put it at 2.6
            Vinf_e = 1.59d-4*Tsurf - 0.89
            vterm = Vinf_e * sqrt(2d0*standard_cgrav*(M)*(1-gamma_edd)/R)/1d5*Z_div_Z_solar**0.20d0
          case ("Hawcroft24")                                                   ! Fits from Atlas I from Hawcroft+ (2024).
            vterm = (9.2d-2*Tsurf-1040)*Z_div_Z_solar**0.22d0
          case ("Kritcka25")                                                    ! Fits from K25.
            vterm = (107+25*logZ_div_Zsun)*(Tsurf/1d3)-1190-430*logZ_div_Zsun
          case default
            ierr = 1 ! Handle unknown string
       end select

       const_k = (clight*Msun*1d5)/(Lsun*3600d0*24d0*365d0)           ! constant to evaluate eta

       eta_trans = 0.75/(1+(vesc**2)/(vterm**2))
       ! eta = (ABS(s% mstar_dot /Msun)*secyer * vterm)/(L/(clight))
       eta = const_k*(ABS(s% mstar_dot/Msun*secyer)*vterm)/(L/Lsun)

       eta = round_3(eta)
       eta_trans = round_3(eta_trans)


       write(*,*)
       write(*,*) "------------------------------------------------------------"
       if ( already_thick ) then
         write(*, '(A, f5.3, A, f5.3)') " Gamma_e: ", gamma_edd, " vs Gamma_switch: ", gamma_edd_switch,  " vs Gamma_trans: ", gamma_trans
       else
         write(*, '(A, f5.3, A, f5.3)') " Gamma_e: ", gamma_edd, " vs Gamma_trans: ", gamma_trans
       end if
       write(*, '(A, f7.3)') " previous log(Mdot): ", log10(abs(s% mstar_dot/Msun*secyer))
       write(*, '(A, f0.3, A, f0.3, A)') " vterm: ", vterm, " km/s vs vesc: ", vesc, " km/s"
       write(*, '(A, f5.3, A, f5.3)') " eta factor: ", eta, " vs eta trans: ", eta_trans
       write(*,*) "------------------------------------------------------------"
       write(*,*)

       ! --------------------------------------------- Check for thick winds ---------------------------------------------------
       thick_met = .false.

       ! Logic: If x_logical_ctrl(8) is .true., the "sticky" behavior is enabled.
       ! If it's .false., already_thick is ignored for the initial thick_met check,
       ! allowing the star to revert to thin winds if physical conditions are no longer met.

       if ( s% x_character_ctrl(1) == 'gamma' ) then
         if ( gamma_edd >= gamma_trans ) thick_met = .true.
       elseif (s% x_character_ctrl(1) == 'eta') then
         if (eta>=eta_trans) thick_met = .true.
       elseif ( s% x_character_ctrl(1) == 'Xsurf' ) then
         if (X<=0.4) thick_met = .true.
       elseif ( s% x_character_ctrl(1) == 'gamma_eta') then
         if ( gamma_edd >= gamma_trans .or. eta>=eta_trans ) thick_met = .true.
       elseif ( s% x_character_ctrl(1) == 'any') then
         if ( gamma_edd >= gamma_trans .or. eta>=eta_trans .or. X<=0.4 ) thick_met = .true.
       end if

       if ( (s% x_character_ctrl(10) == "Xsurf" .and. X <=s%  x_ctrl(3)) .or. &
              (s% x_character_ctrl(10) == "Teff" .and. Tsurf >= s% x_ctrl(3)) ) then
            HPoor_WR_condition = .true.
       else
            HPoor_WR_condition = .false.
       end if

       if (.not. thick_met .and. .not. s% x_logical_ctrl(8)) then
         ! Reset already_thick if we are no longer in the thick regime and sticky is off
         already_thick = .false.
       elsez`
         if (.not. already_thick) then
           ! Capture state for V11 or other transitions at the first crossing
           gamma_edd_switch = gamma_edd
           gamma_edd_switch = round_3(gamma_edd_switch)
           call eval_thin_winds(w)
           L_switch = L
           Mdot_switch = w
           M_switch = M
           already_thick = .true.
         end if
         thick_met = .true.
       end if

     ! --------------------------------------------- Rate Selection ---------------------------------------------------

     if ( Tsurf < s% x_ctrl(7) ) then                                          ! cool supergiant winds

         ! s% max_years_for_timestep = 3d2                                       ! To have more resolution during this phase

         call eval_cool_wind(w)

     elseif (.not. thick_met) then                                             ! Thin winds part
         call eval_thin_winds(w)

     elseif(thick_met) then

         if ((Tsurf/1000 < s% x_ctrl(5) .and. s% x_logical_ctrl(2))) then    ! Cool WR winds, if requested
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
                 gamma_edd > gamma_trans - s%x_ctrl(4) .and. &
                 gamma_edd < gamma_trans + s%x_ctrl(4) .and. &
                 .not. ((eta > eta_trans + s%x_ctrl(4) .and. ("eta" == s%x_character_ctrl(1) .or. &
                 "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))) .or. &
                 (X < 0.4 - s%x_ctrl(4) .and. ("Xsurf" == s%x_character_ctrl(1) .or. &
                 "any" == s%x_character_ctrl(1))))

       eta_condition = ("eta" == s%x_character_ctrl(1) .or. &
       "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1)) .and. &
               eta > eta_trans - s%x_ctrl(4) .and. &
               eta < eta_trans + s%x_ctrl(4) .and. &
               .not. ((gamma_edd > gamma_trans + s%x_ctrl(4) .and. ("gamma" == s%x_character_ctrl(1) .or. &
               "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))) .or. &
               (X < 0.4 - s%x_ctrl(1) .and. ("Xsurf" == s%x_character_ctrl(1) .or. &
               "any" == s%x_character_ctrl(1))))

       xsurf_condition = ("Xsurf" == s%x_character_ctrl(1) .or. &
       "any" == s%x_character_ctrl(1)) .and. &
               X > 0.4 - s%x_ctrl(4) .and. &
               X < 0.4 + s%x_ctrl(4) .and. &
               .not. ((eta > eta_trans + s%x_ctrl(4) .and. ("eta" == s%x_character_ctrl(1) .or. &
               "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))) .or. &
               (gamma_edd > gamma_trans + s%x_ctrl(4) .and. ("gamma" == s%x_character_ctrl(1) .or. &
               "gamma_eta" == s%x_character_ctrl(1) .or. "any" == s%x_character_ctrl(1))))

       ! These conditions above are a way to say "If I am within this range AND
       !   I am not well deep in thick winds due to another transition condition,
       !   then interpolate"

       if (s% x_logical_ctrl(4) .and. Tsurf > s% hot_wind_full_on_T .and. &
       (gamma_condition .or. eta_condition .or. xsurf_condition) .and. s% x_character_ctrl(5) /= "V11") then

           if ( Tsurf/1000 < s% x_ctrl(5) .and. ((gmlogg>3.0d0 .and. s% x_character_ctrl(2) == s% x_character_ctrl(4)) .or. &
           (gmlogg<=3.0d0 .and. s% x_character_ctrl(3) == s% x_character_ctrl(4))) ) then
             write(*,*) "Cool WR winds == Thin winds ; no need for transition"

           else
             write(*,*) "Near optically-thick winds threshold interpolation"
             wind_scheme_interp = wind_scheme

             if (Tsurf/1000 < s% x_ctrl(5) .and. s% x_logical_ctrl(2)) then
               call eval_thin_winds(w1)
             else
               call eval_thick_winds(w1)
             end if

           call eval_thin_winds(w2)

           if(s% x_ctrl(4) == 0) then
              w = 0.5d0*(w1 + w2)
             else
               divisor = 2*s%x_ctrl(4)
             if (gamma_condition) then
               beta = min( (0.5+s%x_ctrl(4) - gamma_edd) / divisor, 1d0)
             else if (eta_condition) then
               beta = min( (eta_trans+s%x_ctrl(4) - eta) / divisor, 1d0)
             else
               beta = min( (X+s%x_ctrl(4) - 0.4) / divisor, 1d0)
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
     s% x_logical_ctrl(5)) then

       write(*,*) "Near cool supergiant winds threshold interpolation"
       wind_scheme_interp = wind_scheme

       call eval_cool_wind(w1)

       if (.not. thick_met) then                                               ! No need to combine the two interpolations because
                                                                               !  if I am close to cool winds, I get at most cool WR
           call eval_thin_winds(w2)                                            !  winds. So even if I am near to the thick winds
                                                                               !  transition, mass loss rates are not about to super-change

      else
         if (Tsurf/1000 < s% x_ctrl(5) .and. s% x_logical_ctrl(2) .and. &
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

     !  -----------------------------LBV Winds---------------------------------

    if ( s% x_logical_ctrl(6) ) then

      call eval_LBV_winds(w)

    end if

     !  -----------------------------------------------------------------------

     ! ------ STOP WINDS WHEN WE REACH THE END OF CORE-C BURNING -------------
     if (s% center_c12<=1d-2 .and. s% center_he4<=1d-3 .and. s% center_h1<=1d-3) then

       w = 0
       write(*,*) "We are past core-C burning. Winds are turned off"


     end if

     old_wind = w

   contains

     ! ====================================================================
     !  SMOOTHING ROUTINE
     ! ====================================================================
     subroutine smooth_wind_log(w, scheme_name)
       real(dp), intent(inout) :: w
       character(len=*), intent(in) :: scheme_name
       logical :: is_smoothed

       is_smoothed = .false.

       ! Pauli26-like logic: prevent the wind from dropping too fast.
       ! This avoids numerical noise and "crazy jumps" downwards.
       ! We allow rapid increases (no check for w > old_wind).
       if (old_wind /= 0d0) then
         if (abs(w) < abs(old_wind) * 0.5d0 .and. s% x_logical_ctrl(7)) then
           w = old_wind * 0.5d0
           is_smoothed = .true.
         end if
       end if

       if (is_smoothed) then
         write(*,*) "Here are ", trim(scheme_name), " winds (smoothed): log(Mdot [Msun/yr]) =", log10(ABS(w))
       else
         write(*,*) "Here are ", trim(scheme_name), " winds: log(Mdot [Msun/yr]) =", log10(ABS(w))
       end if

     end subroutine smooth_wind_log
! ====================================================================

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

   call smooth_wind_log(w, "Reimers RGB cool")

end subroutine eval_Reimers_wind

subroutine eval_Blocker_wind(w)
  ! Bloecker, T. 1995, A&A, 297, 727
   real(dp), intent(out) :: w
   include 'formats'
   s% Blocker_scaling_factor = 0.1
   w = 4.83d-9*(M/Msun)**(-2.1)*(L/Lsun)**2.7*4d-13*(L/Lsun)*(R/Rsun)/(M/Msun)
   w = w*s% Blocker_scaling_factor

   wind_scheme = 11.0
   call smooth_wind_log(w, "Blocker AGB cool")

end subroutine eval_Blocker_wind

  subroutine eval_Krticka24_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: hehratio

    wind_scheme = 25.0

    logMdot=-13.82d0+0.358*logZ_div_Zsun+(1.52d0-0.11*logZ_div_Zsun)*(log10(L/Lsun)-6.d0) &
   + 13.82d0*log10((1.0+0.73*logZ_div_Zsun)*exp(-((Tsurf/1000.0d0-14.16d0)/3.58d0)**2.d0) &
   + 3.84d0*exp(-((Tsurf/1000.0d0-37.9d0)/56.5d0)**2.d0))

    w=10**logMdot
    call smooth_wind_log(w, "K24")

  end subroutine   eval_Krticka24_wind

  subroutine eval_Krticka25_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: t, l, z, dmdt
    real(dp) :: p1, p2, p3, p4, p5, p6, p7
    real(dp) :: c1, c2, c3, c4, c5, c6, c7

    ! Parameters from fit
    real(dp), parameter :: a     = -7.d0
    real(dp), parameter :: fb    = 1.619d0
    real(dp), parameter :: fl1   = -31.333d0
    real(dp), parameter :: fl2   = 119.62d0
    real(dp), parameter :: fl3   = -201.713d0
    real(dp), parameter :: fl4   = 208.466d0
    real(dp), parameter :: fl5   = -140.704d0
    real(dp), parameter :: fl6   = 58.7226d0
    real(dp), parameter :: fl7   = -11.7687d0
    real(dp), parameter :: flz1  = -46.3518d0
    real(dp), parameter :: flz2  = 170.353d0
    real(dp), parameter :: flz3  = -284.575d0
    real(dp), parameter :: flz4  = 292.063d0
    real(dp), parameter :: flz5  = -195.337d0
    real(dp), parameter :: flz6  = 80.3569d0
    real(dp), parameter :: flz7  = -15.3933d0
    real(dp), parameter :: flzz1 = -19.0704d0
    real(dp), parameter :: flzz2 = 70.4817d0
    real(dp), parameter :: flzz3 = -115.793d0
    real(dp), parameter :: flzz4 = 117.131d0
    real(dp), parameter :: flzz5 = -77.1403d0
    real(dp), parameter :: flzz6 = 31.2397d0
    real(dp), parameter :: flzz7 = -5.729d0
    real(dp), parameter :: flzzz1= 0.372222d0
    real(dp), parameter :: flzzz2= 0.119759d0
    real(dp), parameter :: flzzz3= 0.0d0
    real(dp), parameter :: flzzz4= 0.0d0
    real(dp), parameter :: flzzz5= 0.0d0
    real(dp), parameter :: flzzz6= 0.0d0
    real(dp), parameter :: flzzz7= 0.0d0

    wind_scheme = 25.5

    ! --------------------------------------------------------------------
    ! Krticka et al. (2025) - Re-implementation with specific fits
    !   from the Zenodo files https://zenodo.org/records/15965163
    ! --------------------------------------------------------------------

    ! Normalized Temperature: t maps [10kK, 45kK] to [0, 1]
    t = (Tsurf/1000.d0 - 10.d0)/35.d0

    ! Clamping to avoid crash if Tsurf is slightly out of bounds
    if(t > 1.d0) then
        write(*,*) "Warning: K25 Tsurf > 45kK, clamping to fit boundary."
        t = 1.d0
    elseif(t < 0.d0) then
        write(*,*) "Warning: K25 Tsurf < 10kK, clamping to fit boundary."
        t = 0.d0
    endif

    ! Normalized Metallicity
    z = logZ_div_Zsun
    ! Clamping metallicity to valid range [-2, 0]
    if(z < -2.d0) then
        z = -2.d0
    elseif(z > 0.d0) then
        z = 0.d0
    endif

    ! Normalized Luminosity
    l = logL_div_Lsun - 6.d0

    ! --- Polynomials p(t) ---
    ! p1(t)=t
    p1 = t
    ! p2(t)= -0.1d1 / 0.2d1 + 0.3d1 / 0.2d1 * t * t
    p2 = -0.5d0 + 1.5d0*t*t
    ! p3(t)= 0.5d1 / 0.2d1 * t**3 - 0.3d1 / 0.2d1 * t
    p3 = 2.5d0*t**3 - 1.5d0*t
    ! p4(t)= 0.3d1 / 0.8d1 + 0.35d2 / 0.8d1 * t**4 - 0.15d2 / 0.4d1 * t**2
    p4 = 0.375d0 + 4.375d0*t**4 - 3.75d0*t**2
    ! p5(t)= 0.63d2 / 0.8d1 * t**5 - 0.35d2 / 0.4d1 * t**3 + 0.15d2 / 0.8d1 * t
    p5 = 7.875d0*t**5 - 8.75d0*t**3 + 1.875d0*t
    ! p6(t)= -0.5d1 / 0.16d2 + 0.231d3 / 0.16d2 * t**6 - 0.315d3 / 0.16d2 * t**4 + 0.105d3 / 0.16d2 * t**2
    p6 = -0.3125d0 + 14.4375d0*t**6 - 19.6875d0*t**4 + 6.5625d0*t**2
    ! p7(t)= 0.429d3 / 0.16d2 * t**7 - 0.693d3 / 0.16d2 * t**5 + 0.315d3 / 0.16d2 * t**3 - 0.35d2 / 0.16d2 * t
    p7 = 26.8125d0*t**7 - 43.3125d0*t**5 + 19.6875d0*t**3 - 2.1875d0*t

    ! --- Coefficients c(z) ---
    ! c_n(z) = l_n + lz_n*z + lzz_n*z^2 + lzzz_n*z^3
    c1 = fl1 + flz1*z + flzz1*z*z + flzzz1*z**3
    c2 = fl2 + flz2*z + flzz2*z*z + flzzz2*z**3
    c3 = fl3 + flz3*z + flzz3*z*z + flzzz3*z**3
    c4 = fl4 + flz4*z + flzz4*z*z + flzzz4*z**3
    c5 = fl5 + flz5*z + flzz5*z*z + flzzz5*z**3
    c6 = fl6 + flz6*z + flzz6*z*z + flzzz6*z**3
    c7 = fl7 + flz7*z + flzz7*z*z + flzzz7*z**3

    ! --- Final Mass Loss Calculation ---
    ! logMdot = a + b*l + Sum(c_n * p_n)
    logMdot = a + fb*l + c1*p1 + c2*p2 + c3*p3 + c4*p4 + c5*p5 + c6*p6 + c7*p7

    w = 10.d0**logMdot
    call smooth_wind_log(w, "K25")

  end subroutine eval_Krticka25_wind

  subroutine eval_GormazMatamala23_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: hehratio

    wind_scheme = 24.0

    logMdot=-40.314+15.438*lteff+45.838/gmlogg-8.284*lteff/gmlogg+1.0564*gmrstar
    logMdot=logMdot-lteff*gmrstar/2.36-1.1967*gmrstar/gmlogg+11.6*logZ_div_Zsun
    logMdot=logMdot-4.223*lteff*logZ_div_Zsun-16.377*logZ_div_Zsun/gmlogg+(gmrstar*logZ_div_Zsun)/81.735
    !hehratio=0.25*(Y/X)    !Alex said it doesn't change much
    hehratio=0.085
    logMdot=logMdot+0.0475-0.559*hehratio

    w=10**logMdot
    call smooth_wind_log(w, "GM23")

  end subroutine   eval_GormazMatamala23_wind

  subroutine eval_Nieuwenhuijzen_deJager90_wind(w)
    real(dp), intent(inout) :: w

    wind_scheme = 20.0

    w = 9.6d-15*(R/Rsun)**0.81*(L/Lsun)**1.24*(M/Msun)**0.16*Z_div_Z_solar**0.5
    call smooth_wind_log(w, "NdJ90")

  end subroutine   eval_Nieuwenhuijzen_deJager90_wind

  subroutine eval_Vanbeveren98_wind(w,which)
    real(dp), intent(inout) :: w
    character(len=*), intent(in) :: which

    wind_scheme = 90.0

    !*** + metallicity dependence from Vink & de Koeter (2005)

    if ( which == "OB" ) then
      logMdot = 1.67*log10(L/Lsun) -1.55*log10(Tsurf) + 0.85*logZ_div_Zsun -8.29
    elseif ( which == "cool" ) then
      logMdot = 0.8*log10(L/Lsun) -8.7
    elseif ( which == "WR" ) then
      logMdot = log10(L/Lsun) + 0.85*logZ_div_Zsun -10
    end if

    w=10**logMdot
    call smooth_wind_log(w, "Vb98 (" // trim(which) // ")")

  end subroutine   eval_Vanbeveren98_wind


  subroutine eval_Vink01_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

    wind_scheme = 21.0

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

    call smooth_wind_log(w, "V01")

  end subroutine eval_Vink01_wind

  subroutine eval_VinkSander21_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: logMdot, vinf_div_vesc

    wind_scheme = 21.5

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
    call smooth_wind_log(w, "VS21")

  end subroutine eval_VinkSander21_wind

  subroutine eval_Pauli25_wind(w)
     real(dp), intent(inout) :: w
     real(dp) :: logMdot, Meff

     wind_scheme = 26.0

     ! eq 5 from Pauli et al, 2025, A&A, 697 (2025) A114
     logMdot = -3.92 + 4.27*log_gamma_edd + 0.86*logZ_div_Zsun
     w = 10**logMdot
     call smooth_wind_log(w, "P25")

  end subroutine eval_Pauli25_wind

  subroutine eval_Bjorklund23_wind(w)
     real(dp), intent(inout) :: w
     real(dp) :: logMdot, Meff

     wind_scheme = 23.0

     ! "This recipe is valid within the ranges 4.5 <= log L/LSun <= 6.0,
     !    15 <= M/LSun <= 80, 15 000K <= Teff <= 50 000 K, and 0.2 <= Z/ZSun <= 1.0"

     ! electron opacity is constant 0.34 in their models (Eq. 6)
     Meff = M*(1d0 - 0.34d0*L/(pi4*clight*s% cgrav(1)*M))  ! effective mass

     ! eq 7 from Björklund et al, 2023, A&A, Volume 676, id.A109, 14 pp.
     logMdot = - 5.52d0 &
            + 2.39d0 * log10(L/(1d6*Lsun)) &
            - 1.48d0 * log10(M/(4.5d1*Msun)) &
            + 2.12d0 * log10(Tsurf/4.5d4) &
            + (0.75d0 - 1.87d0 * log10(Tsurf/4.5d4)) * logZ_div_Zsun
     w = 10**logMdot
     call smooth_wind_log(w, "Bj23")

  end subroutine eval_Bjorklund23_wind

  subroutine eval_Vink17_wind(w)
     real(dp), intent(inout) :: w
     real(dp) :: logMdot

     wind_scheme = 22.0

     logMdot = - 13.3d0 &
            + 1.36d0 * log10(L/Lsun) &
            + 0.61 * logZ_div_Zsun
     w = 10**logMdot

     call smooth_wind_log(w, "V17")

  end subroutine eval_Vink17_wind

  subroutine eval_Sabhahit25_wind(w)
     real(dp), intent(inout) :: w
     real(dp) :: logMdot
     real(dp) :: a,b,c,T1,T2,G1,G2,Tref,Mref,f_low,f_high,sigmaT

     wind_scheme = 27.0

     a = -5.527
     Tref = 38
     Mref = 60 - 0.521*(Tsurf/1000-Tref)
     f_low = -1.864*log10(M/Msun/Mref)
     f_high = -9.865*log10(M/Msun/Mref)
     b = -2.062 + 5.671*log10(M/Msun/Mref)
     c = 3.974
     G1 = 0.178*exp(9.611*log10(M/Msun/Mref))
     G2 = 0.291*exp(8.658*log10(M/Msun/Mref))
     T2 = 16.941-7.274*log10(M/Msun/Mref)
     T1 = 25
     sigmaT = 2

     logMdot = a + log10(10**f_low + 10**f_high) &
     + b * log10((Tsurf/1000)/Tref) + c * log10((Tsurf/1000)/Tref)**2 &
     -G1 * exp(-((Tsurf/1000-T1)/sigmaT)**2) &
     -G2 * exp(-((Tsurf/1000-T2)/sigmaT)**2)

     w = 10**logMdot

     call smooth_wind_log(w, "Sa25")

  end subroutine eval_Sabhahit25_wind

  subroutine eval_thin_winds(w)
    real(dp), intent(inout) :: w

    if ( X < 0.4 .and. .not. thick_met) then
      if ( s% x_character_ctrl(9) == "V17" ) then
        call eval_Vink17_wind(w)

      elseif ( s% x_character_ctrl(9) == "P25" ) then
        call eval_Pauli25_wind(w)
      end if


    elseif (gmlogg>3.0d0 .and. .not. thick_met) then
      if ( s%x_character_ctrl(2) =='NdJ90' ) then
        call eval_Nieuwenhuijzen_deJager90_wind(w)
      elseif ( s%x_character_ctrl(2) =='Vb98' ) then
        call eval_Vanbeveren98_wind(w,"OB")
      elseif ( s%x_character_ctrl(2) =='V01' ) then
        call eval_Vink01_wind(w)
      elseif (s%x_character_ctrl(2) =='VS21') then
        call eval_VinkSander21_wind(w)
      elseif (s%x_character_ctrl(2) =='GM23') then
        call eval_GormazMatamala23_wind(w)
      elseif (s%x_character_ctrl(2) =='Bj23') then
        call eval_Bjorklund23_wind(w)
      elseif (s%x_character_ctrl(2) =='K24') then
        call eval_Krticka24_wind(w)
      elseif (s%x_character_ctrl(2) =='K25') then
        call eval_Krticka25_wind(w)
      elseif (s%x_character_ctrl(2) =='P25') then
        call eval_Pauli25_wind(w)
      elseif (s%x_character_ctrl(2) =='Sa25') then
        call eval_Sabhahit25_wind(w)
      end if


    else if (gmlogg<=3.0d0 .and. .not. thick_met) then
      if ( s%x_character_ctrl(3) =='NdJ90' ) then
        call eval_Nieuwenhuijzen_deJager90_wind(w)
      elseif ( s%x_character_ctrl(3) =='Vb98' ) then
        call eval_Vanbeveren98_wind(w,"OB")
      elseif ( s%x_character_ctrl(3) =='V01' ) then
        call eval_Vink01_wind(w)
      elseif (s%x_character_ctrl(3) =='VS21') then
        call eval_VinkSander21_wind(w)
      elseif (s%x_character_ctrl(3) =='Bj23') then
        call eval_Bjorklund23_wind(w)
      elseif (s%x_character_ctrl(3) =='K24') then
        call eval_Krticka24_wind(w)
      elseif (s%x_character_ctrl(3) =='K25') then
        call eval_Krticka25_wind(w)
      elseif (s%x_character_ctrl(3) =='P25') then
        call eval_Pauli25_wind(w)
      elseif (s%x_character_ctrl(3) =='Sa25') then
        call eval_Sabhahit25_wind(w)
      end if

    else
      if ( s%x_character_ctrl(4) =='NdJ90' ) then
        call eval_Nieuwenhuijzen_deJager90_wind(w)
      elseif ( s%x_character_ctrl(4) =='Vb98' ) then
        call eval_Vanbeveren98_wind(w,"OB")
      elseif ( s%x_character_ctrl(4) =='V01' ) then
        call eval_Vink01_wind(w)
      elseif (s%x_character_ctrl(4) =='VS21') then
        call eval_VinkSander21_wind(w)
      elseif (s%x_character_ctrl(4) =='Bj23') then
        call eval_Bjorklund23_wind(w)
      elseif (s%x_character_ctrl(4) =='K24') then
        call eval_Krticka24_wind(w)
      elseif (s%x_character_ctrl(4) =='K25') then
        call eval_Krticka25_wind(w)
      elseif (s%x_character_ctrl(4) =='P25') then
        call eval_Pauli25_wind(w)
      elseif (s%x_character_ctrl(4) =='Sa25') then
        call eval_Sabhahit25_wind(w)
      end if

    end if

  end subroutine eval_thin_winds


  subroutine eval_Vink11_wind(w)
      real(dp), intent(inout) :: w
      real(dp) :: logMdot


      if (gamma_edd >= gamma_edd_switch) then
        wind_scheme = 42.5

        logMdot = log10(ABS(Mdot_switch)) + 4.77d0*log10(L/L_switch) - 3.99d0*log10(M/M_switch)
        ! write(*,*) "Mdot ", ABS(Mdot_switch), "Lswitch ", log10(L/L_switch), "M_switch ", log10(M/(M_switch*Msun))
        w = 10**(logMdot)
        call smooth_wind_log(w, "V11")
      else
        call eval_thin_winds(w)

      end if

   end subroutine eval_Vink11_wind

   subroutine eval_Bestenlehner20_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: logMdot

    wind_scheme = 45.0


   !*** Bestenlehner (2020) prescription for hot stars
   !*** with fitting parameters from Brands et al. (2022)
   if(.not. HPoor_WR_condition .and. gamma_edd < 0.4) then                      ! It underestimates Mdot at low gamma, hence this condition
     call eval_thin_winds(w)                                                    !   if HPoor_WR_condition, it just goes to H-poor WR winds

   else
     logMdot = -5.19d0 + 2.69d0 * log10(gamma_edd) - 3.19d0 * log10(1-gamma_edd)
     logMdot = logMdot + (logZ_div_Zsun+0.3d0)*(0.4+15.75d0/M)

     w = 10**(logMdot)

     call smooth_wind_log(w, "B20")

   end if


  end subroutine eval_Bestenlehner20_wind

  subroutine eval_GrafenerHamann08_wind(w)
   real(dp), intent(inout) :: w
   real(dp) :: logMdot
   ! Grafener, G. & Hamann, W.-R. 2008, A&A 482, 945

   wind_scheme = 42.0

   logMdot = 10.046 + 1.727*log10(gamma_edd-0.326) - 3.5*log10(Tsurf) + 0.42*log10(L/Lsun) - 0.45*X

   w = 10**(logMdot)
   call smooth_wind_log(w, "GH08")

 end subroutine eval_GrafenerHamann08_wind

  subroutine eval_Yoon06_wind(w)
   real(dp), intent(inout) :: w
   real(dp) :: logMdot

   wind_scheme = 41.5

   if ( log10(L/Lsun) <= 4.5 ) then
     logMdot = -36.8 + 6.8*log10(L/Lsun)-2.85*X + 0.85*logZ_div_Zsun
   else
     logMdot = -12.95 + 1.5*log10(L/Lsun) - 2.85*X + 0.85*logZ_div_Zsun
   end if

   w = 10**(logMdot)
   call smooth_wind_log(w, "Y06")

 end subroutine eval_Yoon06_wind

  subroutine eval_NugisLamers_wind(w)
   real(dp), intent(inout) :: w
   real(dp) :: logMdot

   wind_scheme = 41.0

    ! If I want an universal NL00 for both WN and WC/WO I take this one below

    if ( .not. s% x_logical_ctrl(3) ) then

      logMdot = log10(1d-11 * (L/Lsun)**1.29d0 * Y**1.7d0 * sqrt(Z))  ! Default MESA setup for everything

     ! Addition of the calibrations from Eldridge & Vink (2006), taken from the study of late-type WN and WC
    !   mass loss predictions of Vink & de Koter (2005). Also this is the GENEC model
   else if (.not. HPoor_WR_condition .or. Z<=0.03d0) then
       ! WN
       logMdot=-13.60d0+1.63d0*logL_div_Lsun+2.22d0*log10(Y)+0.85d0*logZ_div_Zsun
   else
       ! WC + WO
       if (s% kap_rq% Zbase > Zsolar ) then          ! Zinit>Zsolar
         logMdot=-8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.40d0*logZ_div_Zsun
       else if (s% kap_rq% Zbase  >  0.002d0) then
         logMdot=-8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.66d0*logZ_div_Zsun
       else
         if (s% kap_rq% Zbase  <  0.00000001d0) then
              logMdot = -8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.66d0*log10(0.002d0/Zsolar)+ &
                      0.35d0*log10(Z/0.002d0)

        else

         logMdot = -8.30d0+0.84d0*logL_div_Lsun+2.04d0*log10(Y)+1.04d0*log10(Z)+0.66d0*log10(0.002d0/Zsolar)+ &
                     0.35d0*log10(s% kap_rq% Zbase/0.002d0)

       end if
     end if
   end if


 w = 10**(logMdot)
 call smooth_wind_log(w, "NL00")


end subroutine eval_NugisLamers_wind

subroutine eval_Langer89_wind(w)
 real(dp), intent(inout) :: w
 real(dp) :: logMdot

 wind_scheme = 40.0

!*** + metallicity dependence from Vink & de Koeter (2005)

 logMdot = 2.5*log10(M/Msun) + 0.85*logZ_div_Zsun - 7.1

 w = 10**(logMdot)
 call smooth_wind_log(w, "L89")

end subroutine eval_Langer89_wind

subroutine eval_Hamann98_wind(w)
 real(dp), intent(inout) :: w
 real(dp) :: logMdot

 wind_scheme = 40.5

!*** + metallicity dependence from Vink & de Koeter (2005)

 logMdot = 1.5*log10(L/Lsun) + 0.85*logZ_div_Zsun - 13.0

 w = 10**(logMdot)
 call smooth_wind_log(w, "Ha98")


end subroutine eval_Hamann98_wind

subroutine eval_TSK16_wind(w)
 real(dp), intent(inout) :: w
 real(dp) :: logMdot

 wind_scheme = 43.0


 logMdot = 0.85*log10(L/Lsun) + 0.44*log10(Y) + 0.25*logZ_div_Zsun - 9.2

 w = 10**(logMdot)
 call smooth_wind_log(w, "TSK16")

end subroutine eval_TSK16_wind

subroutine eval_Yoon17_wind(w)
 real(dp), intent(inout) :: w
 real(dp) :: f_WR

 f_WR = 1                                                                       ! Should be 1 for clumping factor D = 4 and 1.58 for D = 10
 wind_scheme = 43.5


 w = f_WR*((L/Lsun)**1.18)*(Z_div_Z_solar**0.6)*(10**(-11.32))
 call smooth_wind_log(w, "Y17")

end subroutine eval_Yoon17_wind


subroutine eval_Shenar19_wind(w)
 real(dp), intent(inout) :: w
 real(dp) :: logMdot
 real(dp) :: C1,C2,C3,C4,C5

 wind_scheme = 44.0

 C4 = 1.42

!*** Shenar+ (2019)
!*** erratum from 2020

if ( X >= 0.4 ) then
  C1 = -6.50
  C2 = 0.79
  C3 = -0.37
  C5 = 0.68
elseif ( X > 0.2 ) then
  C1 = -4.02
  C2 = 0.60
  C3 = -0.74
  C5 = 0.43
elseif ( X > 0.05 ) then
  C1 = -3.84
  C2 = 0.76
  C3 = -0.78
  C5 = 0.81
else
  C1 = -8.13
  C2 = 1.01
  C3 = -0.06
  C5 = 0.95
end if

 logMdot = C1 + C2 * log10(L/Lsun) + C3 * log10(Tsurf) + C4* log10(Y) + C5*log10(Z)

 w = 10**(logMdot)
 call smooth_wind_log(w, "Sh19")


end subroutine eval_Shenar19_wind

subroutine eval_Sander19_wind(w)
 real(dp), intent(inout) :: w
 real(dp) :: logMdot,fWN,fWCO,Zscale

 wind_scheme = 44.5

 fWN = -1 + 1.9*tanh(0.58*log10(s% xa(s% net_iso(ife56),0) )+1)
 fWCO = -0.3 + 1.2*tanh(0.5*log10(s% xa(s% net_iso(ife56),0) )+0.5)

 ! metallicity dependence adopted from fits in Costa et al. (2021)
 if (.not. HPoor_WR_condition .or. .not. s% x_logical_ctrl(3)) then
   Zscale = fWN
 else
   Zscale = fWCO

 end if

 logMdot = -8.31 + 0.68*log10(L/Lsun)
 ! logMdot = logMdot+Zscale                                                  ! I do not use Z-calibrations,
                                                                           !  I don't get how to implement them

 w = 10**logMdot
 call smooth_wind_log(w, "S19")


end subroutine eval_Sander19_wind


subroutine eval_SanderVink20_wind(w)
  real(dp), intent(inout) :: w
  real(dp) :: mdg_a,mdg_cbd,mdg_geddb,mdg_logMdotOff,logMdot_breakdown,logMdot_pureWR
  real(dp) :: logMdot

  wind_scheme = 45.5

  !*** Sander & Vink formula for WR winds as a function of Gamma_e
  mdg_a = 2.932
  mdg_geddb = -0.324*logZ_div_Zsun+0.244
  mdg_cbd = -0.44*logZ_div_Zsun+9.15
  mdg_logMdotOff = 0.23*logZ_div_Zsun-2.61
  logMdot_pureWR = mdg_a*(log10(-log10(1.0-gamma_edd))) + mdg_logMdotOff
  logMdot_breakdown = log10(2.0) * (mdg_geddb/gamma_edd)**(mdg_cbd)
  logMdot = logMdot_pureWR - logMdot_breakdown
  !    logMdot = logMdot - 6.0d0*log10(teff/141000.0d0) ! Sander et al. 2023

  w = 10**(logMdot)
  wind_scheme = 5.0

  call smooth_wind_log(w, "SV20")


end subroutine eval_SanderVink20_wind

subroutine eval_thick_winds(w)
     real(dp), intent(inout) :: w
     real(dp) :: logMdot

     include 'formats'
     ! gamma_edd = exp10(-4.813d0)*(1+xsurf)*(L/Lsun)*(Msun/M)

     if (.not. HPoor_WR_condition .or. .not. s% x_logical_ctrl(3)) then         ! H-rich WR winds

       if ( s%x_character_ctrl(5)=='B20' ) then

         call eval_Bestenlehner20_wind(w)

       elseif ( s%x_character_ctrl(5)=='Sh19' ) then

         call eval_Shenar19_wind(w)

       elseif ( s%x_character_ctrl(5)=='S19' ) then

         call eval_Sander19_wind(w)

       elseif ( s%x_character_ctrl(5)=='V11' ) then

         call eval_Vink11_wind(w)

       elseif ( s%x_character_ctrl(5)=='GH08' ) then

         call eval_GrafenerHamann08_wind(w)

       elseif ( s%x_character_ctrl(5)=='Y06' ) then

         call eval_Yoon06_wind(w)

       elseif ( s%x_character_ctrl(5)=='NL00' ) then

         call eval_NugisLamers_wind(w)

       elseif ( s%x_character_ctrl(5)=='Vb98' ) then

         call eval_Vanbeveren98_wind(w,"WR")

       elseif ( s%x_character_ctrl(5)=='L89' ) then

         call eval_Langer89_wind(w)

       elseif ( s%x_character_ctrl(5)=='Ha98' ) then

         call eval_Hamann98_wind(w)
       end if


     else

       if (  s%x_character_ctrl(6)=='SV20' ) then

         call eval_SanderVink20_wind(w)

       elseif ( s%x_character_ctrl(6)=='TSK16' ) then

         call eval_TSK16_wind(w)

       elseif ( s%x_character_ctrl(6)=='Y17' ) then

         call eval_Yoon17_wind(w)

       elseif ( s%x_character_ctrl(6)=='Sh19' ) then

         call eval_Shenar19_wind(w)

       elseif ( s%x_character_ctrl(6)=='S19' ) then

         call eval_Sander19_wind(w)

       elseif( s%x_character_ctrl(6)=='Y06') then

         call eval_Yoon06_wind(w)

      elseif( s%x_character_ctrl(6)=='NL00') then

        call eval_NugisLamers_wind(w)

      elseif ( s%x_character_ctrl(6)=='Vb98' ) then

        call eval_Vanbeveren98_wind(w,"WR")

      elseif ( s%x_character_ctrl(6)=='L89' ) then

        call eval_Langer89_wind(w)

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
       if ( Tsurf > 4000 ) then
         log10w = 2.49536114*log10(L/Lsun) -33.70573385*log10(Tsurf/4000)-19.10219818
       else
         log10w = 2.49536114*log10(L/Lsun) -19.10219818
       end if
     end if

     w = exp10(log10w)

     call smooth_wind_log(w, "A24")

  end subroutine eval_Antoniadis24_wind

  subroutine eval_van_Loon05_wind(w)
     real(dp), intent(out) :: w
     real(dp) :: log10w
     include 'formats'
     wind_scheme = 31.0

     log10w = 1.05 * log10(L/Lsun/1d4) - 6.3 * log10(Tsurf/3500) - 5.56
     w = exp10(log10w)
     call smooth_wind_log(w, "vL05")
  end subroutine eval_van_Loon05_wind

  subroutine eval_Beasor23_wind(w)
     real(dp), intent(out) :: w
     real(dp) :: log10w
     include 'formats'
     wind_scheme = 32.0

     log10w = -0.15 * s% initial_mass + 3.6 * log10(L/Lsun) - 21.5
     w = exp10(log10w)
     call smooth_wind_log(w, "Be23")
  end subroutine eval_Beasor23_wind

  subroutine eval_Yang23_wind(w)
     real(dp), intent(out) :: w
     real(dp) :: log10w
     include 'formats'
     wind_scheme = 33.0

     log10w = 0.45 * log10(L/Lsun)**3 - 5.26 * (log10(L/Lsun))**2 + 20.93 * log10(L/Lsun) - 34.56
     w = exp10(log10w)
     call smooth_wind_log(w, "Ya23")
  end subroutine eval_Yang23_wind

  subroutine eval_Decin24_wind(w)
     real(dp), intent(out) :: w
     real(dp) :: log10w
     include 'formats'
     wind_scheme = 35.0

     log10w = 1.71 -1.63 * s% initial_mass/10 + 3.47 * log10(L/Lsun/1d5) - 5
     w = exp10(log10w)
     call smooth_wind_log(w, "D24")
  end subroutine eval_Decin24_wind

  subroutine eval_de_Jager88_wind(w)
     ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
     real(dp), intent(out) :: w
     real(dp) :: log10w
     include 'formats'
     wind_scheme = 30.0

     log10w = 1.769d0*log10(L/Lsun) - 1.676d0*log10(Tsurf) - 8.158d0
     w = exp10(log10w)
     call smooth_wind_log(w, "dJ88")
  end subroutine eval_de_Jager88_wind

  subroutine eval_cool_wind(w)
    real(dp), intent(inout) :: w
    integer :: which_cool                                     ! Which cool supergiant wind to use (in case of multple initiations)
    real(dp) :: logMdot

    include 'formats'

    if (s% x_ctrl(7) /= s% x_ctrl(8) .and. Tsurf < s% x_ctrl(8) .and. &   ! Second switch to cool supergiant winds. For now these recipes are the only two implemented
       log10(L/Lsun) <= 5.8) then
      which_cool = 7
    else
      which_cool = 8
    end if

    if ( s%x_character_ctrl(which_cool)=='dJ88' ) then

      call eval_de_Jager88_wind(w)

    elseif ( s%x_character_ctrl(which_cool)=='Vb98' ) then

      call eval_Vanbeveren98_wind(w,"cool")

    elseif ( s%x_character_ctrl(which_cool)=='vL05' ) then

      call eval_van_Loon05_wind(w)

    elseif ( s%x_character_ctrl(which_cool)=='Be23' ) then

      call eval_Beasor23_wind(w)

    elseif ( s%x_character_ctrl(which_cool)=='Ya23' ) then

      call eval_Yang23_wind(w)

    elseif ( s%x_character_ctrl(which_cool)=='A24' ) then

      call eval_Antoniadis24_wind(w)

    elseif ( s%x_character_ctrl(which_cool)=='D24' ) then

      call eval_Decin24_wind(w)


    end if


    ! w = w * scaling_factor ! No need for the scaling factor now

  end subroutine eval_cool_wind

  subroutine eval_Hurley00_wLBV(w_standard, w)
    real(dp), intent(in) :: w_standard
    real(dp), intent(out) :: w
    real(dp) :: Mdot_LBV

    ! H00 Eq: Mdot = 0.1 * (10^-5 * R * L^0.5 - 1.0)^3 * (L/6e5 - 1.0)
    ! Conditions: L > 6e5 and 10^-5 * R * L^0.5 > 1.0
    ! Units: R and L are solar units in the formula context (paper says "stellar luminosity, radius... numerical values... in solar units").
    wind_scheme = 90.5
    Mdot_LBV = 0.1d0 * (1.0d-5 * (R/Rsun) * sqrt(L/Lsun) - 1.0d0)**3 * &
               ((L/Lsun)/6.0d5 - 1.0d0)
    ! Convert to Msun/yr. The paper formula result is in Msun/yr.
    w = w_standard + Mdot_LBV
    call smooth_wind_log(w, "H00 LBV")

   return
 end subroutine eval_Hurley00_wLBV

  subroutine eval_Belczynski10_wLBV(w)
     ! Belczynski+2010 LBV2 winds (eq. 8) with factor 1.5
     real(dp), intent(out) :: w
     real(dp) :: log10w
     include 'formats'
     wind_scheme = 91.0

     w  = 1.5d-4
     call smooth_wind_log(w, "Bk10 LBV")
     ! s% max_years_for_timestep = 1d2
  end subroutine eval_Belczynski10_wLBV


  subroutine eval_Cheng24_wLBV(w_standard, w)
    ! Implementation of the Eruptive Mass Loss Model from Cheng et al. (2024)
    ! Logic: Identifies regions where local super-Eddington luminosity creates
    ! excess energy greater than the binding energy of the overlying envelope.

    real(dp), intent(in) :: w_standard
    real(dp), intent(out) :: w

    integer :: k, taucrit_idx, mlost_loc
    real(dp) :: efficiency_xi
    real(dp) :: c_sound_local
    real(dp) :: mlost_final, tdyn_final

    ! Arrays for profile calculations
    real(dp), allocatable :: t_dyn(:), L_excess(:), E_excess_local(:)
    real(dp), allocatable :: E_excess_cumulative(:), E_binding_cumulative(:)
    real(dp), allocatable :: M_above(:), E_diff(:), potential_Mloss(:)
    real(dp) :: integral_excess, integral_mass_denom


    w = 0.0d0
    wind_scheme = 91.5

    ! -----------------------------------------------------------
    ! 1. Configuration
    ! -----------------------------------------------------------
    ! Efficiency parameter xi (Cheng et al. 2024 uses 0.1 to 1.0)
    efficiency_xi = 1.0d0

    ! Allocate arrays based on number of zones
    allocate(t_dyn(s% nz), L_excess(s% nz), E_excess_local(s% nz))
    allocate(E_excess_cumulative(s% nz), E_binding_cumulative(s% nz))
    allocate(M_above(s% nz), E_diff(s% nz), potential_Mloss(s% nz))

    ! Initialize
    E_excess_cumulative = 0d0
    E_binding_cumulative = 0d0
    M_above = 0d0
    E_diff = 0d0
    potential_Mloss = 0d0

    ! -----------------------------------------------------------
    ! 2. Identify Critical Optical Depth (tau_crit)
    ! -----------------------------------------------------------
    ! Find index where tau < c/cs (Definition of inefficient convection region)
    ! We search from surface (k=1) inwards.
    taucrit_idx = 1
    do k = 1, s% nz
       c_sound_local = s% csound(k) ! Adiabatic sound speed usually sufficient
       if (s% tau(k) > (clight / c_sound_local)) then
          taucrit_idx = k
          exit
       end if
    end do

    ! -----------------------------------------------------------
    ! 3. Calculate Local Energetics (Eq 1-6 in Cheng+24)
    ! -----------------------------------------------------------
    do k = 1, s% nz
       ! Local dynamical time: t_dyn = sqrt(r^3 / Gm)
       t_dyn(k) = sqrt(s% r(k)**3 / (standard_cgrav * s% m(k)))

       ! Eddington Luminosity: L_edd = 4*pi*c*G*M / kappa
       ! Excess Luminosity: L_excess = (L_rad) - L_edd
       ! Note: L_rad = L_tot - L_conv
       L_excess(k) = (s% L(k) - s% L_conv(k)) - &
                     (4.0d0 * pi * clight * s% cgrav(k) * s% m_grav(k) / s% opacity(k))

       ! Excess Energy: E = L_excess * t_dyn
       E_excess_local(k) = L_excess(k) * t_dyn(k)

       ! Masking:
       ! 1. Must be super-Eddington (E > 0)
       ! 2. Must be above tau_crit (k < taucrit_idx)
       if (E_excess_local(k) < 0.0d0 .or. k >= taucrit_idx) then
          E_excess_local(k) = 0.0d0
       end if
    end do

    ! -----------------------------------------------------------
    ! 4. Integration (Surface Inwards)
    ! -----------------------------------------------------------
    ! We integrate from surface (k=1) down to current depth k

    integral_excess = 0.0d0
    integral_mass_denom = 0.0d0

    do k = 1, s% nz - 1

       ! Cumulative Binding Energy of envelope above shell k (Eq 8)
       ! Integrating G*m/r * dm
       E_binding_cumulative(k) = E_binding_cumulative(max(1, k-1)) + &
            (standard_cgrav * s% m(k) * s% dm(k) / s% r(k))

       ! Cumulative Excess Energy (Weighted Average) (Eq 7)
       integral_excess = integral_excess + (E_excess_local(k) * s% dm(k))
       integral_mass_denom = integral_mass_denom + s% dm(k)

       if (integral_mass_denom > 0.0d0) then
          E_excess_cumulative(k) = integral_excess / integral_mass_denom
       else
          E_excess_cumulative(k) = 0.0d0
       end if

       ! Mass above current radius (in Solar Masses)
       M_above(k) = (s% m(1) - s% m(k)) / Msun

       ! Energy Difference (Eq 9)
       E_diff(k) = E_excess_cumulative(k) - E_binding_cumulative(k)

       ! Filter for Mass Loss
       ! Logic: If E_diff > 0, this depth *could* be the detachment point.
       if (E_diff(k) > 0.0d0 .and. k < taucrit_idx) then
          potential_Mloss(k) = M_above(k)
       else
          potential_Mloss(k) = 0.0d0
       end if

    end do

    ! -----------------------------------------------------------
    ! 5. Determine Final Mass Loss Rate
    ! -----------------------------------------------------------
    ! Find the innermost radius (maximum mass above) where condition is met
    mlost_final = maxval(potential_Mloss)

    if (mlost_final > 0.0d0) then
       ! Find location of this maximum
       do k = 1, s% nz
          if (potential_Mloss(k) == mlost_final) then
             mlost_loc = k
             exit
          end if
       end do

       ! Get dynamical time at the detachment point (in years)
       tdyn_final = t_dyn(mlost_loc) / secyer

       ! Rate Eq (12): M_dot = xi * DeltaM / t_dyn
       w = efficiency_xi * (mlost_final / tdyn_final)

    else
       w = 0.0d0
    end if

    ! Cleanup
    deallocate(t_dyn, L_excess, E_excess_local)
    deallocate(E_excess_cumulative, E_binding_cumulative)
    deallocate(M_above, E_diff, potential_Mloss)

    w = w_standard + w
    call smooth_wind_log(w, "Ch24 LBV")

 end subroutine eval_Cheng24_wLBV

 subroutine eval_Pauli26_wLBV(w_standard, w)
    ! Pauli+2026 LBV winds
    real(dp), intent(in) :: w_standard
    real(dp), intent(out) :: w
    real(dp) :: log10w,Mdot_inflation
    real(dp) :: R_core,Rcrit_scaling,beta, minRho
    integer :: i ! Explicit declaration for loop variable
    include 'formats'
    wind_scheme = 92.0

    ! check if the envelope is not from a RSG (Pgas is dominating in their envelopes)
    ! Pgas/Ptotal < 0.3 (value arbitrarily chosen from models)

    R_core = s% r(1)
    do i = 1, s% nz, 1
      if((s% Pgas(i)/s% Peos(i)) > 0.15) then
        R_core = s% r(i)
      end if
    end do

    ! Minimum envelope density

    minRho = 1e99
    do i = s % nz, 1, -1
       if(s% rho(i) .le. minRho) then
          minRho = s% rho(i)
       end if
    end do

    ! Petrovic+2006
    Mdot_inflation = 4 * pi * R_core * R_core * minRho *sqrt(standard_cgrav * s% star_mass * Msun / R_core) / (Msun / 31536000d0)
    ! The 31536000d0 part is to convert to years

    ! The model does not apply this mass loss immediately.
    ! It defines a threshold for how much inflation (R​/Rcore​) is allowed before the wind turns on.
    ! This threshold is dynamic and depends on the surface Hydrogen abundance (Xs​).
    Rcrit_scaling = 1d0
    if (X > 0.2d0) then
       Rcrit_scaling = 1d0
    else if (X < 0.2d0 .and. X > 0d0) then
       Rcrit_scaling = Rcrit_scaling + (1-Rcrit_scaling) * X/0.2
    endif

    ! Linear interpolation to smooth and avoid numerical issues
    beta = 0d0
    if (R / R_core > Rcrit_scaling*2) then
       beta = 1d0
    else if (R / R_core > Rcrit_scaling*1.9 .and. R / R_core < Rcrit_scaling*s% x_ctrl(5)) then
       beta = (R/R_core - Rcrit_scaling*1.9)/(Rcrit_scaling*2-Rcrit_scaling*1.9)
    endif

    w = w_standard + beta*Mdot_inflation

    ! From Pauli+ (2026) run_star_extras file:
    !   "soften change in wind to avoid things going bad
    !    I allow a rapid increase in Mdot to properly get the eruptive phase.
    !    to avoid crazy jumps I slowly go down with Mdot. This has only marginal
    !    impact on the final mass but helps with convergence"
    ! This is now handled by the generalized smoothing routine
    call smooth_wind_log(w, "P26 LBV")

 end subroutine eval_Pauli26_wLBV

 subroutine eval_LBV_winds(w)
           real(dp), intent(inout) :: w
           real(dp) :: w_standard
           w_standard = w

           if (s% x_character_ctrl(11) == 'Bk10') then
             if (log10(L/Lsun) > s% x_ctrl(6) .and. 1.0d-5 * R/Rsun * (L/Lsun)**0.5d0 > 1.0d0) then
               call eval_Belczynski10_wLBV(w)
             end if

           elseif (s% x_character_ctrl(11) == 'H00') then
              if (log10(L/Lsun) > s% x_ctrl(6) .and. 1.0d-5 * R/Rsun * (L/Lsun)**0.5d0 > 1.0d0) then
                 call eval_Hurley00_wLBV(w_standard, w)
              end if

            elseif (s% x_character_ctrl(11) == 'Ch24') then
               call eval_Cheng24_wLBV(w_standard, w)

           elseif (s% x_character_ctrl(11) == 'P26') then
              call eval_Pauli26_wLBV(w_standard, w)

           end if
      end subroutine eval_LBV_winds


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
           s% delta_lgT_cntr_limit = 0.05d0                                        ! Stars that enter the RGB need a bit more flexibility for the central temperature
         end if

         extras_check_model = keep_going

        if (s% center_h1 < 1d-3 .and. s% center_he4 < 1d-3 .and. s% center_c12 < 1d-2 &
        .and. s% center_o16 < s% x_ctrl(1) .and. s% x_logical_ctrl(1)) then
         write (*,*) "End of core-O burning. We are happy with that -- Terminate"
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

         if (s% center_c12<=1d-2 .and. s% center_he4<=1d-3 .and. s% center_h1<=1d-3) then


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
