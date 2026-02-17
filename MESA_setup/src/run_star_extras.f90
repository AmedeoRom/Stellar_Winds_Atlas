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

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use binary_def
      use utils_lib



      implicit none

      logical :: already_thick
      real(dp) :: Mdot_switch,L_switch,M_switch,gamma_edd_switch
      real(dp) :: eta,eta_trans,gamma_edd,gamma_edd_old,vterm                   ! gamma_edd_old is to check the previous timestep
      real(dp) :: wind_scheme,wind_scheme_interp                                ! To know which winds model I am using at each timestep

      ! this is used to soften too large changes in wind mass loss rates
      real(dp) :: old_wind = 0
      character(len=10) :: vterm_type = "Hawcroft24"                            ! Options Lamers95, Crowther06, Hawcroft24, Kritcka25


      real(dp) :: max_years_dt_old
      ! these routines are called by the standard run_star check_model


      ! --- Magnetic Variables from Keszthelyi et al. prescription ---
      real(dp) :: Beq             ! Equatorial magnetic field strength (G)
      real(dp) :: etastar         ! Equatorial wind magnetic confinement parameter
      real(dp) :: Ra              ! Alfven radius in stellar radii
      real(dp) :: Rc              ! Closure radius in stellar radii
      real(dp) :: R_init = 0.0d0  ! Initial radius for magnetic field evolution
      real(dp) :: Rk              ! Kepler co-rotation radius in stellar radii
      real(dp) :: Jbrake          ! Total rate of angular momentum loss (dJ/dt)
      real(dp) :: enclosedm       ! Mass in the braking region (Msun)
      real(dp) :: J_surf          ! AM in the surface reservoir
      integer  :: nin             ! Index of the last zone for surface torque

      ! Internal IDs for torque logic
      integer, parameter :: UNIFORM_ID = 1
      integer, parameter :: SURFACE_ID = 2

      contains

      include "my_other_wind.inc"
      include "magnetic_effects.inc"


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
         s% other_torque => magnetic_braking

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

     ! Terminal velocity has its own dedicated function because it is needed for both the winds and the magnetic braking module

     subroutine calc_vterm(M, R, L, Tsurf, X, Z_ratio, v_out)
       real(dp), intent(in) :: M, R, L, Tsurf, X, Z_ratio
       real(dp), intent(out) :: v_out

       real(dp) :: log_gamma_edd, gamma_edd, vesc, vterm_e

       ! 1. Calculate Gamma_Edd (copied from my_other_wind)
       log_gamma_edd = -4.813d0 + log10(1d0+X) + log10(L/Lsun) - log10(M/Msun)
       gamma_edd = 10**log_gamma_edd

       ! 2. Calculate Escape Velocity
       vesc = sqrt(2d0*standard_cgrav*M/R)/1d5

       ! 3. Calculate Vterm based on selected type
       select case (trim(vterm_type))
       case ("Lamers95")
         if ( Tsurf > 25000 ) then
           vterm_e = 2.6
         else
           vterm_e = 1.3
         end if
         v_out = vterm_e * sqrt(2d0*standard_cgrav*(M)*(1-gamma_edd)/R)/1d5 * Z_ratio**0.20d0

       case ("Crowther06")
         vterm_e = 1.59d-4*Tsurf - 0.89
         v_out = vterm_e * sqrt(2d0*standard_cgrav*(M)*(1-gamma_edd)/R)/1d5 * Z_ratio**0.20d0

       case ("Hawcroft24")
         v_out = (9.2d-2*Tsurf-1040) * Z_ratio**0.22d0

       case ("Kritcka25")
         v_out = (107+25*log10(Z_ratio))*(Tsurf/1d3) - 1190 - 430*log10(Z_ratio)

       case default
         v_out = 0d0 ! Safety fallback
       end select

     end subroutine calc_vterm


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

         if ( s% x_logical_ctrl(9) ) then
            how_many_extra_history_columns = 17
         else
           how_many_extra_history_columns = 8
         end if

      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         real(dp) :: Msum,Rsum,E_Bind_G,E_Bind_B,E_Bind_H
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
         names(6) = "E_Bind_G"                                                  ! Summary of the different binding energies from Sgalletta+ (2026)
         names(7) = "E_Bind_B"
         names(8) = "E_Bind_H"

         if ( s% x_logical_ctrl(9) ) then                                       ! In case magnetic braking is active

           names(9) = "Beq"; vals(9) = Beq
           names(10) = "Ra"; vals(10) = Ra
           names(11) = "Rc"; vals(11) = Rc
           names(12) = "Rk"; vals(12) = Rk
           names(13) = "etastar"; vals(13) = etastar
           names(14) = "Jbrake"; vals(14) = Jbrake
           names(15) = "vterm_km_s"; vals(15) = vterm
           names(16) = "enclosedm_braking"; vals(16) = enclosedm
           names(17) = "nin_braking"; vals(17) = real(nin, dp)

         end if


         Msum=0
         Rsum=0
         E_Bind_G = 0
         E_Bind_B = 0
         E_Bind_H = 0

         do i=1, s% nz
           if ((s% center_he4 > 1d-3 .and. s% x(i) >= s% he_core_boundary_h1_fraction) .or. &
           (s% center_h1 < 1d-3 .and. s% center_he4 < 1d-3) .and. s% y(i) >= s% co_core_boundary_he4_fraction) then
             E_Bind_G = E_Bind_G + s% dm(i)*(-standard_cgrav* s% m(i)/s% r(i))
             E_Bind_B = E_Bind_B + s% dm(i)*(-standard_cgrav* s% m(i)/s% r(i) + s% energy(i))
             E_Bind_H = E_Bind_H + s% dm(i)*(-standard_cgrav* s% m(i)/s% r(i) + s% energy(i) + (s% prad(i) + s% pgas(i))/s% rho(i))
           else
             exit
           end if
         end do

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
         vals(6) = E_Bind_G
         vals(7) = E_Bind_B
         vals(8) = E_Bind_H
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
