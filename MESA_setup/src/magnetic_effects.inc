subroutine magnetic_evolution(id,ierr,R_cgs,Z_ratio,B0,B_eq)

  integer, intent(in) :: id
  real(dp), intent(in) :: B0, R_cgs, Z_ratio
  real(dp), intent(out) :: B_eq
  real(dp) :: Deceff, t_ms_calc
  integer :: evolution_id
  character(len=32) :: evolution_str

  integer, intent(out) :: ierr
  type (star_info), pointer :: s

  ierr = 0
  call star_ptr(id, s, ierr)
  if (ierr /= 0) return

  Deceff = s% x_ctrl(12)                                                      ! Decay efficiency (for DEC model)


  if (s% model_number <= 10) then                                             ! WATCH OUT FOR THIS ONE!!!!
      R_init = s% photosphere_r * Rsun                                        ! I didn't set model_number == 1 just in case of relaxation
      ! Calculate MS lifetime from paper (Hurley et al. 2000)
      ! mass in Msun, metallicity Z
      call calc_ms_lifetime(s% star_mass, Z_ratio, t_ms_calc)
  end if

  ! Determine field evolution scenario
  evolution_str = trim(adjustl(s% x_character_ctrl(13)))
  select case (evolution_str)
  case ('const')
      evolution_id = 1
  case ('flux')
      evolution_id = 2
  case ('decay')
      evolution_id = 3
  case default
      evolution_id = 1 ! Default to constant
  end select

  select case (evolution_id)
    case (1) ! CONSTANT
      B_eq = B0
    case (2) ! FLUX CONSERVATION (B \propto R^-2)
      if (R_init > 0.0d0) then
        B_eq = B0 * (R_init / R_cgs)**2
      else
        B_eq = B0
      end if
    case (3) ! DECAY
      ! Scaled to an assumed MS age
      B_eq = B0 * exp(-Deceff * s% star_age / t_ms_calc)
  end select

end subroutine magnetic_evolution

subroutine magnetic_braking(id, ierr)
    integer, intent(in) :: id
    real(dp) :: R_cgs, M_cgs, Mdot_cgs, vterm_cgs, Omega, scale, J_ex, J_lost
    real(dp) :: Z_ratio, X_surf, L_cgs, T_surf, Beff, Deceff, B0, t_ms_calc, B_eq
    integer :: k, k_max, torque_id
    character(len=32) :: method_str

    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! DEBUG PRINTS
    ! write(*,*) "Logical Ctrl 9:", s% x_logical_ctrl(9)
    ! write(*,*) "Doing Relax:", s% doing_relax
    ! write(*,*) "Mdot:", s% mstar_dot
    ! write(*,*) "Omega:", s% omega(1)


    ! Exit if magnetic braking is not enabled or if relaxing
    ! if (s% doing_relax) return

    B0 = s% x_ctrl(10)                                                          ! Initial Equatorial Magnetic Field (Gauss)
    Beff = s% x_ctrl(11)                                                        ! Efficiency factor for torque
    if (Beff == 0) Beff = 1.0d0

    ! Determine torque method from character string (x_character_ctrl(12))
    method_str = trim(adjustl(s% x_character_ctrl(12)))
    select case (method_str)
    case ('uniform')
        torque_id = UNIFORM_ID
    case ('surface')
        torque_id = SURFACE_ID
    case default
        torque_id = SURFACE_ID ! Default fallback
    end select


    s% extra_jdot(:) = 0d0
    s% extra_omegadot(:) = 0d0

    ! Skip if no mass loss or no rotation
    if (s% mstar_dot >= 0 .or. s% omega(1) == 0) return

    M_cgs = s% star_mass * Msun
    R_cgs = s% photosphere_r * Rsun
    L_cgs = s% photosphere_L * Lsun
    T_surf = s% Teff
    X_surf = s% x(s% nz)
    Z_ratio = s% kap_rq% Zbase / s% x_ctrl(9)

    call magnetic_evolution(id,ierr,R_cgs,Z_ratio,B0,B_eq)
    call calc_vterm(M_cgs, R_cgs, L_cgs, T_surf, X_surf, Z_ratio, vterm)

    vterm_cgs = vterm * 1d5
    Mdot_cgs = abs(s% mstar_dot)
    Omega = s% omega(1)

    ! write(*,*) "vterm_cgs TEST =", vterm_cgs


    ! 1. Magnetic Confinement Parameter (eta_star)
    etastar = (B_eq * R_cgs)**2 / (Mdot_cgs * vterm_cgs)

    ! 2. Alfven Radius (Ra) - Depends on geometry (Dipole for Uniform, Quad for Surf)
    if (torque_id == UNIFORM_ID) then
        Ra = 1.0d0 + (etastar + 0.25d0)**(1d0/4d0) - 0.25d0**(1d0/4d0)
    else
        Ra = 1.0d0 + (etastar + 0.25d0)**(1d0/6d0) - 0.25d0**(1d0/6d0)
    end if
    Ra = Ra


    ! 3. Closure and Kepler Radii
    Rc = 1.0d0 + 0.7d0 * (Ra - 1.0d0)
    Rk = (s% v_rot_avg_surf / sqrt(standard_cgrav * M_cgs / R_cgs))**(-2.d0/3.d0)

    ! 4. Total Rate of Angular Momentum Loss
    Jbrake = (2.0d0/3.0d0) * Mdot_cgs * Omega * (Ra * R_cgs)**2

    J_lost = abs(s% angular_momentum_removed)
    J_ex = max(Jbrake * s% dt - J_lost, 0.0d0)

    ! write(*,*) "J_lost:", J_lost
    ! write(*,*) "J_ex:", J_ex
    ! write(*,*) "Mdot_cgs:", Mdot_cgs
    ! write(*,*) "Omega:", Omega
    ! write(*,*) "R_cgs:", R_cgs


    ! --- Torque Application ---
    if (J_ex > 0) then
        select case (torque_id)
        case (UNIFORM_ID)
            scale = Beff * J_ex / s% total_angular_momentum
            ! write(*,*) "scale:", scale

            if (scale < 1.0d0 .and. scale > 1d-15) then
                do k = 1, s% nz
                    s% extra_jdot(k) = - s% j_rot(k) * scale / s% dt
                end do
            end if

        case (SURFACE_ID)
            k_max = 1
            do k = 1, s% nz
                if (s% q(k) <= 0.80d0) then
                    k_max = k
                    exit
                end if
            end do
            nin = k_max
            enclosedm = sum(s% dm_bar(1:k_max)) / Msun

            J_surf = dot_product(s% dm_bar(1:k_max), s% j_rot(1:k_max))
            scale = Beff * J_ex / J_surf

            if (scale < 1.0d0 .and. scale > 1d-15) then
                do k = 1, k_max
                    s% extra_jdot(k) = - s% j_rot(k) * scale / s% dt
                end do
            end if
        end select
    end if

end subroutine magnetic_braking

! Calculate MS lifetime following Hurley, Pols, Tout (2000)
subroutine calc_ms_lifetime(M, Z_ratio, t_ms)
    real(dp), intent(in) :: M, Z_ratio
    real(dp), intent(out) :: t_ms
    real(dp) :: zeta, t_bgb, t_hook, x, mu
    real(dp) :: a(10)

    zeta = Z_ratio

    ! Z-dependent coefficients for t_bgb (Eq. 4 and Appendix A)
    ! Note: These are specific combinations for t_bgb fit
    a(1) = 1.593890d3 + 2.053038d3*zeta + 1.231226d3*zeta**2 + 2.327785d2*zeta**3
    a(2) = 2.706708d3 + 1.483131d3*zeta + 7.411230d1*zeta**2 + 5.772723d2*zeta**3
    a(3) = 1.466143d2 - 1.048442d2*zeta - 6.795374d1*zeta**2 - 1.391127d1*zeta**3
    a(4) = 4.141960d-2 + 4.564888d-2*zeta + 2.958542d-2*zeta**2 + 5.571483d-3*zeta**3
    a(5) = 3.426349d-1

    ! Coefficients for MS hook and lifetime ratios (Eq. 6-7)
    a(6) = 1.949814d1 + 1.758178d0*zeta - 6.008212d0*zeta**2 - 4.470533d0*zeta**3
    a(7) = 4.903830d0
    a(8) = 3.166411d-2 - 2.750074d-3*zeta - 2.271549d-3*zeta**2
    a(9) = 5.212154d-2 + 1.312179d0*zeta + 8.073972d-1*zeta**2
    a(10) = 3.294936d-1 + 9.231860d-2*zeta + 2.610989d-2*zeta**2

    ! Time taken to reach base of giant branch (Eq. 4)
    t_bgb = (a(1) + a(2)*M**4.0d0 + a(3)*M**5.5d0 + M**7.0d0) / &
            (a(4)*M**2.0d0 + a(5)*M**7.0d0)

    ! Fractional timescale parameters (Eq. 6-7)
    x = max(0.95d0, min(0.95d0 - 0.03d0*(zeta + 0.30103d0), 0.99d0))
    mu = max(0.5d0, 1.0d0 - 0.01d0 * max(a(6)/(M**a(7)), a(8) + a(9)/(M**a(10))))

    ! Final MS lifetime calculation (Eq. 5); multiplied by 1d6 to make it in [yr]
    t_ms = x * t_bgb *1.0d6

end subroutine calc_ms_lifetime
