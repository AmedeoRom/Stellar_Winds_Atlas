! Parser that checks history_columns.list with fallback directories
      subroutine init_tracked_isos()
         integer :: iunit, iostat, i_comment, cid, env_stat, i
         character(len=256) :: line, mesa_dir
         character(len=32) :: iso_name

         if (is_tracked_isos_init) return

         tracked_iso_cumulative_yield = 0.0d0
         tracked_iso_step_yield = 0.0d0
         num_tracked_isos = 0

         ! 1. Try local simulation directory
         open(newunit=iunit, file='history_columns.list', status='old', iostat=iostat)

         ! 2. If failed, try parent directory (e.g. if run from src/)
         if (iostat /= 0) then
            open(newunit=iunit, file='../history_columns.list', status='old', iostat=iostat)
         end if

         ! 3. If still failed, try MESA default directory
         if (iostat /= 0) then
            call get_environment_variable("MESA_DIR", mesa_dir, status=env_stat)
            if (env_stat == 0 .and. len_trim(mesa_dir) > 0) then
               open(newunit=iunit, file=trim(mesa_dir)//'/star/defaults/history_columns.list', &
                    status='old', iostat=iostat)
            end if
         end if

         if (iostat == 0) then
            write(*,*) "DynamicYields: Parsing history_columns.list..."
            do
               read(iunit, '(A)', iostat=iostat) line
               if (iostat /= 0) exit

               ! Sanitize: replace tabs (ASCII 9) and Windows CR (ASCII 13) with spaces
               do i = 1, len_trim(line)
                  if (iachar(line(i:i)) == 9 .or. iachar(line(i:i)) == 13) line(i:i) = ' '
               end do

               line = adjustl(line)
               ! Ignore empty lines or comments
               if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

               iso_name = ''
               if (line(1:8) == 'surface ') then
                  iso_name = adjustl(line(9:))
               elseif (line(1:8) == 'surface_') then
                  iso_name = adjustl(line(9:))
               end if

               if (len_trim(iso_name) > 0) then
                  ! Strip out inline comments
                  i_comment = index(iso_name, '!')
                  if (i_comment > 0) iso_name = iso_name(1:i_comment-1)
                  iso_name = trim(adjustl(iso_name))

                  ! Force lowercase to ensure map_string_to_chem_id works perfectly
                  do i = 1, len_trim(iso_name)
                     if (iachar(iso_name(i:i)) >= 65 .and. iachar(iso_name(i:i)) <= 90) then
                        iso_name(i:i) = achar(iachar(iso_name(i:i)) + 32)
                     end if
                  end do

                  ! Map string to official chem_def parameter constants as per EOS docs
                  cid = map_string_to_chem_id(iso_name)

                  if (cid > 0) then
                     num_tracked_isos = num_tracked_isos + 1
                     tracked_iso_names(num_tracked_isos) = iso_name
                     tracked_iso_net_ids(num_tracked_isos) = cid ! Store global chem ID

                     write(*,*) "  -> Tracking: ", trim(iso_name), " (Chem ID: ", cid, ")"
                     if (num_tracked_isos == max_tracked_isos) exit
                  end if
               end if
            end do
            close(iunit)
            write(*,*) "DynamicYields: Total tracked isotopes = ", num_tracked_isos
         else
            write(*,*) "Warning: Could not locate a history_columns.list file for dynamic yield tracking."
         end if

         is_tracked_isos_init = .true.
      end subroutine init_tracked_isos


      ! Calculator for mass loss yield at each timestep
      subroutine update_tracked_yields(s)
         type (star_info), pointer :: s
         real(dp) :: dm_lost, mass_to_strip, cell_mass, stripped_from_cell
         integer :: k, cid, local_id, i, i_spec

         ! s% mstar_dot is internally in grams/sec, s% dt is in seconds.
         if (s% mstar_dot < 0.0d0) then
            dm_lost = (abs(s% mstar_dot) * s% dt) / Msun
         else
            dm_lost = 0.0d0
         end if

         ! Initialize step yields to zero for the accumulation
         tracked_iso_step_yield = 0.0d0

         if (dm_lost > 0.0d0) then
            mass_to_strip = dm_lost

            ! Dig deeper into the grid starting from the surface (i=1)
            do i = 1, s% nz
               ! s% dm(i) evaluates the cell mass perfectly in grams.
               cell_mass = s% dm(i) / Msun

               if (cell_mass <= mass_to_strip) then
                  ! The entire cell is ejected
                  stripped_from_cell = cell_mass
                  mass_to_strip = mass_to_strip - cell_mass
               else
                  ! Only a fraction of this cell is ejected
                  stripped_from_cell = mass_to_strip
                  mass_to_strip = 0.0d0
               end if

               ! Add this cell's contribution to each tracked isotope (LINEAR ACCUMULATION)
               do k = 1, num_tracked_isos
                  cid = tracked_iso_net_ids(k)

                  ! ROBUST MAPPING: Search the current network for the global ID
                  local_id = 0
                  do i_spec = 1, s% species
                     if (s% chem_id(i_spec) == cid) then
                        local_id = i_spec
                        exit
                     end if
                  end do

                  if (local_id > 0) then
                     tracked_iso_step_yield(k) = tracked_iso_step_yield(k) + &
                                                 (stripped_from_cell * s% xa(local_id, i))
                  end if
               end do

               ! Stop digging if we've found enough mass
               if (mass_to_strip <= 0.0d0) exit
            end do
         end if

         ! AFTER linear accumulation is completely finished, safely convert to log10
         do k = 1, num_tracked_isos
            if (tracked_iso_step_yield(k) > 1.0d-99) then
               tracked_iso_step_yield(k) = log10(tracked_iso_step_yield(k))
            else
               tracked_iso_step_yield(k) = -99.0d0 ! Safe floor if mass lost is exactly 0
            end if

            ! Unconditional print statement OUTSIDE the digging loop
            if (s% mstar_dot > 0.0d0) then
              write(*,*) "Timestep Yield [log10(Msun)] for ", trim(tracked_iso_names(k)), " : ", tracked_iso_step_yield(k)
            endif
         end do

         ! Cumulative yields kept commented out per user request
         ! do k = 1, num_tracked_isos
         !    tracked_iso_cumulative_yield(k) = tracked_iso_cumulative_yield(k) + tracked_iso_step_yield(k)
         ! end do
      end subroutine update_tracked_yields


      integer function map_string_to_chem_id(iso_name)
         character(len=*), intent(in) :: iso_name
         select case(trim(iso_name))
         case('h1'); map_string_to_chem_id = ih1
         case('h2'); map_string_to_chem_id = ih2
         case('he3'); map_string_to_chem_id = ihe3
         case('he4'); map_string_to_chem_id = ihe4
         case('li7'); map_string_to_chem_id = ili7
         case('be9'); map_string_to_chem_id = ibe9
         case('b11'); map_string_to_chem_id = ib11
         case('c12'); map_string_to_chem_id = ic12
         case('c13'); map_string_to_chem_id = ic13
         case('n14'); map_string_to_chem_id = in14
         case('n15'); map_string_to_chem_id = in15
         case('o16'); map_string_to_chem_id = io16
         case('o17'); map_string_to_chem_id = io17
         case('o18'); map_string_to_chem_id = io18
         case('f19'); map_string_to_chem_id = if19
         case('ne20'); map_string_to_chem_id = ine20
         case('ne21'); map_string_to_chem_id = ine21
         case('ne22'); map_string_to_chem_id = ine22
         case('na23'); map_string_to_chem_id = ina23
         case('mg24'); map_string_to_chem_id = img24
         case('mg25'); map_string_to_chem_id = img25
         case('mg26'); map_string_to_chem_id = img26
         case('al27'); map_string_to_chem_id = ial27
         case('si28'); map_string_to_chem_id = isi28
         case('p31'); map_string_to_chem_id = ip31
         case('s32'); map_string_to_chem_id = is32
         case('ar36'); map_string_to_chem_id = iar36
         case('ca40'); map_string_to_chem_id = ica40
         case('ti44'); map_string_to_chem_id = iti44
         case('cr48'); map_string_to_chem_id = icr48
         case('fe52'); map_string_to_chem_id = ife52
         case('fe54'); map_string_to_chem_id = ife54
         case('fe56'); map_string_to_chem_id = ife56
         case('co56'); map_string_to_chem_id = ico56
         case('ni56'); map_string_to_chem_id = ini56
         case default; map_string_to_chem_id = 0
         end select
      end function map_string_to_chem_id
