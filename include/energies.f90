module energies
    use types
    use box_pot_parameters

    implicit none

    private
    public force, force_ll

contains

    subroutine force(r, f, ener)
        real(dp), intent(in) :: r(:, :)
        real(dp), intent(inout) :: f(:, :)
        real(dp), intent(out) :: ener

        ! Local variables
        integer, parameter :: n = 3 ! Dimensión espacial
        integer :: i, j
        real(dp) :: rij, rposij(n), uij
        real(dp) :: fij(n), fxij, fyij, fzij
        real(dp) :: dnrm2 ! Norma euclidiana de BLAS

        ! Inicializar arreglos y variables
        ener = 0.0_dp
        f = 0.0_dp
        fij = 0.0_dp

        do i = 1, np - 1
            do j = i + 1, np
                uij = 0._dp
                fxij = 0._dp
                fyij = 0._dp
                fzij = 0._dp

                ! xij = x(i)-x(j)
                ! yij = y(i)-y(j)
                ! zij = z(i)-z(j)
                rposij = r(:, i) - r(:, j)

                ! xij = xij-boxl*idnint(xij/boxl)
                ! yij = yij-boxl*idnint(yij/boxl)
                ! zij = zij-boxl*idnint(zij/boxl)
                rposij = rposij - boxl*idnint(rposij/boxl)

                ! rij = norm2( [xij, yij, zij] )
                ! Mejor rendimiento con BLAS
                ! rij = dnrm2( 3,[xij, yij, zij],1 )
                rij = dnrm2(n, rposij, 1)

                if (rij < rc) then
                    call hardsphere(rij, uij, rposij, fxij, fyij, fzij)
                    fij = [fxij, fyij, fzij]

                    ener = ener + uij
                    ! fx(i) = fx(i) + fxij
                    ! fy(i) = fy(i) + fyij
                    ! fz(i) = fz(i) + fzij
                    f(:, i) = f(:, i) + fij

                    ! fx(j) = fx(j) - fxij
                    ! fy(j) = fy(j) - fyij
                    ! fz(j) = fz(j) - fzij
                    f(:, j) = f(:, j) - fij
                end if
            end do
        end do
    end subroutine force

    SUBROUTINE force_ll(rpos, f, ener)
        USE, INTRINSIC :: iso_fortran_env, ONLY: output_unit, error_unit
        USE link_list_module, ONLY: make_list, sc, head, list
        implicit none

        ! Variables de entrada
        real(dp), intent(inout) :: rpos(:, :)
        real(dp), intent(inout) :: f(:, :)
        real(dp), intent(out) :: ener
        ! Variables locales
        INTEGER               :: i, j, ci1, ci2, ci3, k
        INTEGER, DIMENSION(3) :: ci, cj
        REAL(dp) :: r_cut_box, uij
        REAL(dp), DIMENSION(3) :: fij, rposij
        real(dp) :: fxij, fyij, fzij, rij
        real(dp) :: dnrm2 ! Norma euclidiana de BLAS

        ! Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
        ! The cells are chosen so that if (d1,d2,d3) appears, then (-d1,-d2,-d3) does not.
        INTEGER, PARAMETER :: nk = 13
        INTEGER, DIMENSION(3, 0:nk), PARAMETER :: d = RESHAPE([ &
            &                0, 0, 0, 1, 0, 0, &
            &    1, 1, 0, -1, 1, 0, 0, 1, 0, &
            &    0, 0, 1, -1, 0, 1, 1, 0, 1, &
            &   -1, -1, 1, 0, -1, 1, 1, -1, 1, &
            &   -1, 1, 1, 0, 1, 1, 1, 1, 1], [3, nk + 1])

        r_cut_box = rc/boxl

        ! Periodic boundary conditions
        rpos(:, :) = rpos(:, :)/boxl
        rpos(:, :) = rpos(:, :) - dnint(rpos(:, :))
        ! do i = 1, np
        !     print *, rpos(:, i)/boxl
        ! end do
        ! stop
        CALL make_list(np, r_cut_box, rpos(:, :))

        ! Initialize
        f = 0.0_dp

        ! Triple loop over cells
        DO ci1 = 0, sc - 1
            DO ci2 = 0, sc - 1
                DO ci3 = 0, sc - 1

                    ci(:) = [ci1, ci2, ci3] ! 3D index of i-cell
                    i = head(ci1, ci2, ci3)     ! First i-atom in cell

                    DO ! Begin loop over i-atoms in list
                        IF (i == 0) EXIT ! End of link list

                        DO k = 0, nk ! Loop over neighbouring cells

                            IF (k == 0) THEN
                                j = list(i) ! First j-atom is downlist from i in current cell
                            ELSE
                                cj(:) = ci(:) + d(:, k)          ! Neighbour j-cell 3D index
                                cj(:) = MODULO(cj(:), sc)    ! Periodic boundary correction
                                j = head(cj(1), cj(2), cj(3)) ! First j-atom in neighbour cell
                            END IF

                            DO ! Begin loop over j-atoms in list
                                IF (j == 0) EXIT ! End of link list

                                IF (j == i) THEN ! This should never happen
                                    WRITE (unit=error_unit, fmt='(a,2i15)') 'Index error', i, j
                                    STOP 'Impossible error in force'
                                END IF

                                rposij(:) = rpos(:, i) - rpos(:, j)
                                rposij(:) = rposij(:) - dnint(rposij(:))
                                rposij(:) = rposij(:)*boxl
                                rij = dnrm2(3, rposij(:), 1)
                                ! print *, rij

                                IF ((rij > 0.5) .and. (rij < rc)) THEN
                                    ! if (rij < rc) then
                                    ! print *, rij
                                    call hardsphere(rij, uij, rposij(:), fxij, fyij, fzij)
                                    fij = [fxij, fyij, fzij]

                                    ener = ener + uij
                                    ! if (any(fij .ne. 0.0_dp)) print *, fij
                                    ! print *, fij
                                    f(:, i) = f(:, i) + fij
                                    f(:, j) = f(:, j) - fij

                                END IF ! End check within cutoff

                                j = list(j) ! Next j-atom
                            END DO ! End loop over j-atoms in list

                        END DO ! End loop over neighbouring cells

                        i = list(i) ! Next i-atom
                    END DO ! End loop over i-atoms in list

                END DO
            END DO
        END DO

        rpos(:, :) = rpos(:, :)*boxl
    END SUBROUTINE force_ll

    subroutine hardsphere(rij, uij, rposij, fxij, fyij, fzij)
        real(dp), intent(in) :: rij, rposij(:)
        real(dp), intent(inout) :: fxij, fyij, fzij, uij

        ! Local variables
        real(dp) ::  fij

        if (rij < bpot) then
            uij = (a2/dT)*((1.0_dp/rij)**dlr - (1.0_dp/rij)**dla)
            uij = uij + 1.0_dp/dT
            fij = dlr*(1.0_dp/rij)**(dlr + 1.0_dp) - dla*(1.0_dp/rij)**(dla + 1.0_dp)
            fij = a2*fij/dT
        else
            uij = 0.0_dp
            fij = 0.0_dp
        end if

        fxij = fij*rposij(1)/rij
        fyij = fij*rposij(2)/rij
        fzij = fij*rposij(3)/rij
    end subroutine hardsphere
end module energies
