!! Copyright, David Strubbe, 30 Nov 2015
!! Based on octopus 5.0.0, www.tddft.org/programs/octopus
!! Code which is GPL-2 licensed, originally by X. Andrade
!! Files src/grid/io_function_inc.F90, src/basic/openscad.F90

!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.

module bxsf2scad_m

  use polyhedron_m

  implicit none

contains

  !---------------------------------------------------------------------------
  ! bring ii to allowed indexing interval [1, nk], with periodic boundary conditions
  integer function iadjust_bc(ii, nk)
    integer, intent(in) :: ii, nk

    if(ii <= 0) then
      iadjust_bc = mod(ii, nk) + nk
    else if(ii <= nk) then
      iadjust_bc = ii
    else
      iadjust_bc = mod(ii, nk)
    endif
  end function iadjust_bc

  !---------------------------------------------------------------------------
  ! bring ii to allowed indexing interval [1, nk], with periodic boundary conditions
  real*8 function dadjust_bc(xx, nk)
    real*8, intent(in) :: xx
    integer, intent(in) :: nk

    if(xx <= 0) then
      dadjust_bc = xx + nk
    else if(xx <= nk) then
      dadjust_bc = xx
    else
      dadjust_bc = xx - nk
    endif
  end function dadjust_bc
  
  !---------------------------------------------------------------------------
  function interpolate_isolevel(energies, fermi_energy, ik, jk, nk) result(pos)
    real*8,          intent(in) :: energies(:,:,:)
    real*8,          intent(in) :: fermi_energy
    integer,         intent(in) :: ik(3), jk(3), nk(3)
    real*8                      :: pos(1:3)

    real*8 :: v1, v2, x1(3), x2(3), xdiff(3), xsum(3)
    integer :: ii
  
    v1 = energies(ik(1), ik(2), ik(3))
    v2 = energies(jk(1), jk(2), jk(3))
    x1(1:3) = ik(1:3)
    x2(1:3) = jk(1:3)

    if((v1 - fermi_energy)*(v2 - fermi_energy) > 0) then
      write(0,*) 'warning, not straddling'
    endif
    ! must be ef is btw v1, v2
    ! note: pos should be close to x1, so do not impose PBCs on the difference vector

    ! take shortest distance between two points, in interval [-nk/2 + 1, nk/2]
    do ii = 1, 3
!      xdiff(ii) = adjust_bc(jk(ii) - ik(ii), nk(ii))
      xdiff(ii) = jk(ii) - ik(ii)
      if(abs(v2 - v1) > 1d-10) then
        xdiff(ii) = (fermi_energy - v1)*(xdiff(ii))/(v2 - v1)
      else
        xdiff(ii) = xdiff(ii) * 0.5d0
      endif
!      if(xdiff(ii) >  nk(ii)*0.5d0) xdiff(ii) = xdiff(ii) - nk(ii)
!      if(xdiff(ii) < -nk(ii)*0.5d0) xdiff(ii) = xdiff(ii) + nk(ii)
    enddo

!      pos(1:3) = x1(1:3) + (fermi_energy - v1)*(xdiff(1:3))/(v2 - v1)
      ! this should never happen, but just to be sure -- why not?
!      pos(1:3) = xsum(1:3)

    pos(1:3) = x1(1:3) + xdiff(1:3)

!    if(abs(v2 - v1) > 1d-10) then
!      write(0,*) 'x1 = ', ik(1:3), ' x2 = ', jk(1:3)
!      write(0,*) 'diff', xdiff(1:3)
!      write(0,*) 'pos = ', pos(1:3)
!    endif

  end function interpolate_isolevel

  ! build surface LESS THAN Fermi energy, i.e. electron
  !---------------------------------------------------------------------------
  logical function inside_isolevel(energies, ik, nk, fermi_energy)
    real*8, intent(in) :: energies(:,:,:)
    integer, intent(in) :: ik(3), nk(3)
    real*8, intent(in) :: fermi_energy

    inside_isolevel = energies(iadjust_bc(ik(1), nk(1)), iadjust_bc(ik(2), nk(2)), iadjust_bc(ik(3), nk(3))) < fermi_energy
  end function inside_isolevel

  !---------------------------------------------------------------------------
  subroutine openscad_file_polyhedron(iunit, poly)
    integer,            intent(in) :: iunit
    type(polyhedron_t), intent(in) :: poly
    
    integer :: ii, minmap, maxmap, jj
    integer, allocatable :: map(:)

    minmap = minval(poly%point_indices(1:poly%npoints))
    maxmap = maxval(poly%point_indices(1:poly%npoints))

    allocate(map(minmap:maxmap))

    map = -1

    write(iunit, '(a)') ' polyhedron( points = [ '

    jj = 0
    do ii = 1, poly%npoints
      if(map(poly%point_indices(ii)) == -1) then !skip duplicated points
        write(iunit, '(a,f12.6,a,f12.6,a,f12.6,a,i6)') &
          '[', poly%points(1, ii), ',', poly%points(2, ii), ',',  poly%points(3, ii), '], //', jj
        map(poly%point_indices(ii)) = jj
        jj = jj + 1
      end if
    end do

    write(iunit, '(a)') '], faces = [ '

    do ii = 1, poly%ntriangles
      write(iunit, '(a,i10,a,i10,a,i10,a)') &
        '[', map(poly%triangles(3, ii)), ',', map(poly%triangles(2, ii)), ',',  map(poly%triangles(1, ii)), '],'
    end do

    write(iunit, '(a)') ']);'

    deallocate(map)

  end subroutine openscad_file_polyhedron

  !---------------------------------------------------------------------------
  subroutine out_openscad(energies, fermi_energy, nk, bvec)
    real*8,  intent(in) :: energies(:, :, :)
    real*8,  intent(in) :: fermi_energy
    integer, intent(in) :: nk(3)
    real*8,  intent(in) :: bvec(3, 3)

    integer :: ip, ii, jp, jj, kk, ll, mm, npoly, curve_res, iunit_scad
    type(polyhedron_t) :: poly

    integer, allocatable :: edges(:), triangles(:, :)
    integer :: iunit, cubeindex, cube_point(0:7, 1:3)
    real*8 :: vertlist(1:3, 0:11), vertlist_cart(1:3, 0:11), vertlist0_old(1:3), offset(1:3)

    allocate(edges(0:255))
    allocate(triangles(1:16, 0:255))

    iunit = 60
    open(unit = iunit, file= "marching_cubes_edges.data", action='read', status='old')
    do ii = 0, 255
      read(iunit, *) edges(ii)
    end do
    close(iunit)

    open(unit = iunit, file= "marching_cubes_triangles.data", action='read', status='old')
    do ii = 0, 255
      read(iunit, *) (triangles(jj, ii), jj = 1, 16)
    end do
    close(iunit)

    iunit_scad = 61
    open(unit = iunit_scad, file='bxsf.scad', action = 'write', status='replace')
    curve_res = 50
    write(iunit_scad, '(a,i10,a)') '$fn = ', curve_res, ';'

    npoly = 0
    do kk = 1, nk(3)
    do jj = 1, nk(2)
    do ii = 1, nk(1)   

      ! these are allowed to go out of bounds, i.e. be 0 or nk(ii) + 1
      cube_point(0, :) = (/ii    , jj    , kk    /)
      cube_point(1, :) = (/ii    , jj + 1, kk    /)
      cube_point(2, :) = (/ii + 1, jj + 1, kk    /)
      cube_point(3, :) = (/ii + 1, jj    , kk    /)
      cube_point(4, :) = (/ii    , jj    , kk + 1/)
      cube_point(5, :) = (/ii    , jj + 1, kk + 1/)
      cube_point(6, :) = (/ii + 1, jj + 1, kk + 1/)
      cube_point(7, :) = (/ii + 1, jj    , kk + 1/)

      ! all non-straddling warnings seemed to come from when this was not satisfied
!      if(any(cube_point(:,:) <= 0)) cycle
!      do ll = 1, 3
!        if(any(cube_point(:,ll) > nk(ll))) cycle
!      enddo
      
      cubeindex = 0
      if(inside_isolevel(energies, cube_point(0, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 1
      if(inside_isolevel(energies, cube_point(1, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 2
      if(inside_isolevel(energies, cube_point(2, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 4
      if(inside_isolevel(energies, cube_point(3, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 8
      if(inside_isolevel(energies, cube_point(4, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 16
      if(inside_isolevel(energies, cube_point(5, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 32
      if(inside_isolevel(energies, cube_point(6, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 64
      if(inside_isolevel(energies, cube_point(7, :), nk(1:3), fermi_energy)) cubeindex = cubeindex + 128

      if(edges(cubeindex) == 0) cycle
      npoly = npoly + 1
      
      vertlist = 3.333333333333333333d0

      if(iand(edges(cubeindex), 1) /= 0) then
        vertlist(1:3, 0) = interpolate_isolevel(energies, fermi_energy, cube_point(0, :), cube_point(1, :), nk)
      end if
      if(iand(edges(cubeindex), 2) /= 0) then
        vertlist(1:3, 1) = interpolate_isolevel(energies, fermi_energy, cube_point(1, :), cube_point(2, :), nk)
      end if
      if(iand(edges(cubeindex), 4) /= 0) then
        vertlist(1:3, 2) = interpolate_isolevel(energies, fermi_energy, cube_point(2, :), cube_point(3, :), nk)
      end if
      if(iand(edges(cubeindex), 8) /= 0) then
        vertlist(1:3, 3) = interpolate_isolevel(energies, fermi_energy, cube_point(3, :), cube_point(0, :), nk)
      end if
      if(iand(edges(cubeindex), 16) /= 0) then
        vertlist(1:3, 4) = interpolate_isolevel(energies, fermi_energy, cube_point(4, :), cube_point(5, :), nk)
      end if
      if(iand(edges(cubeindex), 32) /= 0) then
        vertlist(1:3, 5) = interpolate_isolevel(energies, fermi_energy, cube_point(5, :), cube_point(6, :), nk)
      end if
      if(iand(edges(cubeindex), 64) /= 0) then
        vertlist(1:3, 6) = interpolate_isolevel(energies, fermi_energy, cube_point(6, :), cube_point(7, :), nk)
      end if
      if(iand(edges(cubeindex), 128) /= 0) then
        vertlist(1:3, 7) = interpolate_isolevel(energies, fermi_energy, cube_point(7, :), cube_point(4, :), nk)
      end if
      if(iand(edges(cubeindex), 256) /= 0) then
        vertlist(1:3, 8) = interpolate_isolevel(energies, fermi_energy, cube_point(0, :), cube_point(4, :), nk)
      end if
      if(iand(edges(cubeindex), 512) /= 0) then
        vertlist(1:3, 9) = interpolate_isolevel(energies, fermi_energy, cube_point(1, :), cube_point(5, :), nk)
      end if
      if(iand(edges(cubeindex), 1024) /= 0) then
        vertlist(1:3, 10) = interpolate_isolevel(energies, fermi_energy, cube_point(2, :), cube_point(6, :), nk)
      end if
      if(iand(edges(cubeindex), 2048) /= 0) then
        vertlist(1:3, 11) = interpolate_isolevel(energies, fermi_energy, cube_point(3, :), cube_point(7, :), nk)
      end if

      ! transform to Cartesian coordinates
      vertlist_cart(1:3, 0:11) = matmul(bvec(1:3, 1:3), vertlist(1:3, 0:11))

      call polyhedron_init(poly)
      ll = 1
      do
        if(triangles(ll, cubeindex) == -1) exit

        call polyhedron_add_point(poly, triangles(ll    , cubeindex), vertlist_cart(1:3, triangles(ll    , cubeindex)))
        call polyhedron_add_point(poly, triangles(ll + 1, cubeindex), vertlist_cart(1:3, triangles(ll + 1, cubeindex)))
        call polyhedron_add_point(poly, triangles(ll + 2, cubeindex), vertlist_cart(1:3, triangles(ll + 2, cubeindex)))
        call polyhedron_add_triangle(poly, triangles(ll:ll + 2, cubeindex))

        ll = ll + 3
      end do

      call openscad_file_polyhedron(iunit_scad, poly)
      call polyhedron_end(poly)

    end do
    enddo
    enddo

    deallocate(edges)
    deallocate(triangles)
    
    close(iunit_scad)
    write(6,'(a,i9,a,a)') ' Wrote ', npoly, ' polyhedra to ', "scad"

    if(npoly == 0) then
      write(0,*) "There were no points inside the isosurface for OpenSCAD output."
    endif

  end subroutine out_openscad

end module bxsf2scad_m

!---------------------------------------------------------------------------
program bxsf2scad

  use bxsf2scad_m

  implicit none

  ! question: how can you close the pores on the edges of the Fermi surface?

  integer :: iunit_bxsf, iunit_scad, nbands, nk(3), ii, iband, ix, iy, iz, ix_, iy_, iz_
  real*8 :: dk(3), bvec(3, 3), fermi_energy, dummy
  real*8, allocatable :: energies(:, :, :), energies_general(:, :, :)
  character*256 :: line, line_trim, word
  logical :: ready
  integer :: center(1:3)

  iunit_bxsf = 70
  open(unit = iunit_bxsf, file = 'bxsf', form = 'formatted', status = 'old', action = 'read')

  ! XCrySDen BXSF format reference:
  ! http://www.xcrysden.org/doc/XSF.html#__toc__14
  ! note: BXSF uses general grid, including edge point twice.
    
  ready = .false.
  do while(.not. ready)
    read(iunit_bxsf, '(a)') line
    line_trim = adjustl(line)
    if(line_trim(1:12) == 'Fermi Energy') then
      read(line_trim, *) word, word, fermi_energy
    endif
    ready = (trim(adjustl(line)) == 'BEGIN_BLOCK_BANDGRID_3D')
  enddo

  read(iunit_bxsf, '(a)') line ! band_energies
  read(iunit_bxsf, '(a)') line ! ..BANDGRID_3D..
  
  read(iunit_bxsf, *) nbands
  read(iunit_bxsf, *) nk(1:3)
  read(iunit_bxsf, *) dk(1:3)
  do ii = 1, 3
    read(iunit_bxsf, *) bvec(1:3, ii)
  enddo

  ! center on Gamma, typically
  center = nk(1:3) / 2 + nint(dk(1:3))
  
  ! option to select band, Fermi level, filenames?
  allocate(energies(nk(1), nk(2), nk(3)))
  allocate(energies_general(nk(1), nk(2), nk(3)))

  write(6,*) 'Fermi energy = ', fermi_energy

  do iband = 1, nbands
    write(6, *) 'Reading band ', iband, ' of ', nbands
    read(iunit_bxsf, '(a)') line ! BAND:  iband

    read(iunit_bxsf, *) (((energies_general(ix, iy, iz), iz=1,nk(3)), iy=1,nk(2)), ix=1,nk(1))

!    mod(ii, nk(i)-1) is in [0, nk-2]
    do ix = 1, nk(1)-1
      ix_ = mod(ix + center(1) - 1, nk(1)-1) + 1
      do iy = 1, nk(2)-1
        iy_ = mod(iy + center(2) - 1, nk(2)-1) + 1
        do iz = 1, nk(3)-1
          iz_ = mod(iz + center(3) - 1, nk(3)-1) + 1
          energies(ix_, iy_, iz_) = energies_general(ix, iy, iz)
        enddo
      enddo
    enddo

    ! restore to general grid, after translation
    energies(nk(1), :, :) = energies(1, :, :)
    energies(:, nk(2), :) = energies(:, 1, :)
    energies(:, :, nk(3)) = energies(:, :, 1)

    write(6,*) 'Minimum energy = ', minval(energies), ' Maximum energy = ', maxval(energies)
    if(fermi_energy < maxval(energies) .and. fermi_energy > minval(energies)) then
      write(6,*) 'Fermi energy is in range of energies for this band. Plotting.'
      ! do not use extra points in general grid, use nk-1
      call out_openscad(energies(:, :, :), fermi_energy, nk - 1, bvec)
      ! FIXME: so far, plotting only first relevant band
      exit
    endif
  enddo

  read(iunit_bxsf, '(a)') line ! END_BANDGRID_3D
  read(iunit_bxsf, '(a)') line ! END_BLOCK_BANDGRID_3D

  close(iunit_bxsf)

  deallocate(energies_general)
  deallocate(energies)

end program bxsf2scad
