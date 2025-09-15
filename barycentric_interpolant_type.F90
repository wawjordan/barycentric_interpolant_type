module set_precision
  use iso_fortran_env, only : real64
  implicit none
  private
  public :: dp
  integer, parameter :: dp  = real64
end module set_precision

module set_constants
  use set_precision, only : dp
  implicit none
  private
  public :: zero, one, two, three, four
  public :: half, third, fourth, large, pi
  real(dp), parameter :: zero      = 0.0_dp
  real(dp), parameter :: one       = 1.0_dp
  real(dp), parameter :: two       = 2.0_dp
  real(dp), parameter :: three     = 3.0_dp
  real(dp), parameter :: four      = 4.0_dp
  real(dp), parameter :: third     = one / three
  real(dp), parameter :: fourth    = 0.25_dp
  real(dp), parameter :: half      = 0.50_dp
  real(dp), parameter :: large  = huge(one)
  real(dp), parameter :: pi     = acos(-one)
end module set_constants

module project_inputs
  implicit none
  private
  public :: verbose_level
  integer :: verbose_level = 0
end module project_inputs

module message

  implicit none

  private

  public :: error_message, warning_message
  public :: WARN_ALWAYS, WARN_SOMETIMES, WARN_RARELY

  integer, parameter :: WARN_ALWAYS    = 0
  integer, parameter :: WARN_SOMETIMES = 1
  integer, parameter :: WARN_RARELY    = 2

contains

!================================ error_message ==============================80
!>
!! Description: Writes an error message to the screen.
!!
!! Inputs:      routine_name: Routine in which error is occuring
!!              message:      Error message to print to the screen
!!
!! Outputs:     err:          Flag indicating an error
!<
!=============================================================================80
  function error_message( routine_name, message ) result( err )

    use ISO_FORTRAN_ENV, only : error_unit

    implicit none

    character(*), intent(in) :: routine_name
    character(*), intent(in) :: message
    logical                  :: err

    continue

    err = .true.

    write(error_unit,*)
    write(error_unit,*) ' ERROR: In ' // trim(routine_name)
    write(error_unit,*) '   ', trim(message)
    write(error_unit,*) ' Stopping ...'
    call abort
    stop

  end function error_message

!=============================== warning_message =============================80
!>
!! Description: Writes a warning message to the screen.
!!
!! Inputs:      warn_level:   Important level for warning output
!!              routine_name: Routine in which warning is occuring
!!              message:      Warning message to print to the screen
!!
!! Outputs:     warn:         Flag indicating an warning
!<
!=============================================================================80
  function warning_message( warn_level, routine_name, message ) result( warn )

    use ISO_FORTRAN_ENV, only : error_unit
    use project_inputs,  only : verbose_level

    implicit none

    integer,      intent(in) :: warn_level
    character(*), intent(in) :: routine_name
    character(*), intent(in) :: message
    logical                  :: warn

    continue

    ! Setup
    warn = .true.

    ! Print Warning Message
    if ( warn_level <= verbose_level ) then
      write(error_unit,*)
      write(error_unit,*) ' WARNING: In ' // trim(routine_name)
      write(error_unit,*) '   ', trim(message)
    end if

  end function warning_message

end module message

module timer_derived_type

  use set_precision, only : dp
  use set_constants, only : zero

  implicit none

  private

  public :: basic_timer_t

  type :: basic_timer_t
    private
    real(dp)         :: time_start   = zero
    real(dp), public :: time_elapsed = zero
  contains
    private
    procedure, public, pass :: tic => timer_tick
    procedure, public, pass :: toc => timer_tock
  end type basic_timer_t

contains


function get_time()
  integer(kind=8) :: ticks, ticks_per_sec, max_ticks
  real(dp) :: get_time

  call system_clock( count      = ticks,                                     &
                     count_rate = ticks_per_sec,                             &
                     count_max  = max_ticks )

  if ( ticks_per_sec == 0 ) then
    get_time = zero
  else
    get_time = real(ticks,dp) / real(ticks_per_sec,dp)
  end if
end function get_time

  subroutine timer_tick( this )
    class(basic_timer_t), intent(inout) :: this
    this%time_elapsed = zero
    this%time_start   = get_time()
  end subroutine timer_tick

  function timer_tock( this )
    class(basic_timer_t), intent(in) :: this
    real(dp)                         :: timer_tock
    timer_tock = get_time() - this%time_start
  end function timer_tock

end module timer_derived_type

module index_conversion
  implicit none
  private
  public :: global2local, global2local_bnd, global2local_ghost
  public :: local2global, local2global_bnd, local2global_ghost
  public :: in_bound, cell_face_nbors
  public :: get_face_idx_from_id

  interface cell_face_nbors
    module procedure cell_face_nbors_lin
    module procedure cell_face_nbors_sub
  end interface cell_face_nbors
contains
  
  pure function global2local(iG,nSub) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: nSub
    integer, dimension(size(nSub)) :: iSub
    integer :: i, nDims, p, iGtmp, iTmp
    nDims = size(nSub)
    if (nDims==1) then
      iSub(1) = iG
      return
    end if
    p = product(nSub)
    iGtmp = iG
    do i = nDims,1,-1
      p = p/nSub(i)
      iTmp = mod(iGtmp-1,p) + 1
      iSub(i) = (iGtmp-iTmp)/p + 1
      iGtmp = iTmp
    end do
  end function global2local

  pure function local2global(iSub,nSub) result(iG)
    integer, dimension(:), intent(in) :: iSub, nSub
    integer :: iG
    integer :: nDims, p, i
    nDims = size(iSub)
    p = 1
    iG = 1
    do i = 1,nDims
        iG = iG + ( iSub(i) - 1 )*p
        p = p*nSub(i)
    end do
  end function local2global

  pure function global2local_ghost(iG,nSub,nGhost) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: nSub, nGhost
    integer, dimension(size(nSub)) :: iSub, nSub2
    nSub2 = nSub + 2*nGhost
    iSub = global2local(iG,nSub2)
    iSub = iSub - nGhost
  end function global2local_ghost

  pure function local2global_ghost(iSub,nSub,nGhost) result(iG)
    integer, dimension(:), intent(in) :: iSub, nSub, nGhost
    integer, dimension(size(nSub)) :: iSub2, nSub2
    integer :: iG
    iSub2 = iSub + nGhost
    nSub2 = nSub + 2*nGhost
    iG = local2global(iSub2,nSub2)
  end function local2global_ghost

  pure function global2local_bnd(iG,lo,hi) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: lo, hi
    integer, dimension(size(lo)) :: iSub, nSub
    nSub = hi - lo + 1
    iSub = global2local(iG,nSub)
    iSub = iSub + lo - 1
  end function global2local_bnd

  pure function local2global_bnd(iSub,lo,hi) result(iG)
    integer, dimension(:), intent(in) :: iSub, lo, hi
    integer, dimension(size(iSub)) :: idx, nSub
    integer :: iG
    idx  = iSub - lo + 1
    nSub = hi - lo + 1
    iG   = local2global(idx,nSub)
  end function local2global_bnd

  pure function in_bound( dim, idx, bnd_min, bnd_max )
    integer,                 intent(in) :: dim
    integer, dimension(dim), intent(in) :: idx, bnd_min, bnd_max
    logical                             :: in_bound
    in_bound =     all(idx>=bnd_min).and.all(idx<=bnd_max)                       &
              .or. all(idx<=bnd_min).and.all(idx>=bnd_max)
  end function in_bound

  pure subroutine cell_face_nbors_sub( dim, idx, bnd_min, bnd_max, nbor_cell_idx, nbor_face_id, n_int )
    integer,                       intent(in) :: dim
    integer, dimension(dim),       intent(in) :: idx, bnd_min, bnd_max
    integer, dimension(dim,2*dim), intent(out) :: nbor_cell_idx
    integer, dimension(2*dim),     intent(out) :: nbor_face_id
    integer,                       intent(out) :: n_int
    integer, dimension(dim,2*dim) :: nbor_cell_idx_tmp
    integer, dimension(2*dim) :: nbor_face_id_tmp
    integer, dimension(dim) :: idx_tmp
    integer :: s, j, n_ext, cnt
    cnt   = 0
    n_int = 0
    n_ext = 0
    do j = 1,dim
      do s = -1,1,2
        cnt = cnt + 1
        idx_tmp = idx
        idx_tmp(j) = idx_tmp(j) + s
        if ( in_bound(dim,idx_tmp,bnd_min,bnd_max) ) then
            n_int = n_int + 1
            nbor_cell_idx(:,n_int) = idx_tmp
            nbor_face_id(n_int) = cnt
        else
          n_ext = n_ext + 1
          nbor_cell_idx_tmp(:,n_ext) = idx_tmp
          nbor_face_id_tmp(n_ext) = cnt
        end if
      end do
    end do
    do j = 1,n_ext
      nbor_cell_idx(:,n_int+j) = nbor_cell_idx_tmp(:,j)
      nbor_face_id(n_int+j) = nbor_face_id_tmp(j)
    end do
  end subroutine cell_face_nbors_sub

  pure subroutine cell_face_nbors_lin( dim, lin_idx, bnd_min, bnd_max, &
                                       nbor_cell_idx, nbor_face_id, n_int )
    integer,                       intent(in) :: dim, lin_idx
    integer, dimension(dim),       intent(in) :: bnd_min, bnd_max
    integer, dimension(2*dim), intent(out) :: nbor_cell_idx
    integer, dimension(2*dim), intent(out) :: nbor_face_id
    integer,                       intent(out) :: n_int
    integer, dimension(dim,2*dim) :: nbor_idx
    integer, dimension(dim) :: idx
    integer :: s, j, n_ext, cnt
    idx = global2local_bnd(lin_idx,bnd_min,bnd_max)
    call cell_face_nbors_sub( dim, idx, bnd_min, bnd_max, nbor_idx, nbor_face_id, n_int )
    do j = 1,2*dim
      nbor_cell_idx(j) = local2global_bnd(nbor_idx(:,j),bnd_min,bnd_max)
    end do
  end subroutine cell_face_nbors_lin

  

  pure elemental subroutine get_face_info_from_id(face_id,dir,offset)
    integer, intent(in)  :: face_id
    integer, intent(out) :: dir, offset
    dir    = (face_id-1)/2 + 1
    offset = mod(face_id+1,2)
  end subroutine get_face_info_from_id

  pure subroutine get_face_idx_from_id(idx,face_id,dir,face_idx)
    integer, dimension(:),         intent(in) :: idx
    integer,                       intent(in)  :: face_id
    integer,                       intent(out) :: dir
    integer, dimension(size(idx)), intent(out) :: face_idx

    integer, dimension(size(idx)) :: face_offset
    integer :: offset
    call get_face_info_from_id(face_id,dir,offset)
    face_offset = 0
    face_offset(dir) = offset
    face_idx = idx + face_offset
  end subroutine get_face_idx_from_id

end module index_conversion

module combinatorics
  implicit none
  private
  public :: nchoosek
contains

  pure function nchoosek( n, k ) result( c )
    integer, intent(in) :: n, k
    integer             :: c
    integer :: i
    c = 0
    if (k>n) return

    c = 1
    do i = 1, min(n-k,k)
      c = c * ( n - (i-1) )
      c = c / i
    end do
  end function nchoosek

end module combinatorics

module math
  use set_precision, only : dp
  implicit none
  private
  public :: cross_product, vector_norm
contains

  pure function cross_product( vec1, vec2 )
    real(dp), dimension(3), intent(in) :: vec1, vec2
    real(dp), dimension(3)             :: cross_product
    cross_product(1) =  ( vec1(2)*vec2(3) - vec1(3)*vec2(2) )
    cross_product(2) = -( vec1(1)*vec2(3) - vec1(3)*vec2(1) )
    cross_product(3) =  ( vec1(1)*vec2(2) - vec1(2)*vec2(1) )
  end function cross_product

  pure function vector_norm( vector )
    use set_precision, only : dp
    use set_constants, only : zero
    real(dp), dimension(:), intent(in) :: vector
    real(dp)                           :: vector_norm
    integer :: i
    vector_norm = zero
    do i = 1, size(vector)
      vector_norm = vector_norm + vector(i)**2
    end do
    vector_norm = sqrt( vector_norm )
  end function vector_norm

end module math


module vector_derived_type
  use set_precision, only : dp
  use set_constants, only : zero
  implicit none
  private
  public :: face_vec
  public :: face_vec_ptr_3D
  type face_vec
    integer :: n
    real(dp), allocatable, dimension(:,:) :: v
  contains
    private
    procedure, public, pass :: create  => allocate_face_vec
    procedure, public, pass :: destroy => deallocate_face_vec
  end type face_vec

  type face_vec_ptr_3D
    type(face_vec), dimension(:,:,:), pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_face_vec_ptr_3D
  end type face_vec_ptr_3D

contains

  subroutine allocate_face_vec( this, n )
    class(face_vec), intent(inout) :: this
    integer,       intent(in)      :: n
    continue
    this%n = n
    allocate( this%v(3,n) )
    this%v = zero
  end subroutine allocate_face_vec

  pure elemental subroutine deallocate_face_vec( this )
    class(face_vec), intent(inout) :: this
    continue
    this%n = 0
    if( allocated( this%v  ) ) deallocate( this%v )
  end subroutine deallocate_face_vec

  pure elemental subroutine destroy_face_vec_ptr_3D( this )
    class(face_vec_ptr_3D), intent(inout) :: this
    this%p => null()
  end subroutine destroy_face_vec_ptr_3D
end module vector_derived_type

module pointers
  use set_precision, only : dp
  implicit none
  private
  public :: array_ptr_3D_real, array_ptr_4D_real

  type array_ptr_3D_real
    real(dp), dimension(:,:,:),     pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_real_3D
  end type array_ptr_3D_real

  type array_ptr_4D_real
    real(dp), dimension(:,:,:,:),   pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_real_4D
  end type array_ptr_4D_real

contains

  pure elemental subroutine destroy_real_3D( this )
    class(array_ptr_3D_real), intent(inout) :: this
    this%p => null()
  end subroutine destroy_real_3D

  pure elemental subroutine destroy_real_4D( this )
    class(array_ptr_4D_real), intent(inout) :: this
    this%p => null()
  end subroutine destroy_real_4D
end module pointers

module linspace_helper
  use set_precision, only : dp
  private
  public :: unit_cartesian_mesh_cat
  public :: linspace, meshgrid2, meshgrid3
contains
  pure function unit_cartesian_mesh_cat(nx,ny,nz) result(xyz)
    integer, intent(in) :: nx, ny, nz
    real(dp), dimension(3,nx,ny,nz) :: xyz
    real(dp), dimension(nx,ny,nz) :: tmp_x, tmp_y, tmp_z
    integer :: i, j, k

    call unit_cartesian_mesh(nx,ny,nz,tmp_x,tmp_y,tmp_z)

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          xyz(1,i,j,k) = tmp_x(i,j,k)
          xyz(2,i,j,k) = tmp_y(i,j,k)
          xyz(3,i,j,k) = tmp_z(i,j,k)
        end do
      end do
    end do
  end function unit_cartesian_mesh_cat

  pure subroutine unit_cartesian_mesh(nx,ny,nz,x,y,z)
    use set_constants, only : zero, one
    integer, intent(in) :: nx, ny, nz
    real(dp), dimension(nx,ny,nz), intent(out) :: x, y, z

    call meshgrid3( linspace(nx,zero,one), &
                    linspace(ny,zero,one), &
                    linspace(nz,zero,one), x,y,z)
  end subroutine unit_cartesian_mesh

  pure function linspace(N,x1,x2) result(array)
    integer,  intent(in)   :: N
    real(dp), intent(in)   :: x1, x2
    real(dp), dimension(N) :: array
    real(dp) :: range_den
    integer :: i
    if (N==0) return
    if (N==1) then
      array(1) = x1
      return
    end if
    range_den = (x2-x1)/real(N-1,dp)
    do i = 1,N
      array(i) = x1 + range_den*real(i-1,dp)
    end do
  end function linspace

  pure subroutine meshgrid2(x1,x2,x1_array,x2_array)
    real(dp), dimension(:),   intent(in)  :: x1, x2
    real(dp), dimension(:,:), intent(out) :: x1_array, x2_array
    integer :: N1, N2
    N1 = size(x1)
    N2 = size(x2)
    x1_array = spread(x1,2,N2)
    x2_array = spread(x2,1,N1)
  end subroutine meshgrid2

  pure subroutine meshgrid3(x1,x2,x3,x1_array,x2_array,x3_array)
    real(dp), dimension(:),     intent(in)  :: x1, x2, x3
    real(dp), dimension(:,:,:), intent(out) :: x1_array, x2_array, x3_array
    real(dp), dimension(size(x1),size(x2)) :: x1_tmp
    real(dp), dimension(size(x2),size(x3)) :: x2_tmp
    real(dp), dimension(size(x2),size(x3),size(x1)) :: x2_tmp2
    real(dp), dimension(size(x3),size(x1)) :: x3_tmp
    real(dp), dimension(size(x3),size(x1),size(x2)) :: x3_tmp2
    integer, parameter, dimension(3) :: o2 = [2,3,1], o3 = [3,1,2]
    integer :: N1, N2, N3
    N1 = size(x1)
    N2 = size(x2)
    N3 = size(x3)

    x1_tmp   = spread(x1,2,N2)
    x2_tmp   = spread(x2,2,N3)
    x3_tmp   = spread(x3,2,N1)
    x1_array = spread(x1_tmp,3,N3)
    x2_tmp2  = spread(x2_tmp,3,N1)
    x3_tmp2  = spread(x3_tmp,3,N2)
    x2_array = reshape(x2_tmp2,shape(x2_array),order=o2)
    x3_array = reshape(x3_tmp2,shape(x3_array),order=o3)
  end subroutine meshgrid3

end module linspace_helper

module lagrange_interpolation
  use set_precision, only : dp
  use set_constants, only : zero, one, two, half
  implicit none
  private
  ! public :: xb, wb, Dmat
  public :: generate_1D_barycentric_info, destroy_1D_barycentric_info
  public :: lagbary, lagbary_wderiv
  public :: lagbary_2D, lagbary_2D_wgrad
  public :: lagbary_3D, lagbary_3D_wgrad
  public :: calc_grid_metrics, jacobian_determinant
  public :: normal_vectors
  interface jacobian_determinant
    module procedure jacobian_determinant_2D
    module procedure jacobian_determinant_3D
  end interface jacobian_determinant

  integer :: Nmax = 10
  real(dp), dimension(:,:,:,:), allocatable :: Dmat
  real(dp), dimension(:,:),   allocatable :: xb, wb
contains

  elemental logical function almost_equal(a,b)
    real(dp), intent(in) :: a, b
    logical :: test1, test2, test3
    test1 = ( (a==zero) .or. (b==zero) )
    test2 = ( abs(a-b) <= two*epsilon(one) )
    test3 = ( ( abs(a-b) <= epsilon(abs(a)) ) .and. &
            ( abs(a-b) <= epsilon(abs(b)) ) )
    almost_equal = ( ( test1 .and. test2 ) .or. ( (.not. test1) .and. test3 ) )
  end function almost_equal

  function barycentric_weights(x) result(w)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x))       :: w
    integer :: j, k, N
    N = size(x)
    w = one
    do j = 2,N
      do k = 1,j-1
        w(k) = w(k) * ( x(k) - x(j) )
        w(j) = w(j) * ( x(j) - x(k) )
      end do
    end do
    w = one/w
  end function barycentric_weights

  pure function polynomial_derivative_matrix(x,w) result(D)
    real(dp), dimension(:), intent(in) :: x, w
    real(dp), dimension(size(x),size(x)) :: D
    integer :: i, j, N
    D = zero
    N = size(x)
    do i = 1,N
      do j = 1,N
        if (j/=i) then
          D(i,j) = w(j)/w(i) * one / ( x(i) - x(j) )
          D(i,i) = D(i,i) - D(i,j)
        end if
      end do
    end do
  end function polynomial_derivative_matrix

  pure function mth_order_polynomial_derivative_matrix(x,w,M) result(D)
    real(dp), dimension(:), intent(in) :: x, w
    integer,                intent(in) :: M
    real(dp), dimension(size(x),size(x),M) :: D
    integer :: i, j, k, N
    D = zero
    N = size(x)
    D(:,:,1) = polynomial_derivative_matrix(x,w)
    do k = 2,M
      do i = 1,N
        D(i,i,k) = zero
        do j = 1,N
          if (j/=i) then
            D(i,j,k) = ( real(k,dp) / (x(i) - x(j)) )                          &
                     * ( w(j)/w(i)*D(i,i,k-1) - D(i,j,k-1) )
            D(i,i,k) = D(i,i,k) - D(i,j,k)
          end if
        end do
      end do
    end do
  end function mth_order_polynomial_derivative_matrix

  subroutine generate_1D_barycentric_info(N)
    use linspace_helper, only : linspace
    integer, intent(in), optional :: N
    integer :: j
    real(dp) :: x1 = -one, x2 = one
    if ( present(N) ) then
      if ( N > Nmax ) then
        Nmax = N
        continue
      end if
    end if
    if ( allocated(Dmat) ) deallocate(Dmat)
    if ( allocated(xb) )   deallocate(xb)
    if ( allocated(wb) )   deallocate(wb)
    allocate( Dmat(Nmax,Nmax,Nmax,2), xb(Nmax,Nmax), wb(Nmax,Nmax) )
    Dmat = zero; xb = zero; wb = zero
    wb(1,1) = one
    do j = 2,Nmax
      ! call linspace(x1,x2,xb(1:j,j))
      xb(1:j,j) = linspace(j,x1,x2)
      wb(1:j,j) = barycentric_weights( xb(1:j,j) )
      Dmat(1:j,1:j,j,:) = mth_order_polynomial_derivative_matrix( xb(1:j,j), wb(1:j,j), 2 )
    end do
  end subroutine generate_1D_barycentric_info

  subroutine destroy_1D_barycentric_info
    if ( allocated(Dmat) ) deallocate(Dmat)
    if ( allocated(xb) ) deallocate(xb)
    if ( allocated(wb) ) deallocate(wb)
  end subroutine destroy_1D_barycentric_info

  pure subroutine lagbary(x,dir,fval,Npts,val)
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val
    real(dp) :: A, F
    real(dp) :: x1, t1
    integer :: j, N
    A = zero
    F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val = fval(j)
        return
      end if
      t1 = wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
    end do
    val = A/F
  end subroutine lagbary

  pure subroutine lagbary_wderiv(x,dir,fval,Npts,val,dval)
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val, dval
    real(dp) :: A, B, C, F
    real(dp) :: x1, t1, t2, FF, AC
    integer :: j, N
    A = zero
    B = zero
    C = zero
    F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val = fval(j)
        dval = dot_product( Dmat(j,1:N,N,1), fval )
        return
      end if
      t1 = wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
      t2 = t1/x1
      B = B + t2 * fval(j)
      C = C + t2
    end do
    val = A/F
    FF = F*F
    AC = A*C
    dval = (B * F - AC)/FF
  end subroutine lagbary_wderiv

  pure subroutine lagbary_wderiv2(x,dir,fval,Npts,val,dval,d2val)
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val, dval, d2val
    real(dp) :: A, B, C, D, E, F
    real(dp) :: x1, t1, t2, t3, FF, AC
    integer :: j, N
    A = zero
    B = zero
    C = zero
    D = zero
    E = zero
    F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val   = fval(j)
        dval  = dot_product( Dmat(j,1:N,N,1), fval )
        d2val = dot_product( Dmat(j,1:N,N,2), fval )
        return
      end if
      t1 = wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
      t2 = t1/x1
      B = B + t2 * fval(j)
      C = C + t2
      t3 = t2/x1
      D = D + t3 * fval(j)
      E = E + t3
    end do
    val = A/F
    FF = F*F
    AC = A*C
    dval = (B * F - AC)/FF
    d2val = ( two * D ) / F - ( two * E * A ) / FF - ( two * B * C ) / FF      &
          + ( two * C * AC ) / ( FF * F )
  end subroutine lagbary_wderiv2

  pure subroutine lagbary_2D(x,fval,Npts,val)
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(size(fval,2)) :: tmp
    integer :: j
    do j = 1,Npts(2)
      call lagbary( x(1), 1, fval(:,j), Npts, tmp(j) )
    end do
    call lagbary( x(2), 2, tmp, Npts, val )
  end subroutine lagbary_2D

  pure subroutine lagbary_2D_wgrad(x,fval,Npts,val,grad)
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(2),   intent(out) :: grad
    real(dp), dimension(size(fval,2)) :: tmp, gtmp
    integer :: j
    do j = 1,Npts(2)
      call lagbary_wderiv( x(1), 1, fval(:,j), Npts, tmp(j), gtmp(j) )
    end do
    call lagbary_wderiv( x(2), 2,  tmp, Npts, val, grad(2) )
    call lagbary(        x(2), 2, gtmp, Npts,      grad(1) )
  end subroutine lagbary_2D_wgrad

  pure subroutine lagbary_2D_whess(x,fval,Npts,val,grad,hess)
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(2),   intent(out) :: grad
    real(dp), dimension(3),   intent(out) :: hess
    real(dp), dimension(size(fval,2)) :: tmp, gtmp, htmp
    integer :: j
    do j = 1,Npts(2)
      call lagbary_wderiv2( x(1), 1, fval(:,j), Npts, tmp(j), gtmp(j), htmp(j) )
    end do
    call lagbary_wderiv2( x(2), 2,  tmp, Npts, val, grad(2), hess(3) )
    call lagbary_wderiv(  x(2), 2, gtmp, Npts,      grad(1), hess(2) )
    call lagbary(         x(2), 2, htmp, Npts,               hess(1) )
  end subroutine lagbary_2D_whess

  pure subroutine lagbary_3D(x,fval,Npts,val)
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp
    real(dp), dimension(size(fval,3)) :: tmp2
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call lagbary( x(1), 1, fval(:,j,k), Npts, tmp(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call lagbary( x(2), 2, tmp(:,k), Npts, tmp2(k) )
    end do
    call lagbary( x(3), 3, tmp2, Npts, val )
  end subroutine lagbary_3D

  pure subroutine lagbary_3D_wgrad(x,fval,Npts,val,grad)
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(3),     intent(out) :: grad
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp, gtmp0
    real(dp), dimension(size(fval,3)) :: tmp2, gtmp1, gtmp2
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call lagbary_wderiv( x(1), 1, fval(:,j,k), Npts, tmp(j,k), gtmp0(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call lagbary_wderiv( x(2), 2,   tmp(:,k), Npts, tmp2(k), gtmp2(k) )
      call lagbary(        x(2), 2, gtmp0(:,k), Npts, gtmp1(k) )
    end do
    call lagbary_wderiv( x(3), 3,  tmp2, Npts, val, grad(3) )
    call lagbary(        x(3), 3, gtmp2, Npts,      grad(2) )
    call lagbary(        x(3), 3, gtmp1, Npts,      grad(1) )
  end subroutine lagbary_3D_wgrad

  pure subroutine lagbary_3D_whess(x,fval,Npts,val,grad,hess)
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(3),     intent(out) :: grad
    real(dp), dimension(6),     intent(out) :: hess
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp, gtmp, htmp
    real(dp), dimension(size(fval,3)) :: tmp1, gtmp1, gtmp2, htmp1, htmp2, htmp3
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call lagbary_wderiv2( x(1), 1, fval(:,j,k), Npts, tmp(j,k), gtmp(j,k), htmp(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call lagbary_wderiv2( x(2), 2,  tmp(:,k), Npts, tmp1(k), gtmp2(k), htmp3(k) )
      call lagbary_wderiv(  x(2), 2, gtmp(:,k), Npts,          gtmp1(k), htmp2(k) )
      call lagbary(         x(2), 2, htmp(:,k), Npts,                    htmp1(k) )
    end do
    call lagbary_wderiv2( x(3), 3,  tmp1, Npts, val, grad(3), hess(6) )
    call lagbary_wderiv(  x(3), 3, gtmp2, Npts,      grad(2), hess(5) )
    call lagbary(         x(3), 3, htmp3, Npts,               hess(4) )
    call lagbary_wderiv(  x(3), 3, gtmp1, Npts,      grad(1), hess(3) )
    call lagbary(         x(3), 3, htmp2, Npts,               hess(2) )
    call lagbary(         x(3), 3, htmp1, Npts,               hess(1) )
  end subroutine lagbary_3D_whess

  pure function covariant_base_vectors_3D(point,X1,X2,X3) result(a)
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: a
    real(dp) :: junk
    integer, dimension(3) :: Npts
    Npts = shape(X1)
    call lagbary_3D_wgrad(point,X1,Npts,junk,a(1,:))
    call lagbary_3D_wgrad(point,X2,Npts,junk,a(2,:))
    call lagbary_3D_wgrad(point,X3,Npts,junk,a(3,:))
  end function covariant_base_vectors_3D

  function contravariant_base_vectors_3D(point,X1,X2,X3) result(a)
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: a
    real(dp) :: J
    a = calc_grid_metrics(point,X1,X2,X3)
    J = jacobian_determinant_3D(point,X1,X2,X3)
    a = a / J
  end function contravariant_base_vectors_3D

  function calc_grid_metrics(point,X1,X2,X3) result(Ja)
    use pointers, only : array_ptr_3D_real
    real(dp), dimension(3), intent(in) :: point
    real(dp), dimension(:,:,:), intent(in), target :: X1, X2, X3
    real(dp), dimension(3,3) :: Ja
    real(dp), dimension(size(X1,1),size(X1,2),size(X1,3)) :: X_l, X_m, tmp
    type(array_ptr_3D_real), dimension(3) :: X
    integer, dimension(3) :: Npts
    integer :: i
    real(dp) :: junk
    real(dp), dimension(3) :: dX_l, dX_m, dd1, dd2, dd3
    Ja = zero
    X(1)%p => X1
    X(2)%p => X2
    X(3)%p => X3
    Npts = shape(X1)
    do i = 1,3
      X_l = X(mod(4+i,3)+1)%p
      X_m = X(mod(3+i,3)+1)%p
      call lagbary_3D_wgrad( point, X_l, Npts, junk, dX_l )
      call lagbary_3D_wgrad( point, X_m, Npts, junk, dX_m )
      tmp = X_l*dX_m(1) - X_m*dX_l(1)
      call lagbary_3D_wgrad( point, tmp, Npts, junk, dd1 )
      tmp = X_l*dX_m(2) - X_m*dX_l(2)
      call lagbary_3D_wgrad( point, tmp, Npts, junk, dd2 )
      tmp = X_l*dX_m(3) - X_m*dX_l(3)
      call lagbary_3D_wgrad( point, tmp, Npts, junk, dd3 )
      Ja(i,1) = -half*( dd3(2) - dd2(3) );
      Ja(i,2) = -half*( dd1(3) - dd3(1) );
      Ja(i,3) = -half*( dd2(1) - dd1(2) );
    end do
    nullify( X(1)%p, X(2)%p, X(3)%p ) ! may not be necessary
  end function calc_grid_metrics


  pure function normal_vectors(point,X1,X2,X3) result(Nvec)
    use math, only : cross_product, vector_norm
    implicit none
    real(dp), dimension(3), intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: Nvec
    integer, dimension(3) :: Npts
    real(dp), dimension(3,3) :: A
    real(dp) :: junk
    Npts = shape(X1)
    call lagbary_3D_wgrad(point,X1,Npts,junk,A(:,1))
    call lagbary_3D_wgrad(point,X2,Npts,junk,A(:,2))
    call lagbary_3D_wgrad(point,X3,Npts,junk,A(:,3))

    Nvec(:,1) = cross_product(A(:,2),A(:,3))
    Nvec(:,2) = cross_product(A(:,1),A(:,3))
    Nvec(:,3) = cross_product(A(:,1),A(:,2))

    Nvec(:,1) = Nvec(:,1)/vector_norm(Nvec(:,1))
    Nvec(:,2) = Nvec(:,2)/vector_norm(Nvec(:,2))
    Nvec(:,3) = Nvec(:,3)/vector_norm(Nvec(:,3))
  end function normal_vectors

  pure function jacobian_determinant_2D(point,X1,X2,X3) result(Jac)
    implicit none
    real(dp), dimension(2), intent(in) :: point
    real(dp), dimension(:,:), intent(in) :: X1, X2, X3
    real(dp) :: Jac, junk
    integer, dimension(2) :: Npts
    real(dp), dimension(2,3) :: A
    Npts = shape(X1)
    call lagbary_2D_wgrad(point,X1,Npts,junk,A(:,1))
    call lagbary_2D_wgrad(point,X2,Npts,junk,A(:,2))
    call lagbary_2D_wgrad(point,X3,Npts,junk,A(:,3))
    Jac = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  end function jacobian_determinant_2D

  pure function jacobian_determinant_3D(point,X1,X2,X3) result(Jac)
    implicit none
    real(dp), dimension(3), intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp) :: Jac, junk
    integer, dimension(3) :: Npts
    real(dp), dimension(3,3) :: A
    Npts = shape(X1)
    call lagbary_3D_wgrad(point,X1,Npts,junk,A(:,1))
    call lagbary_3D_wgrad(point,X2,Npts,junk,A(:,2))
    call lagbary_3D_wgrad(point,X3,Npts,junk,A(:,3))
    Jac = A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(3,2)*A(1,3) + A(3,1)*A(1,2)*A(2,3) &
        - A(3,1)*A(2,2)*A(1,3) - A(2,1)*A(1,2)*A(3,3) - A(1,1)*A(3,2)*A(2,3)
  end function jacobian_determinant_3D

end module lagrange_interpolation

program main
  use lagrange_interpolation, only : generate_1D_barycentric_info, destroy_1D_barycentric_info
  call generate_1D_barycentric_info()
  call destroy_1D_barycentric_info()

end program main