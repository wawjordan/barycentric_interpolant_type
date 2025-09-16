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
  public :: cross_product
  public :: det_3x3
  public :: LegendrePolynomialAndDerivative
  public :: LegendreGaussNodesAndWeights
contains

  pure function cross_product( vec1, vec2 )
    real(dp), dimension(3), intent(in) :: vec1, vec2
    real(dp), dimension(3)             :: cross_product
    cross_product(1) =  ( vec1(2)*vec2(3) - vec1(3)*vec2(2) )
    cross_product(2) = -( vec1(1)*vec2(3) - vec1(3)*vec2(1) )
    cross_product(3) =  ( vec1(1)*vec2(2) - vec1(2)*vec2(1) )
  end function cross_product

  pure function det_3x3( mat )
    real(dp), dimension(3,3), intent(in) :: mat
    real(dp)                             :: det_3x3
    continue
    det_3x3 = mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
            - mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1)) &
            + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))
  end function det_3x3

  elemental subroutine LegendrePolynomialAndDerivative(N,x,LN,dLN)
    use set_constants, only : zero, one, two
    integer, intent(in) :: N
    real(dp), intent(in) :: x
    real(dp), intent(out) :: LN, dLN
    real(dp) :: LNm2, LNm1, dLNm2, dLNm1
    integer :: j
    if (N == 0) then
      LN = one
      dLN = zero
    elseif (N == 1) then
      LN = x
      dLN = one
    else
      LNm2 = one
      LNm1 = x
      dLNm2 = zero
      dLNm1 = one
      do j = 2,N
        LN = real(2*j-1,dp)/real(j,dp) * x * LNm1 &
          - real(j-1,dp)/real(j,dp) * LNm2
        dLN = dLNm2 + real(2*j-1,dp) * LNm1
        LNm2 = LNm1
        LNm1 = LN
        dLNm2 = dLNm1
        dLNm1 = dLN
      end do
    end if
  end subroutine LegendrePolynomialAndDerivative

  pure subroutine LegendreGaussNodesAndWeights(N,x,w)
    use set_constants, only : zero, one, two, four, third, pi
    integer,                  intent(in)  :: N
    real(dp), dimension(N+1), intent(out) :: x, w
    integer :: j, k
    real(dp) :: eps4, delta, LNp1, dLNp1
    integer, parameter :: quad_n_iter = 10
    eps4 = four*epsilon(one)
    x = zero
    w = zero

    if (N == 0) then
      x(1) = zero
      w(1) = two
    elseif (N == 1) then
      x(1) = -sqrt(third)
      w(1) = one
      x(2) = -x(1)
      w(2) = w(1)
    else
      do j = 0,(N+1)/2 - 1
        x(j+1) = -cos( ( real(2*j+1,dp)/real(2*N+2,dp) )*pi )
        do k = 1,quad_n_iter
          call LegendrePolynomialAndDerivative(N+1,x(j+1),LNp1,dLNp1)
          delta = -LNp1/dLNp1
          x(j+1) = x(j+1) + delta
          if ( abs(delta) <= eps4*abs(x(j+1)) ) then
            exit
          end if
        end do
        call LegendrePolynomialAndDerivative(N+1,x(j+1),LNp1,dLNp1)
        x(N+1-j) = -x(j+1)
        w(j+1) = two/( (one-x(j+1)**2)*dLNp1**2)
        w(N+1-j) = w(j+1)
      end do
      if (mod(N,2) == 0) then
        call LegendrePolynomialAndDerivative(N+1,zero,LNp1,dLNp1)
        x(N/2+1) = zero
        w(N/2+1) = two/dLNp1**2
      end if
    end if
  end subroutine LegendreGaussNodesAndWeights

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


module barycentric_interpolant_derived_type
  use set_precision, only : dp
  use set_constants, only : zero, one, two, half
  implicit none
  private
  public :: interpolant_t

  type :: interpolant_t
    integer :: Nmax
    real(dp), dimension(:,:),     allocatable :: xb, wb
    real(dp), dimension(:,:,:,:), allocatable :: Dmat
  contains
    private
    procedure, public, nopass :: constructor
    procedure, public, pass   :: destroy       => destroy_interpolant
    procedure,         pass   :: lagbary, lagbary_wderiv, lagbary_wderiv2
    procedure,         pass   :: lagbary_2D, lagbary_2D_wgrad, lagbary_2D_whess
    procedure,         pass   :: lagbary_3D, lagbary_3D_wgrad, lagbary_3D_whess
    procedure, public, pass   :: calc_grid_metrics, calc_grid_metrics_alt
    procedure, public, pass   :: normal_vectors
    procedure, public, pass   :: map_point => map_point_3D_curve, map_point_3D_surface, map_point_3D_volume
  end type interpolant_t

  interface interpolant_t
    procedure constructor
  end interface interpolant_t

  ! interface map_point
  !   module procedure map_point_3D_curve
  !   module procedure map_point_3D_surface
  !   module procedure map_point_3D_volume
  ! end interface map_point

contains

  pure elemental subroutine destroy_interpolant(this)
    class(interpolant_t), intent(inout) :: this
    if ( allocated(this%Dmat) ) deallocate(this%Dmat)
    if ( allocated(this%xb) )   deallocate(this%xb)
    if ( allocated(this%wb) )   deallocate(this%wb)
    this%Nmax = 0
  end subroutine destroy_interpolant

  pure elemental function constructor(N) result(this)
    use linspace_helper, only : linspace
    use set_constants,   only : zero, one
    integer, optional, intent(in) :: N
    type(interpolant_t)           :: this
    integer :: j
    call this%destroy()
    if ( present(N) ) this%Nmax = max(N,2)
    allocate( this%Dmat(this%Nmax,this%Nmax,this%Nmax,2) )
    allocate(   this%xb(this%Nmax,this%Nmax), this%wb(this%Nmax,this%Nmax) )
    this%Dmat = zero; this%xb = zero; this%wb = zero
    this%wb(1,1) = one
    do j = 2,this%Nmax
      this%xb(1:j,j) = linspace(j,-one,one)
      this%wb(1:j,j) = barycentric_weights( this%xb(1:j,j) )
      this%Dmat(1:j,1:j,j,:) = mth_order_polynomial_derivative_matrix( this%xb(1:j,j), this%wb(1:j,j), 2 )
    end do
  end function constructor

  pure elemental logical function almost_equal(a,b)
    real(dp), intent(in) :: a, b
    logical :: test1, test2, test3
    test1 = ( (a==zero) .or. (b==zero) )
    test2 = ( abs(a-b) <= two*epsilon(one) )
    test3 = ( ( abs(a-b) <= epsilon(abs(a)) ) .and. &
              ( abs(a-b) <= epsilon(abs(b)) ) )
    almost_equal = ( ( test1 .and. test2 ) .or. ( (.not. test1) .and. test3 ) )
  end function almost_equal

  pure function barycentric_weights(x) result(w)
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
            D(i,j,k) = ( real(k,dp) / (x(i) - x(j)) )        &
                     * ( w(j)/w(i)*D(i,i,k-1) - D(i,j,k-1) )
            D(i,i,k) = D(i,i,k) - D(i,j,k)
          end if
        end do
      end do
    end do
  end function mth_order_polynomial_derivative_matrix

  pure subroutine lagbary(this,x,dir,fval,Npts,val)
    class(interpolant_t),    intent(in)  :: this
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val
    real(dp) :: A, F
    real(dp) :: x1, t1
    integer :: j, N
    A = zero; F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = this%xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val = fval(j)
        return
      end if
      t1 = this%wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
    end do
    val = A/F
  end subroutine lagbary

  pure subroutine lagbary_wderiv(this,x,dir,fval,Npts,val,dval)
    class(interpolant_t),    intent(in)  :: this
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val, dval
    real(dp) :: A, B, C, F
    real(dp) :: x1, t1, t2, FF, AC
    integer :: j, N
    A = zero; B = zero; C = zero; F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = this%xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val = fval(j)
        dval = dot_product( this%Dmat(j,1:N,N,1), fval )
        return
      end if
      t1 = this%wb(j,N)/x1
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

  pure subroutine lagbary_wderiv2(this,x,dir,fval,Npts,val,dval,d2val)
    class(interpolant_t),    intent(in)  :: this
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val, dval, d2val
    real(dp) :: A, B, C, D, E, F
    real(dp) :: x1, t1, t2, t3, FF, AC
    integer :: j, N
    A = zero; B = zero; C = zero; D = zero; E = zero; F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = this%xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val   = fval(j)
        dval  = dot_product( this%Dmat(j,1:N,N,1), fval )
        d2val = dot_product( this%Dmat(j,1:N,N,2), fval )
        return
      end if
      t1 = this%wb(j,N)/x1
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
    d2val = ( two * D      ) / F          &
          - ( two * E * A  ) / FF         &
          - ( two * B * C  ) / FF         &
          + ( two * C * AC ) / ( FF * F )
  end subroutine lagbary_wderiv2

  pure subroutine lagbary_2D(this,x,fval,Npts,val)
    class(interpolant_t),      intent(in)  :: this
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(size(fval,2)) :: tmp
    integer :: j
    do j = 1,Npts(2)
      call this%lagbary( x(1), 1, fval(:,j), Npts, tmp(j) )
    end do
    call this%lagbary( x(2), 2, tmp, Npts, val )
  end subroutine lagbary_2D

  pure subroutine lagbary_2D_wgrad(this,x,fval,Npts,val,grad)
    class(interpolant_t),      intent(in)  :: this
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(2),   intent(out) :: grad
    real(dp), dimension(size(fval,2)) :: tmp, gtmp
    integer :: j
    do j = 1,Npts(2)
      call this%lagbary_wderiv( x(1), 1, fval(:,j), Npts, tmp(j), gtmp(j) )
    end do
    call this%lagbary_wderiv( x(2), 2,  tmp, Npts, val, grad(2) )
    call this%lagbary(        x(2), 2, gtmp, Npts,      grad(1) )
  end subroutine lagbary_2D_wgrad

  pure subroutine lagbary_2D_whess(this,x,fval,Npts,val,grad,hess)
    class(interpolant_t),      intent(in)  :: this
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(2),   intent(out) :: grad
    real(dp), dimension(3),   intent(out) :: hess
    real(dp), dimension(size(fval,2)) :: tmp, gtmp, htmp
    integer :: j
    do j = 1,Npts(2)
      call this%lagbary_wderiv2( x(1), 1, fval(:,j), Npts, tmp(j), gtmp(j), htmp(j) )
    end do
    call this%lagbary_wderiv2( x(2), 2,  tmp, Npts, val, grad(2), hess(3) )
    call this%lagbary_wderiv(  x(2), 2, gtmp, Npts,      grad(1), hess(2) )
    call this%lagbary(         x(2), 2, htmp, Npts,               hess(1) )
  end subroutine lagbary_2D_whess

  pure subroutine lagbary_3D(this,x,fval,Npts,val)
    class(interpolant_t),        intent(in)  :: this
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp
    real(dp), dimension(size(fval,3)) :: tmp2
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call this%lagbary( x(1), 1, fval(:,j,k), Npts, tmp(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call this%lagbary( x(2), 2, tmp(:,k), Npts, tmp2(k) )
    end do
    call this%lagbary( x(3), 3, tmp2, Npts, val )
  end subroutine lagbary_3D

  pure subroutine lagbary_3D_wgrad(this,x,fval,Npts,val,grad)
    class(interpolant_t),        intent(in)  :: this
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
        call this%lagbary_wderiv( x(1), 1, fval(:,j,k), Npts, tmp(j,k), gtmp0(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call this%lagbary_wderiv( x(2), 2,   tmp(:,k), Npts, tmp2(k), gtmp2(k) )
      call this%lagbary(        x(2), 2, gtmp0(:,k), Npts, gtmp1(k) )
    end do
    call this%lagbary_wderiv( x(3), 3,  tmp2, Npts, val, grad(3) )
    call this%lagbary(        x(3), 3, gtmp2, Npts,      grad(2) )
    call this%lagbary(        x(3), 3, gtmp1, Npts,      grad(1) )
  end subroutine lagbary_3D_wgrad

  pure subroutine lagbary_3D_whess(this,x,fval,Npts,val,grad,hess)
    class(interpolant_t),        intent(in)  :: this
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
        call this%lagbary_wderiv2( x(1), 1, fval(:,j,k), Npts, tmp(j,k), gtmp(j,k), htmp(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call this%lagbary_wderiv2( x(2), 2,  tmp(:,k), Npts, tmp1(k), gtmp2(k), htmp3(k) )
      call this%lagbary_wderiv(  x(2), 2, gtmp(:,k), Npts,          gtmp1(k), htmp2(k) )
      call this%lagbary(         x(2), 2, htmp(:,k), Npts,                    htmp1(k) )
    end do
    call this%lagbary_wderiv2( x(3), 3,  tmp1, Npts, val, grad(3), hess(6) )
    call this%lagbary_wderiv(  x(3), 3, gtmp2, Npts,      grad(2), hess(5) )
    call this%lagbary(         x(3), 3, htmp3, Npts,               hess(4) )
    call this%lagbary_wderiv(  x(3), 3, gtmp1, Npts,      grad(1), hess(3) )
    call this%lagbary(         x(3), 3, htmp2, Npts,               hess(2) )
    call this%lagbary(         x(3), 3, htmp1, Npts,               hess(1) )
  end subroutine lagbary_3D_whess

  pure function calc_grid_metrics(this,point,X1,X2,X3) result(Ja)
    use set_constants, only : zero
    class(interpolant_t),       intent(in) :: this
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: Ja
    real(dp), dimension(size(X1,1),size(X1,2),size(X1,3),3) :: X
    real(dp), dimension(size(X1,1),size(X1,2),size(X1,3))   :: tmp
    real(dp), dimension(3) :: dX_l, dX_m, dd1, dd2, dd3
    real(dp) :: junk
    integer, dimension(3), parameter :: ijk = [1,2,3]
    integer, dimension(3), parameter :: kij = cshift(ijk,1)
    integer, dimension(3), parameter :: jki = cshift(kij,1)
    integer, dimension(3) :: Npts
    integer :: i
    Ja = zero
    Npts = shape(X1)
    call this%lagbary_3D_wgrad( point, X3, Npts, junk, dX_l )
    call this%lagbary_3D_wgrad( point, X2, Npts, junk, dX_m )
    tmp = X3*dX_m(1) - X2*dX_l(1)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd1 )
    tmp = X3*dX_m(2) - X2*dX_l(2)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd2 )
    tmp = X3*dX_m(3) - X2*dX_l(3)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd3 )
    Ja(1,1) = -half*( dd3(2) - dd2(3) );
    Ja(1,2) = -half*( dd1(3) - dd3(1) );
    Ja(1,3) = -half*( dd2(1) - dd1(2) );

    call this%lagbary_3D_wgrad( point, X1, Npts, junk, dX_l )
    call this%lagbary_3D_wgrad( point, X3, Npts, junk, dX_m )
    tmp = X1*dX_m(1) - X3*dX_l(1)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd1 )
    tmp = X1*dX_m(2) - X3*dX_l(2)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd2 )
    tmp = X1*dX_m(3) - X3*dX_l(3)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd3 )
    Ja(2,1) = -half*( dd3(2) - dd2(3) );
    Ja(2,2) = -half*( dd1(3) - dd3(1) );
    Ja(2,3) = -half*( dd2(1) - dd1(2) );

    call this%lagbary_3D_wgrad( point, X2, Npts, junk, dX_l )
    call this%lagbary_3D_wgrad( point, X1, Npts, junk, dX_m )
    tmp = X2*dX_m(1) - X1*dX_l(1)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd1 )
    tmp = X2*dX_m(2) - X1*dX_l(2)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd2 )
    tmp = X2*dX_m(3) - X1*dX_l(3)
    call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd3 )
    Ja(3,1) = -half*( dd3(2) - dd2(3) );
    Ja(3,2) = -half*( dd1(3) - dd3(1) );
    Ja(3,3) = -half*( dd2(1) - dd1(2) );
  end function calc_grid_metrics

  pure function calc_grid_metrics_alt(this,point,X1,X2,X3) result(Ja)
    use set_constants, only : zero
    class(interpolant_t),       intent(in) :: this
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: Ja
    real(dp), dimension(size(X1,1),size(X1,2),size(X1,3),3) :: X
    real(dp), dimension(size(X1,1),size(X1,2),size(X1,3))   :: tmp
    real(dp), dimension(3) :: dX_l, dX_m, dd1, dd2, dd3
    real(dp) :: junk
    integer, dimension(3), parameter :: ijk = [1,2,3]
    integer, dimension(3), parameter :: kij = cshift(ijk,1)
    integer, dimension(3), parameter :: jki = cshift(kij,1)
    integer, dimension(3) :: Npts
    integer :: i
    Ja = zero
    X(:,:,:,1) = X1
    X(:,:,:,2) = X2
    X(:,:,:,3) = X3
    Npts = shape(X1)
    do i = 1,3
      associate( l => kij(i), m => jki(i) )
        associate( X_l => X(:,:,:,l), X_m => X(:,:,:,m) )
          call this%lagbary_3D_wgrad( point, X_l, Npts, junk, dX_l )
          call this%lagbary_3D_wgrad( point, X_m, Npts, junk, dX_m )
          tmp = X_l*dX_m(1) - X_m*dX_l(1)
          call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd1 )
          tmp = X_l*dX_m(2) - X_m*dX_l(2)
          call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd2 )
          tmp = X_l*dX_m(3) - X_m*dX_l(3)
          call this%lagbary_3D_wgrad( point, tmp, Npts, junk, dd3 )
          Ja(i,1) = -half*( dd3(2) - dd2(3) );
          Ja(i,2) = -half*( dd1(3) - dd3(1) );
          Ja(i,3) = -half*( dd2(1) - dd1(2) );
        end associate
      end associate
    end do
  end function calc_grid_metrics_alt

  pure function normal_vectors(this,point,X1,X2,X3) result(Nvec)
    use math, only : cross_product
    class(interpolant_t),       intent(in) :: this
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: Nvec
    Nvec = this%calc_grid_metrics(point,X1,X2,X3)
    Nvec(:,1) = Nvec(:,1)/norm2(Nvec(:,1))
    Nvec(:,2) = Nvec(:,2)/norm2(Nvec(:,2))
    Nvec(:,3) = Nvec(:,3)/norm2(Nvec(:,3))
  end function normal_vectors

  pure subroutine map_point_3D_curve(this,point,X1,X2,X3,xyz,dS)
    class(interpolant_t),   intent(in)  :: this
    real(dp), dimension(1), intent(in)  :: point ! [t]
    real(dp), dimension(:), intent(in)  :: X1, X2, X3
    real(dp), dimension(3), intent(out) :: xyz
    real(dp),               intent(out) :: dS
    integer,  dimension(1) :: Npts
    real(dp), dimension(3) :: dval
    Npts = shape(X1)
    call this%lagbary_wderiv(point(1),1,X1,Npts,xyz(1),dval(1))
    call this%lagbary_wderiv(point(1),1,X2,Npts,xyz(2),dval(2))
    call this%lagbary_wderiv(point(1),1,X3,Npts,xyz(3),dval(3))
    dS = norm2(dval)
  end subroutine map_point_3D_curve

  pure subroutine map_point_3D_surface(this,point,X1,X2,X3,xyz,dA)
    use math, only : cross_product
    class(interpolant_t),     intent(in)  :: this
    real(dp), dimension(2),   intent(in)  :: point ! [u,v]
    real(dp), dimension(:,:), intent(in)  :: X1, X2, X3
    real(dp), dimension(3),   intent(out) :: xyz
    real(dp),                 intent(out) :: dA
    integer,  dimension(2) :: Npts
    real(dp), dimension(2) :: tmp
    real(dp), dimension(3) :: drdu, drdv
    Npts = shape(X1)
    call this%lagbary_2D_wgrad(point,X1,Npts,xyz(1),tmp)
    drdu(1) = tmp(1); drdv(1) = tmp(2)
    call this%lagbary_2D_wgrad(point,X2,Npts,xyz(2),tmp)
    drdu(2) = tmp(1); drdv(2) = tmp(2)
    call this%lagbary_2D_wgrad(point,X3,Npts,xyz(3),tmp)
    drdu(3) = tmp(1); drdv(3) = tmp(2)
    dA = norm2( cross_product(drdu,drdv) )
  end subroutine map_point_3D_surface

  pure subroutine map_point_3D_volume(this,point,X1,X2,X3,xyz,dV)
    use math, only : det_3x3
    class(interpolant_t),       intent(in)  :: this
    real(dp), dimension(3),     intent(in)  :: point ! [xi,eta,zeta]
    real(dp), dimension(:,:,:), intent(in)  :: X1, X2, X3
    real(dp), dimension(3),     intent(out) :: xyz
    real(dp),                   intent(out) :: dV
    integer, dimension(3) :: Npts
    real(dp), dimension(3,3) :: A
    Npts = shape(X1)
    call this%lagbary_3D_wgrad(point,X1,Npts,xyz(1),A(:,1))
    call this%lagbary_3D_wgrad(point,X2,Npts,xyz(2),A(:,2))
    call this%lagbary_3D_wgrad(point,X3,Npts,xyz(3),A(:,3))
    dV = det_3x3(A)
  end subroutine map_point_3D_volume

end module barycentric_interpolant_derived_type

module quadrature_derived_type

  use set_precision,       only : dp
  use set_constants,       only : zero
  implicit none
  private
  public :: quad_t
  public :: quad_ptr, quad_ptr_3D
  public :: create_quad_ref_1D, create_quad_ref_2D, create_quad_ref_3D
  ! public :: map_quad_ref_to_physical_1D
  ! public :: map_quad_ref_to_physical_2D
  ! public :: map_quad_ref_to_physical_3D
  type quad_t
    integer :: n_quad = 0
    real(dp), allocatable, dimension(:,:) :: quad_pts
    real(dp), allocatable, dimension(:)   :: quad_wts
  contains
    private
    procedure, public, pass :: create  => allocate_quad
    procedure, public, pass :: destroy => deallocate_quad
    generic,   public :: integrate => integrate_scalar, integrate_vector
    procedure :: integrate_scalar
    procedure :: integrate_vector
  end type quad_t

  type quad_ptr
    type(quad_t), pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_quad_ptr
  end type quad_ptr

  type quad_ptr_3D
    type(quad_t), dimension(:,:,:), pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_quad_ptr_3D
  end type quad_ptr_3D

contains

  pure elemental subroutine destroy_quad_ptr_3D( this )
    class(quad_ptr_3D), intent(inout) :: this
    this%p => null()
  end subroutine destroy_quad_ptr_3D

  pure elemental subroutine destroy_quad_ptr( this )
    class(quad_ptr), intent(inout) :: this
    this%p => null()
  end subroutine destroy_quad_ptr

  pure elemental subroutine allocate_quad( this, n_quad )
    use set_constants, only : zero
    class(quad_t), intent(inout) :: this
    integer,       intent(in)    :: n_quad
    this%n_quad = n_quad
    allocate( this%quad_pts(3,n_quad) )
    this%quad_pts = zero
    allocate( this%quad_wts(n_quad) )
    this%quad_wts = zero
  end subroutine allocate_quad

  pure elemental subroutine deallocate_quad( this )
    class(quad_t), intent(inout) :: this
    this%n_quad = 0
    if( allocated( this%quad_wts  ) ) deallocate( this%quad_wts  )
    if( allocated( this%quad_pts  ) ) deallocate( this%quad_pts  )
  end subroutine deallocate_quad

  pure function integrate_scalar( this, f ) result( integral )
    use set_precision, only : dp
    class(quad_t),                    intent(in) :: this
    real(dp), dimension(this%n_quad), intent(in) :: f
    real(dp)                                     :: integral
    integral = dot_product(f,this%quad_wts)
  end function integrate_scalar

  pure function integrate_vector( this, neq, f ) result( integral )
    use set_precision, only : dp
    class(quad_t),                        intent(in) :: this
    integer,                              intent(in) :: neq
    real(dp), dimension(neq,this%n_quad), intent(in) :: f
    real(dp), dimension(neq)                         :: integral
    integer :: n
    do n = 1, neq
      integral(n) = dot_product(f(n,:),this%quad_wts)
    end do
  end function integrate_vector

  pure function gauss_1D_size( polynomial_order ) result( N_quad )
    use set_constants, only : half
    integer, intent(in) :: polynomial_order
    integer             :: N_quad
    N_quad = ceiling( half*(polynomial_order + 1) )
  end function gauss_1D_size

  pure subroutine gauss_1D( n_quad, pts_1D, wts_1D )
    use math, only : LegendreGaussNodesAndWeights
    integer,                       intent(in)  :: n_quad
    real(dp), dimension( n_quad ), intent(out) :: pts_1D
    real(dp), dimension( n_quad ), intent(out) :: wts_1D
    call LegendreGaussNodesAndWeights(n_quad-1, pts_1D, wts_1D)
  end subroutine gauss_1D

  pure subroutine create_quad_ref_1D( quad_order, quad_ref )
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: xtmp
    integer :: n_quad
    n_quad = gauss_1D_size( quad_order )
    call quad_ref%destroy()
    call quad_ref%create( n_quad )
    call gauss_1D( n_quad, xtmp, quad_ref%quad_wts )
    quad_ref%quad_pts(1,:) = xtmp
  end subroutine create_quad_ref_1D

  pure subroutine create_quad_ref_2D( quad_order, quad_ref )
    use set_constants, only : zero
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    integer :: n_quad
    integer :: i, j, cnt
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: pts_1D
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: wts_1D
    n_quad = gauss_1D_size( quad_order )
    call gauss_1D(n_quad, pts_1D, wts_1D)
    call quad_ref%destroy()
    call quad_ref%create( n_quad**2 )
    cnt = 0
    do j = 1, n_quad
      do i = 1, n_quad
        cnt = cnt + 1
        quad_ref%quad_pts(:,cnt) = [ pts_1D(i), pts_1D(j), zero ]
        quad_ref%quad_wts(cnt) = wts_1D(i)*wts_1D(j)
      end do
    end do
  end subroutine create_quad_ref_2D

  pure subroutine create_quad_ref_3D( quad_order, quad_ref )
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    integer :: n_quad
    integer :: i, j, k, cnt
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: pts_1D
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: wts_1D
    n_quad = gauss_1D_size( quad_order )
    call gauss_1D(n_quad, pts_1D, wts_1D)
    call quad_ref%destroy()
    call quad_ref%create( n_quad**3 )
    cnt = 0
    do k = 1, n_quad
      do j = 1, n_quad
        do i = 1, n_quad
          cnt = cnt + 1
          quad_ref%quad_pts(:,cnt) = [ pts_1D(i), pts_1D(j), pts_1D(k) ]
          quad_ref%quad_wts(cnt) = wts_1D(i)*wts_1D(j)*wts_1D(k)
        end do
      end do
    end do
  end subroutine create_quad_ref_3D

  pure subroutine map_quad_ref_to_physical_1D( x1nodes, x2nodes, x3nodes, mask, dir,    &
                                               quad_ref, interpolant, quad_physical )
    use set_constants,           only : zero, half
    use math,                    only : vector_norm
    use barycentric_interpolant_derived_type, only : interpolant_t
    ! use reshape_array, only : extract_2D_slice_from_3D_array

    real(dp), dimension(:,:,:), intent(in)  :: x1nodes, x2nodes, x3nodes
    logical,  dimension(:,:,:), intent(in)  :: mask
    integer,                    intent(in)  :: dir
    type(quad_t),               intent(in)  :: quad_ref
    type(interpolant_t),        intent(in)  :: interpolant
    type(quad_t),               intent(out) :: quad_physical
    integer  :: n, i
    integer, dimension(1) :: Npts
    real(dp) :: dS
    real(dp), dimension(3) :: node_diff

    real(dp), dimension(3) :: tangent, normal, binormal
    real(dp) :: kappa

    continue
    Npts = shape(x1nodes)

    if( quad_ref%n_quad /= quad_physical%n_quad ) then
      call quad_physical%destroy()
      call quad_physical%create(quad_ref%n_quad)
    end if

    quad_physical%quad_pts = zero
    do n = 1,quad_ref%n_quad
      ! call interpolant%map_point(quad_ref%quad_pts(dir,n), x1nodes, x2nodes, x3nodes, quad_physical%quad_pts(1,n)
      quad_physical%quad_wts(n) = det_jac * quad_ref%quad_wts(n)

    end do
  end subroutine map_quad_ref_to_physical_1D

  ! pure subroutine map_quad_ref_to_physical_2D( x1nodes, x2nodes, x3nodes,    &
  !                                                 quad_ref, quad_physical )
  !   use set_constants,           only : zero
  !   use lagrange_interpolation,  only : lagbary_2D, jacobian_determinant

  !   real(dp),     dimension(:,:), intent(in)  :: x1nodes, x2nodes, x3nodes
  !   type(quad_t),                 intent(in)  :: quad_ref
  !   type(quad_t),                 intent(out) :: quad_physical
  !   integer  :: n
  !   real(dp) :: det_jac
  !   integer,      dimension(2) :: Npts
  !   continue

  !   Npts = shape(x1nodes)
  !   if( quad_ref%n_quad /= quad_physical%n_quad ) then
  !     call quad_physical%destroy()
  !     call quad_physical%create(quad_ref%n_quad)
  !   end if

  !   quad_physical%quad_pts = zero
  !   do n = 1,quad_ref%n_quad
  !     call lagbary_2D( quad_ref%quad_pts(1:2,n), x1nodes, Npts, &
  !                 quad_physical%quad_pts(  1,n) )
  !     call lagbary_2D( quad_ref%quad_pts(1:2,n), x2nodes, Npts, &
  !                 quad_physical%quad_pts(  2,n) )
  !     call lagbary_2D( quad_ref%quad_pts(1:2,n), x3nodes, Npts, &
  !                 quad_physical%quad_pts(  3,n) )
  !     det_jac = jacobian_determinant( quad_ref%quad_pts(1:2,n), &
  !                 x1nodes, x2nodes, x3nodes )
  !     quad_physical%quad_wts(n) = abs(det_jac) * quad_ref%quad_wts(n)
  !   end do
  ! end subroutine map_quad_ref_to_physical_2D

  ! pure subroutine map_quad_ref_to_physical_3D( x1nodes, x2nodes, x3nodes,    &
  !                                                 quad_ref, quad_physical )
  !   use set_constants,           only : zero
  !   use lagrange_interpolation,  only : lagbary_3D, jacobian_determinant
  !   real(dp), dimension(:,:,:), intent(in)  :: x1nodes, x2nodes, x3nodes
  !   type(quad_t),               intent(in)  :: quad_ref
  !   type(quad_t),               intent(out) :: quad_physical
  !   integer  :: n
  !   real(dp) :: det_jac
  !   integer,  dimension(3) :: Npts
  !   continue

  !   Npts = shape(x1nodes)
  !   if( quad_ref%n_quad /= quad_physical%n_quad ) then
  !     call quad_physical%destroy()
  !     call quad_physical%create(quad_ref%n_quad)
  !   end if

  !   quad_physical%quad_pts = zero
  !   do n = 1,quad_ref%n_quad
  !     call lagbary_3D( quad_ref%quad_pts(:,n), x1nodes, Npts, &
  !                     quad_physical%quad_pts(  1,n) )
  !     call lagbary_3D( quad_ref%quad_pts(:,n), x2nodes, Npts, &
  !                     quad_physical%quad_pts(  2,n) )
  !     call lagbary_3D( quad_ref%quad_pts(:,n), x3nodes, Npts, &
  !                     quad_physical%quad_pts(  3,n) )
  !     det_jac = jacobian_determinant( quad_ref%quad_pts(:,n), &
  !                                     x1nodes, x2nodes, x3nodes )
  !     quad_physical%quad_wts(n) = det_jac * quad_ref%quad_wts(n)
  !   end do
  ! end subroutine map_quad_ref_to_physical_3D

  ! pure subroutine map_cell_quad_points( dim, n_skip, quad_order, coords, mask, ref_quads, quad )
  !   use quadrature_derived_type, only : map_quad_ref_to_physical_1D,           &
  !                                       map_quad_ref_to_physical_2D,           &
  !                                       map_quad_ref_to_physical_3D
  !   integer,                                                        intent(in)    :: quad_order, dim
  !   integer,      dimension(3),                                     intent(in)    :: n_skip
  !   real(dp),     dimension(n_skip(1)+1,n_skip(1)+1,n_skip(1)+1,3), intent(in)    :: coords
  !   logical,      dimension(n_skip(1)+1,n_skip(1)+1,n_skip(1)+1),   intent(in)    :: mask
  !   type(quad_t), dimension(3),                                     intent(in)    :: ref_quads
  !   type(quad_t),                                                   intent(inout) :: quad
  !   real(dp),     dimension(product(n_skip+1)) :: Xtmp, Ytmp, Ztmp
  !   integer :: n_mask, n_quad

  !   n_mask = count(mask)

  !   n_quad = ref_quads(1)%n_quad

  !   Xtmp(1:n_mask) = pack(coords(:,:,:,1),mask)
  !   Ytmp(1:n_mask) = pack(coords(:,:,:,2),mask)
  !   Ztmp(1:n_mask) = pack(coords(:,:,:,3),mask)

  !   select case(dim)
  !   case(1)
  !     call map_quad_ref_to_physical_1D( Xtmp(1:n_mask), &
  !                                       Ytmp(1:n_mask), &
  !                                       Ztmp(1:n_mask), &
  !                                       ref_quads(1), quad )
  !   case(2)
  !     call map_quad_ref_to_physical_2D( reshape( Xtmp(1:n_mask), [n_quad, n_quad ] ), &
  !                                       reshape( Ytmp(1:n_mask), [n_quad, n_quad ] ), &
  !                                       reshape( Ztmp(1:n_mask), [n_quad, n_quad ] ), &
  !                                       ref_quads(2), quad )
  !   case(3)
  !     call map_quad_ref_to_physical_3D( reshape( Xtmp(1:n_mask), [n_quad, n_quad, n_quad ] ), &
  !                                       reshape( Ytmp(1:n_mask), [n_quad, n_quad, n_quad ] ), &
  !                                       reshape( Ztmp(1:n_mask), [n_quad, n_quad, n_quad ] ), &
  !                                       ref_quads(3), quad )
  !   end select
  ! end subroutine map_cell_quad_points
end module quadrature_derived_type

module reshape_array
  use set_precision, only : dp
  implicit none
  private
  public :: extract_2D_slice_from_3D_array
contains


  pure function extract_2D_slice_from_3D_array(A,mask,sz) result(slice)
    use set_constants, only : zero
    real(dp), dimension(:,:,:), intent(in) :: A
    logical,  dimension(:,:,:), intent(in) :: mask
    integer,  dimension(2),     intent(in) :: sz
    real(dp), dimension(sz(1),sz(2)) :: slice, field
    logical,  dimension(sz(1),sz(2)) :: mask2
    field = zero
    mask2 = .true.
    slice = unpack( pack( A, mask ), mask2, field )
  end function extract_2D_slice_from_3D_array

end module reshape_array

program main
  use barycentric_interpolant_derived_type, only : interpolant_t

  type(interpolant_t) :: interp

  interp = interpolant_t(10)

  write(*,*) 'Here'

  call interp%destroy()

end program main