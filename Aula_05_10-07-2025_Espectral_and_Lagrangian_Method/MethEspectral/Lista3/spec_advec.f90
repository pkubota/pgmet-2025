program advecao_espectral_manual
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: N = 64
  real(dp), parameter :: pi = 3.141592653589793d0
  real(dp), parameter :: L = 2.0d0 * pi, dx = L / N
  real(dp), parameter :: c = 1.0d0
  real(dp), parameter :: dt = 0.01d0, tmax = 4.0d0
  integer, parameter :: ic=1
  integer :: i, j, k, nsteps,irec,lrec
  real(dp) :: x(N), u(N),um(N),up(N),pdfdt(N),t
  real     :: us(N)
  complex(dp), parameter :: II = (0.0_dp, 1.0_dp)
  real(dp) :: k_norm
  complex(dp), dimension(N) :: filter_coeff ! Filter coefficients
  complex(dp), dimension(N) :: k_wn       ! Real part of wavenumbers

  ! Filter parameters
  double precision, parameter :: filter_alpha = 36.0d0! Controls filter strength (higher = stronger)
  double precision, parameter :: filter_s     = 106.0d0! Controls filter order    (higher = sharper cutoff)
  real(dp),parameter :: gamma=0.1_dp  ! filter const.

  inquire(iolength=lrec)us
  open(1,file='espectral_method.bin',ACCESS='direct',FORM='unformatted',action='write',&
       status='unknown',recl=lrec)


  ! Initial condition 
  if(ic==1)then
  ! Initial condition (e.g., a seno function )
     do i = 1, N
        x(i) = (i - 1) * dx
        u(i) = sin(x(i))
     end do
  else if(ic==2)then
  ! Initial condition (e.g., a square wave)
     do i = 1, N
        x(i) = (i - 1) * L / N
       if (x(i) >= L/4.0d0 .and. x(i) <= L/2.0d0) then
          u(i) = 1.0d0
       else
          u(i) = 0.0d0
       end if
     end do
  else
  ! Inicializa domínio
     do i = 1, N
        x(i) = (i - 1) * dx
        u(i) = sin(x(i))
     end do
  endif

  um=u

  ! spectral filter (exponential filter)

  call  filter_epectral()
  
  irec=1
  us=u
  write(1,rec=irec) us

  ! Loop temporal
  t = 0.0_dp
  nsteps = int(tmax / dt)
  do k = 1, nsteps
      irec=irec+1

      call TendSpectral(pdfdt,N)

      call leapfrog(pdfdt,N)
      us=u
      write(1,rec=irec) us

      t = t + dt
  end do

  ! Output final (print simples)
  print *, "x      u(x)"
  do i = 1, N
     print *, x(i), u(i)
  end do

contains

  subroutine leapfrog(pdfdt,N)
    integer , intent(in ) :: N
    real(dp), intent(in ) :: pdfdt(N)
    integer :: i
      ! Atualiza solução com esquema de Euler
      do i = 1, N
!     a) add tendency to old (t-1) (filtered) value
      
         up(i)=um(i)+ 2*dt*pdfdt(i)
!
!     b) compute the new (filtered) (t-1) value 
!
         um(i)=u(i)+gamma*(um(i)-2.*u(i)+up(i))
!
!     c) move the (t+1) value to the new (t) value
!
         u(i)=up(i)
      end do

  end subroutine leapfrog


  subroutine TendSpectral(pdfdt,N)
    integer , intent(in ) :: N
    real(dp), intent(out ) :: pdfdt(N)
    complex(dp) :: u_hat(N)
    complex(dp) :: dudx_hat(N)
    real(dp)    :: du(N)
    integer     :: i,j

      ! FFT direta: u(x) => u_hat(k)
      call fft_manual_forward(u, u_hat, N)

      ! Apply exponential filter
      do j = 1, N
         u_hat(j) = u_hat(j) * filter_coeff(j)
      end do

      ! Derivada espectral: du/dx = i * k * u_hat
      do j = 1, N
        if (j-1 <= N/2) then
           dudx_hat(j) = II * real(j-1    , dp) * u_hat(j)
        else
           dudx_hat(j) = II * real(j-1 - N, dp) * u_hat(j)
        end if
      end do

      ! FFT inversa: u'(x) => du_hat(k)

      call fft_manual_backward(dudx_hat, du, N)

      do i = 1, N
         pdfdt(i) = - c * du(i) ! input tendency
      end do
  end subroutine TendSpectral

  subroutine fft_manual_forward(f, f_hat, N)
    real(dp), intent(in) :: f(N)
    complex(dp), intent(out) :: f_hat(N)
    integer, intent(in) :: N
    integer :: k, nn
    complex(dp) :: sum, w

    do k = 0, N-1
       sum = (0.0_dp, 0.0_dp)
       do nn = 0, N-1
          w = exp(-2.0_dp * pi * II * real(k * nn, dp) / real(N, dp))
          sum = sum + f(nn+1) * w
       end do
       f_hat(k+1) = sum / real(N, dp)
    end do
  end subroutine

  subroutine fft_manual_backward(f_hat, f, N)
    complex(dp), intent(in) :: f_hat(N)
    real(dp), intent(out) :: f(N)
    integer, intent(in) :: N
    integer :: nn, k
    complex(dp) :: sum, w

    do nn = 0, N-1
       sum = (0.0_dp, 0.0_dp)
       do k = 0, N-1
          w = exp(2.0_dp * pi * II * real(k * nn, dp) / real(N, dp))
          sum = sum + f_hat(k+1) * w
       end do
       f(nn+1) = real(sum)
    end do
  end subroutine

  subroutine filter_epectral()
    ! Initialize wavenumbers
    k_wn=0.0d0
    ! For N even, k = [0, 1, ..., N/2-1, -N/2, ..., -1]
    ! For this spectral method, we need k_real and k_imag for complex multiplication
    do j = 1, N/2 - 1
      k_wn(j        ) = CMPLX( 0.0d0, (j - 1) * 2.0 * 3.141592653589793d0 / L)
      k_wn(N - j + 1) = CMPLX( 0.0d0,(-j   ) * 2.0 * 3.141592653589793d0 / L)
    end do
    k_wn(N/2)     = CMPLX( 0.0d0,(N/2 - 1) * 2.0 * 3.141592653589793d0 / L)
    k_wn(N/2 + 1) = CMPLX( 0.0d0,(N/2    ) * 2.0 * 3.141592653589793d0 / L) ! Nyquist frequency for even N
    ! Calculate filter coefficients (exponential filter)
    do i = 1, N
      k_norm = abs(k_wn(i)) / (N/2.0d0 * 2.0 * 3.141592653589793d0 / L) ! Normalized wavenumber
      filter_coeff(i) = exp(-filter_alpha * (k_norm**filter_s))
    end do
  end subroutine filter_epectral

end program
