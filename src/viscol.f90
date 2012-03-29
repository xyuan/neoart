!--------------------------------------------------------------------
!> This subroutine calculates the energy dependent collision
!> frequency and the energy dependent energy scattering.
!>
!> input    nsm  : the maximum number of species
!>          ncm  : the maximum number of charge states
!>          ns   : the actual number of species
!>          is   : the species number for which the collision
!>                 frequency is to be calculated
!>          ic   : the charge state of the species for which
!>                 the collision frequency is to be calculated.
!>          tau  : normalized collision frequency (n_i m_i /
!>                 tau_ij)  array(nsm,nsm)
!>          m    : array(nsm) the mass of the species in kg
!>          t    : array(nsm) the temperature of the species
!>                 in keV
!>          xi   : array(nsm,ncm) the relative weight of every
!>                 charge state
!>          nc   : array(nsm) the number of charge states for
!>                 every species
!>          x    : normalized (to thermal) velocity for which
!>                 the energy dependent collision frequency is
!>                 to be calculated
!>          den  : array(nsm,ncm) density of every component
!>                 in 10^19 m^-3
!> output   nud  : the pitch angle scattering frequency
!>          nue  : the energy scattering frequency
!>
!--------------------------------------------------------------------
subroutine viscol(nsm, ncm, ns, is, ic, tau, m, t, xi, nc, x, &
                   & den, nud, nue)

  implicit none

  integer nsm, ncm, is, ic, nc, i, ns
  real tau, m, t, xi, x, nud, nue, den
  real ph, g, pi, xab
  dimension tau(nsm,nsm), m(nsm), t(nsm), xi(nsm, ncm), nc(nsm), &
          & den(nsm, ncm)

  pi = 4. * atan(1.)

  nud = 0.
  nue = 0.

  ! the loop over species
  do i = 1, ns
    ! calculate xab = vthb / vtha
    xab = sqrt(m(is) * t(i) / (m(i) * t(is)))

    ph = erf(x / xab)
    g  = (ph - 2 * x * exp(-(x / xab)**2) / (xab * sqrt(pi))) &
    &    / (2 * (x / xab)**2)

    nud = nud + tau(is,i) * ( ph - g ) / x**3
    nue = nue + tau(is,i) * (-2.*ph/x**3 + 4.*(t(is)/t(i) &
    &       + 1./xab**2) *g/x)
  enddo

  nud = nud * 0.75 * sqrt(pi) * xi(is,ic) / (1.e19*den(is,ic) &
  &     * m(is))
  nue = nue * 0.75 * sqrt(pi) * xi(is,ic) / (1.e19*den(is,ic) &
  &     * m(is))

  return
end
