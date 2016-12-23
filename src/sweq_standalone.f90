

SUBROUTINE rsweqintegrate(h,u,r, elev, nsteps, kens, ndim, noiseswitch, kh, ku, kr, &
                            noise_freq, alpha, beta, h_c, h_r, uniforms1, uniforms2)
    IMPLICIT NONE
    INTEGER, PARAMETER:: dp=kind(0.d0)   
    REAL(dp), PARAMETER :: g = 10
    REAL(dp), PARAMETER :: pi = 3.141592653589793_dp

    !number of ensemble, size of the domain (in # of grid points) and number of time steps
    INTEGER :: kens, ndim, nsteps       

    ! input fields h, u and r:
    ! If you pass an ensemble, each member should be appended together.
    ! I implemented it in this way for convenience with my R interface, but you can
    ! easily redefine the input as matrices with the dimension as for the 
    ! subaltern fields below.
    REAL(dp), DIMENSION(ndim*kens)  :: h, u, r 
    
    ! Dummy fields arrays for the ensembles: new and m are for the update and the RAW filter
    REAL(dp), DIMENSION(ndim, kens) :: h_a, u_a, r_a, hnew, unew, rnew, elev_array, hm, um, rm
    REAL(dp), DIMENSION(ndim, kens) :: phi, beta_mask    !for rain dynamic
    REAL(dp), DIMENSION(ndim, kens) :: unoise            !for random perturbations

    REAL(dp), DIMENSION(ndim)       :: elev ! not tested for value different than uniformly 0
    INTEGER, DIMENSION(ndim)        :: fwd, bwd !some dummy variable to implement shifting operators

    ! random noise:
    ! the uniforms1 are for the random plume positions in the spatial domain
    ! the uniforms2 are for the probability to get a noise term in this time step
    REAL(dp), DIMENSION(kens*nsteps)  :: uniforms1, uniforms2 !the necessary random numbers (passed from outside)
    REAL(dp), DIMENSION(nsteps, kens)        :: uniforms_a, uniforms_b ! the same but rearranged
    
    INTEGER :: nx, i, j, l 
    INTEGER :: noiseswitch  !0=no noise, 1=noise
    REAL(dp) :: noise_freq  !probability to get a noise term in the step (often set to 1)

    !for noise:
    REAL(dp) :: pos, sig, temp, znoise(ndim)
    REAL(dp) :: amplitude= 0.005_dp

    !for RAW filter:
    REAL(dp), DIMENSION(ndim) :: d
 
    !diffusion parameters: 
    REAL(dp) :: kh, ku, kr ! typically: kh=25000, ku=25000, kr=200

    REAL(dp) :: gam        !influence of rain on the potential phi
    REAL(dp) :: dx, dt     !increments (could be passed as arguments)
    REAL(dp) :: alpha      !terminal velocity factor (default = 1/4000.)
    REAL(dp) :: beta       !rain formation factor    (default = 3.)
    REAL(dp) :: h_c        !threshold for updraft (Wuersch = 90.02)
    REAL(dp) :: h_r        !critical height for the building of rain (Wuersch = 90.4)
    REAL(dp) :: hmax       !Value for the replaced height in the calculation of the h grad. 89.977
    REAL(dp) :: phic       !potential



    ! stack ensemble vectors in arrays:
    DO j = 1,kens
        h_a(:,j) = h( ((j-1)*ndim + 1) : j*ndim)
        u_a(:,j) = u( ((j-1)*ndim + 1) : j*ndim)
        r_a(:,j) = r( ((j-1)*ndim + 1) : j*ndim)
    END DO

    DO j = 1,kens
        uniforms_a(:,j) = uniforms1( ((j-1)*nsteps + 1) : j*nsteps)
        uniforms_b(:,j) = uniforms2( ((j-1)*nsteps + 1) : j*nsteps)
    END DO


    ! some parameter values (could be passed as argument if you want to change)
    gam = 1._dp
    hmax = 89.977_dp
    phic = hmax * g

    !increments:
    dx=500._dp  ! resolution = 500 meters
    dt=1._dp

    !number of cells
    nx = ndim

    !duplicate elev kens times:
    elev_array = SPREAD(elev, 2, kens)


    fwd = (/ (/ (i, i=2,nx) /), 1 /)
    bwd = (/ nx,  (/ (i, i=1,(nx-1)) /) /)  

    !set initial intermediary fields (for RAW)
    hm = h_a
    um = u_a
    rm = r_a

    DO i=1, nsteps
        ! Transform h into potential:
        phi = g * (h_a + elev_array)                                    !without updraft
        WHERE ( h_a + elev_array > h_c ) phi = phic + g * elev_array    !updraft over h_cloud

        phi = phi + gam*r_a         !influence of rain

        !random perturbations:
        unoise = 0._dp
        if (noiseswitch == 1) THEN
            DO j=1,kens !loop through ensemble members
                sig = 4.
                pos = uniforms_a(i, j)
                pos = pos * nx

                DO l = 1,nx
                    znoise(l) = exp( -0.5_dp * ( ( ABS( REAL(nx, dp)/2 - ABS( ABS(l - pos) - REAL(nx, dp)/2 ) ) )**2 / sig**2 ) )
                END DO

                znoise = znoise - znoise(bwd)
                znoise = znoise/maxval(znoise)
                unoise(:,j) = unoise(:,j) + amplitude * znoise
                
                !add noise only in 'prob' cases:
                temp = uniforms_b(i, j)
                IF (temp <= noise_freq) THEN
                    u_a(:,j) = u_a(:,j) + unoise(:,j)
                END IF
            END DO
        END IF

        !--------------------------------------------------------------
        !dynamic:
        unew = um - (dt/2/dx) * ( u_a(fwd,:)**2 - u_a(bwd,:)**2 ) &
                - (dt*2/dx) * ( phi - phi(bwd,:) ) &
                + (ku*dt/4/(dx*dx)) * ( u_a(fwd,:) - 2*u_a + u_a(bwd,:) )


        hnew = hm - dt/dx * ( u_a(fwd,:) * ( h_a + h_a(fwd,:) )  - u_a * ( h_a(bwd,:) + h_a ) )  &
            + (kh*dt/4/(dx*dx)) * ( h_a(fwd,:) - 2*h_a + h_a(bwd,:) )

        beta_mask = 0._dp
        WHERE ( (h_a + elev_array) > h_r  .AND. ( u_a(fwd,:) - u_a ) < 0._dp) beta_mask = beta


        rnew = rm - 2 * alpha * dt * r_a &
            - 2 * beta_mask * dt/dx * ( u_a(fwd,:) - u_a ) &
            - dt/2/dx * ( u_a + u_a(fwd,:)) * (r_a(fwd,:) - r_a(bwd,:) ) &
            + (kr*dt/4/(dx*dx)) * ( r_a(fwd,:) - 2*r_a + r_a(bwd,:) ) 

        !--------------------------------------------------------------
        !periodic boundary conditions:
        DO j=1,kens !loop through ensemble members
            unew(1,j) = unew(nx-1,j)
            unew(nx,j) = unew(2,j)
            
            hnew(1,j) = hnew(nx-1,j)
            hnew(nx,j) = hnew(2,j)

            rnew(1,j) = rnew(nx-1,j)
            rnew(nx,j) = rnew(2,j)
        END DO

        WHERE(r_a < 0) r_a=0

        !--------------------------------------------------------------
        !RAW time smoother:
        DO j=1,kens  !loop through ensemble members
            d = .1_dp * .5_dp * ( unew(:,j) - 2._dp *  u_a(:,j) + um(:,j))
            um(:,j) =  u_a(:,j) + 0.53_dp * d
            u_a(:,j) = unew(:,j) - 0.47_dp * d

            d = .1_dp * .5_dp * ( hnew(:,j) - 2._dp *  h_a(:,j) + hm(:,j))
            hm(:,j) =  h_a(:,j) + 0.53_dp * d
            h_a(:,j) = hnew(:,j) - 0.47_dp * d

            d = .1_dp * .5_dp * ( rnew(:,j) - 2._dp *  r_a(:,j) + rm(:,j))
            rm(:,j) =  r_a(:,j) + 0.53_dp * d
            r_a(:,j) = rnew(:,j) - 0.47_dp * d
        END DO



    END DO


!--------------------------------------------------------------
    ! put back the array in big vectors:
    DO j = 1,kens
        h( ((j-1)*ndim + 1) : j*ndim) = h_a(:,j)
        u( ((j-1)*ndim + 1) : j*ndim) = u_a(:,j)
        r( ((j-1)*ndim + 1) : j*ndim) = r_a(:,j)
    END DO

END SUBROUTINE
