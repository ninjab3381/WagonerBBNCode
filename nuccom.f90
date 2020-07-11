!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     Changes (to run inder DEC unix f77): 
!     -----------------------------------
!     COMMON /bessel/ -> COMMON /besselcb/
!     COMMON /therm/  -> COMMON /thermcb/
!     COMMON /time/   -> COMMON /ttime/
!     ir=1 -> ir=5
!     iw=1 -> iw=6
!     All `entry' routines removed

!========================IDENTIFICATION DIVISION================================

SUBROUTINE driver

    !----------LINKAGES.
    !     CALLED BY - [subroutine] run
    !     CALLS     - [subroutine] start, derivs, accum

    !----------REMARKS.
    !     Runge-Kutta computational routine

    !----------PARAMETERS.
    PARAMETER (nvar = 29)          !Number of variables to be evolved.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.
    PARAMETER (cl = 1.e-16)        !Lower limit on size of time step.

    !----------COMMON AREAS.
    COMMON /evolp1/ t9, hv, phie, y                   !Evolution parameters.
    COMMON /evolp2/ dt9, dhv, dphie, dydt             !Evolution parameters.
    COMMON /evolp3/ t90, hv0, phie0, y0               !Evolution parameters.
    COMMON /compr/  cy, ct, t9i, t9f, ytmin, inc        !Computation parameters.
    COMMON /ttime/   t, dt, dlt9dt                    !Time variables.
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags,counters.
    COMMON /tcheck/  itime                          !Computation location.
    COMMON /runopt/ irun, isize, jsize               !Run options.


    !==========================DECLARATION DIVISION=================================

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature (in units of 10**9 K).
    REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
    REAL    phie                 !Chemical potential for electron.
    REAL    y(nnuc)              !Relative number abundances.

    !----------EVOLUTION PARAMETERS (DERIVATIVES).
    REAL    dydt(nnuc)           !Change in rel number abundances.

    !----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
    REAL    y0(nnuc)             !Rel # abund at beginning of iteration.

    !----------COMPUTATION PARAMETERS.
    REAL    cy                   !Time step limiting constant on abundances.
    REAL    ct                   !Time step limiting constant on temperature.
    REAL    t9f                  !Final temperature (in 10**9 K).
    REAL    ytmin                !Smallest abundances allowed.
    INTEGER inc                  !Accumulation increment.

    !----------TIME AND TIME STEP VARIABLES.
    REAL    t                    !Time.
    REAL    dt                   !Time step.
    REAL    dlt9dt               !(1/t9)*d(t9)/d(t).

    !----------COUNTERS AND FLAGS.
    INTEGER loop                 !Counts which Runge-Kutta loop.
    INTEGER ltime                !Indicates termination status.
    INTEGER is                   !# total time steps for particular run.
    INTEGER ip                   !# time steps after outputting a line.

    !----------COMPUTATION LOCATION.
    INTEGER itime                !Time check.

    !----------RUN OPTION.
    INTEGER isize                !Number of nuclides in computation.

    !----------TIME AND TIME STEP VARIABLES.
    REAL    dtmin                !Mininum time step.
    REAL    dtl                  !Time step from limitation on abund changes.

    !----------LABELS FOR VARIABLES TO BE TIME EVOLVED.
    INTEGER mvar                 !Total number of variables to be evolved.
    REAL    v(nvar)              !Variables to be time evolved.
    REAL    dvdt(nvar)           !Time derivatives.
    REAL    v0(nvar)             !Value of variables at original point.
    REAL    dvdt0(nvar)          !Value of derivatives at original point.

    !----------EQUIVALENCE STATEMENTS.
    EQUIVALENCE (v(4), y(1)), (dvdt(4), dydt(1)), (v0(4), y0(1))


    !===========================PROCEDURE DIVISION==================================

    !10--------INPUT INITIALIZATION INFORMATION, RELABEL----------------------------

    ltime = 0                    !Set termination indicator to zero.
    CALL start                   !Input initialization information.
    mvar = isize + 3            !Total number of variables to be evolved.

    !20--------LOOP ONE-------------------------------------------------------------

    200  continue                     !Begin Runge-Kutta looping.
    loop = 1                     !Loop indicator.
    !..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
    CALL derivs(loop)
    itime = 4                    !Time = 1st R-K loop.
    CALL check                   !Check interface subroutine.
    !..........ACCUMULATE.
    IF ((t9.le.t9f).or.                         !Low temp.&
        (dt.lt.abs(cl/dlt9dt)).or.              !Small dt.&
        (ip.eq.inc)) CALL accum                 !Enough iterations.
        !..........POSSIBLY TERMINATE COMPUTATION.
        IF (ltime.eq.1) THEN         !Return to run selection.
            RETURN
        END IF
        !..........RESET COUNTERS.
        IF (ip.eq.inc) THEN          !Reset iteration counters.
            ip = 0
        END IF
        ip = ip + 1
        is = is + 1
        !..........ADJUST TIME STEP.
        IF (is.gt.3) THEN            !Adjust time step after 3 iterations.
            dtmin = abs(1. / dlt9dt) * ct  !Trial value for minimum time step (Ref 1).
            DO i = 1, isize             !Go through all abundance changes.
                IF ((dydt(i).ne.0.).and.(y(i).gt.ytmin)) THEN
                    dtl = abs(y(i) / dydt(i)) * cy&
                            * (1. + (alog10(y(i)) / alog10(ytmin))**2)  !(Ref 2).
                    IF (dtl.lt.dtmin) dtmin = dtl         !Find smallest time step.
                END IF
            END DO
            IF (dtmin.gt.1.5 * dt) dtmin = 1.5 * dt       !Limit change in time step.
            dt = dtmin                 !Set new time step.
        END IF
        t = t + dt                   !Increment time.
        !..........STORE AND INCREMENT VALUES (Ref 3).
        DO i = 1, mvar
            v0(i) = v(i)
            dvdt0(i) = dvdt(i)
            v(i) = v0(i) + dvdt0(i) * dt
            IF ((i.ge.4).and.(v(i).lt.ytmin)) v(i) = ytmin  !Set at minimum value.
        END DO

        !30--------LOOP TWO-------------------------------------------------------------

        loop = 2                     !Step up loop counter.
        !..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
        CALL derivs(loop)
        itime = 7                    !Time = 2nd R-K loop.
        CALL check                   !Check interface subroutine.
        !..........INCREMENT VALUES.
        DO i = 1, mvar
            v(i) = v0(i) + .5 * (dvdt(i) + dvdt0(i)) * dt
            IF ((i.ge.4).and.(v(i).lt.ytmin)) v(i) = ytmin  !Set at minimum value.
        END DO
        GO TO 200

        !----------REFERENCES-----------------------------------------------------------
        !     1)  Constraint on dt from the requirement that
        !                (d(t9)/dt)*(dt/t9) < ct
        !           Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 293, equation C6.
        !     2)  Constraint on dt from
        !                dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2)
        !          Wagoner, R.V. 1969, page 293, equation C7 but with log term squared.
        !     3)  Wagoner, R.V. 1969, page 292, equations C1, C2.

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE start

    !----------LINKAGES.
    !     CALLED BY - [subroutine] driver
    !     CALLS     - [subroutine] rate1, bessel, rate0
    !               - [function] ex

    !----------REMARKS.
    !     Sets initial conditions.

    !----------PARAMETERS.
    PARAMETER (nrec = 88)          !Number of nuclear reactions.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.
    PARAMETER (const1 = 0.09615)   !Relation between time and temperature.
    PARAMETER (const2 = 6.6700e-8) !Gravitational constant.

    !----------COMMON AREAS.
    COMMON /rates/  f, r                            !Reaction rates.
    COMMON /evolp1/ t9, hv, phie, y                   !Evolution parameters.
    COMMON /evolp2/ dt9, dhv, dphie, dydt(nnuc)       !Evolution parameters.
    COMMON /evolp3/ t90, hv0, phie0, y0               !Evolution parameters.
    COMMON /compr/  cy, ct, t9i, t9f, ytmin, inc        !Computation parameters.
    COMMON /modpr/  g, tau, xnu, c, cosmo, xi           !Model parameters.
    COMMON /varpr/  dt1, eta1                       !Variational parameters.
    COMMON /ttime/   t, dt, dlt9dt                    !Time variables.
    COMMON /endens/ rhone0, rhob0, rhob, rnb          !Energy densities.
    COMMON /besselcb/ bl1, bl2, bl3, bl4, bl5, !Eval of function bl(z).&
    bm1,bm2, bm3, bm4, bm5, !Eval of function bm(z).&
    bn1,bn2, bn3, bn4, bn5            !Eval of function bn(z).
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags,counters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Neutrino parameters.
    COMMON /runopt/ irun, isize, jsize               !Run options.


    !==========================DECLARATION DIVISION=================================

    !----------REACTION RATES.
    REAL    f(nrec)              !Forward reaction rate coefficients.

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature (in units of 10**9 K).
    REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
    REAL    phie                 !Chemical potential of electron.
    REAL    y(nnuc)              !Relative number abundances.

    !----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
    REAL    y0(nnuc)             !Rel # abund at start of iteration.

    !----------COMPUTATION SETTINGS.
    REAL    t9i                  !Initial temperature (in 10**9 K).
    REAL    ytmin                !Smallest abundances allowed.
    INTEGER inc                  !Accumulation increment.

    !----------EARLY UNIVERSE MODEL PARAMETERS.
    REAL    g                    !Gravitational constant.
    REAL    tau                  !Neutron lifetime.
    REAL    xnu                  !Number of neutrino species.
    REAL    c(3)                 !c(1) is variation of grav. constant.&
    !c(2) is neutron lifetime (sec).&
    !c(3) is number of neutrino species.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------VARIATIONAL PARAMETERS.
    REAL    dt1                  !Initial time step.
    REAL    eta1                 !Baryon-to-photon ratio.

    !----------TIME VARIABLES.
    REAL    t                    !Time.
    REAL    dt                   !Time step.

    !----------ENERGY DENSITIES.
    REAL    rhone0               !Initial electron neutrino mass density.
    REAL    rhob0                !Initial baryon mass density.

    !----------EVALUATION OF FUNCTIONS bl,bm,bn.
    REAL    bl1, bl2, bl3, bl4, bl5  !Evaluation of function bl(z).

    !----------COUNTERS AND FLAGS.
    INTEGER ltime                !Indicates if output buffer printed.
    INTEGER is                   !# total time steps for particular run.
    INTEGER ip                   !# time steps after outputting a line.
    INTEGER it                   !# times accumulated in output buffer.
    INTEGER mbad                 !Indicates if gaussian elimination fails.

    !----------NEUTRINO PARAMETERS.
    REAL    tnu                  !Neutrino temperature.
    REAL    cnorm                !Normalizing constant.

    !----------RUN OPTION.
    INTEGER isize                !Number of nuclides in computation.

    !----------LOCAL VARIABLES.
    REAL    z                    !Defined by z = m(electron)*c**2/k*t9.


    !===========================PROCEDURE DIVISION==================================

    !10--------INITIALIZE FLAGS AND COUNTERS----------------------------------------

    ltime = 0                    !No output yet.
    is = 1                    !First iteration coming up.
    ip = inc                  !Set to maximum allowed # of iteration.
    it = 0                    !No accumulation yet.
    mbad = 0                    !No computational errors.

    !20--------SETTINGS-------------------------------------------------------------

    !..........COMPUTATIONAL SETTINGS.
    t9 = t9i                    !Initial temperature.
    tnu = t9                     !Initial neutrino temperature.
    t = 1 / (const1 * t9)**2       !Initial time (Ref 1).
    dt = dt1                    !Initial time step.
    !..........MODEL SETTINGS.
    g = const2 * c(1)            !Modify gravitational constant.
    tau = c(2)                   !Convert n half-life (min) to lifetime (sec).
    tau = tau / 0.98               !Coulomb correction (Ref 2).
    xnu = c(3)                   !Number of neutrino species.

    !30--------COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON--------------------

    IF ((15.011 / t9 + xi(1)).gt.58.) THEN      !Overabundance of antineutrinos.
        y(1) = 1.e-25              !Very little of neutrons.
        y(2) = 1.                  !Essentially all protons.
    ELSE
        IF ((15.011 / t9 + xi(1)).lt.-58.) THEN   !Overabundance of neutrinos.
            y(1) = 1.                !Essentially all neutrons.
            y(2) = 1.e-25            !Very little of protons.
        ELSE
            y(1) = 1. / (ex(15.011 / t9 + xi(1)) + 1.)  !Initial n abundance (Ref 3).
            y(2) = 1. / (ex(-15.011 / t9 - xi(1)) + 1.) !Initial p abundance (Ref 3).
        END IF
    END IF
    IF (xi(1).ne.0.) THEN        !Electron neutrino degeneracy.
        cnorm = 1.
        tnu = .00001             !Low temperature.
        CALL rate1(0.00001)        !Find normalization constant at low temp.
        cnorm = 1 / tau / f(1)
    END IF
    y0(1) = y(1)
    y0(2) = y(2)

    !40--------FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED--------------------

    z = 5.930 / t9            !Inverse of temperature.
    CALL bessel(z)
    hv = 3.3683e+4 * eta1 * 2.75 !(Ref 4 but with final eta).
    phie = hv * (1.784e-5 * y(2))  !Chemical potential of electron (Ref 5).&
    /(.5*z**3*(bl1-2.*bl2+3.*bl3-4.*bl4+5.*bl5))
    rhob0 = hv * t9**3            !Baryon density.
    IF ((xi(1).eq.0.).and.(xi(2).eq.0.).and.(xi(3).eq.0)) THEN  !Nondegen.
        rhone0 = 7.366 * t9**4       !Electron neutrino density (Ref 6).
    END IF

    !50--------SET ABUNDANCES FOR REST OF NUCLIDES----------------------------------

    y(3) = y(1) * y(2) * rhob0 * ex(25.82 / t9) / (.471e+10 * t9**1.5)  !(Ref 7).
    y0(3) = y(3)
    DO i = 4, isize
        y(i) = ytmin              !Set rest to minimum abundance.
        y0(i) = y(i)               !Init abundances at beginning of iteration.
    END DO
    CALL rate0                   !Compute weak decay rates.
    RETURN

    !----------REFERENCES-----------------------------------------------------------
    !     1) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
    !          page 44, equation A15.
    !     2) Coulomb correction obtained by dividing by correction factor Fp(t9)
    !               Fp(t9) = 1 - 0.5(pi/(137<v>/c))
    !          Wagoner, R.V. 1973, Ap. J. 179, page 358.
    !     3) For the nondegenerate case:
    !          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
    !          page 4, equation 3.
    !        For the case with neutrino degeneracy:
    !          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49,
    !          page 417, equation 9.
    !     4) Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 250, equation 4.
    !          3.3683e+4 = Mu(ng/t9**3) with Mu the atomic mass, ng the
    !          photon density.  2.75 is for the 11/4 factor difference
    !          between the initial and final values of eta.
    !     5) Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
    !          Kellogg Radiation Lab preprint OAP-714.
    !          equation D.2.
    !     6) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
    !          page 43, equation A4.
    !          7.366 is used instead of 14.73 as the latter is the sum total
    !          for 2 neutrino species.
    !     7) Initial deuterium abundance from nuclear statistical equilibrium
    !          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
    !          page 19, equation 17.

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE derivs(loop)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] driver
    !     CALLS     - [subroutine] therm, rate1, rate4, rate3, rate2, sol

    !----------REMARKS.
    !     Computes derivatives of
    !       - Temperature
    !       - hv
    !       - Chemical potential
    !       - abundances

    !----------PARAMETERS.
    PARAMETER (nvar = 29)          !Number of variables to be evolved.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.
    PARAMETER (pi = 3.141593)

    !----------COMMON AREAS.
    COMMON /evolp1/ t9, hv, phie, y                   !Evolution parameters.
    COMMON /evolp2/ dt9, dhv, dphie, dydt             !Evolution parameters.
    COMMON /evolp3/ t90, hv0, phie0, y0               !Evolution parameters.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi(3)     !Model parameters.
    COMMON /ttime/   t, dt, dlt9dt                   !Time variables.
    COMMON /thermcb/  thm, hubcst                   !Dynamic variables.
    COMMON /endens/ rhone0, rhob0, rhob, rnb          !Energy densities.
    COMMON /nucdat/ am(nnuc), zm, dm                 !Nuclide data.
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags,counters.
    COMMON /runopt/ irun, isize, jsize               !Run options.


    !==========================DECLARATION DIVISION=================================

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature (in units of 10**9 K).
    REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
    REAL    phie                 !Chemical potential for electron.
    REAL    y(nnuc)              !Relative number abundances.

    !----------EVOLUTION PARAMETERS (DERIVATIVES).
    REAL    dt9                  !Change in temperature.
    REAL    dhv                  !Change in hv.
    REAL    dphie                !Change in chemical potential.
    REAL    dydt(nnuc)           !Change in rel number abundances.

    !----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
    REAL    y0(nnuc)             !Rel # abund at beginning of iteration.

    !----------MODEL PARAMETERS.
    REAL    g                    !Gravitational constant.
    REAL    cosmo                !Cosmological constant.

    !----------TIME VARIABLES.
    REAL    dlt9dt               !(1/t9)*d(t9)/d(t).

    !----------DYNAMIC VARIABLES.
    REAL    thm(14)              !Thermodynamic variables.
    REAL    hubcst               !Expansion rate.

    !----------ENERGY DENSITIES.
    REAL    rhob0                !Initial baryon mass density.
    REAL    rhob                 !Baryon mass density.
    REAL    rnb                  !Baryon mass density (ratio to init value).

    !----------NUCLIDE DATA.
    REAL    zm(nnuc)             !Charge of nuclide.
    REAL    dm(nnuc)             !Mass excess of nuclide.

    !----------COUNTERS AND FLAGS.
    INTEGER mbad                 !Indicates if gaussian elimination fails.

    !----------RUN OPTION.
    INTEGER irun                 !Run network size.
    INTEGER isize                !Number of nuclides in computation.

    !----------SUMS.
    REAL    sumy                 !Sum of abundances.
    REAL    sumzy                !Sum of charge*abundances.
    REAL    sumdy                !Sum of abundance flows.
    REAL    summdy               !Sum of (mass excess)*(abundance flows).
    REAL    sumzdy               !Sum of (charge)*(abundance flows).

    !----------DERIVATIVES.
    REAL    dphdt9               !d(phi e)/d(t9).
    REAL    dphdln               !d(phi e)/d(h).
    REAL    dphdzy               !d(phi e)/d(sumzy).
    REAL    dlndt9               !(1/h)*d(h)/d(t9).
    REAL    bar                  !Baryon density and pressure terms.

    !----------LOCAL VARIABLES.
    INTEGER loop                 !Counts which Runge-Kutta loop.


    !===========================PROCEDURE DIVISION==================================

    !10--------COMPUTE DERIVATIVES FOR ABUNDANCES-----------------------------------

    rnb = hv * t9 * t9 * t9 / rhob0   !Baryon mass density (ratio to init value).
    !..........VARIOUS THERMODYNAMIC QUANTITIES.
    CALL therm
    hubcst = sqrt((8. / 3.) * pi * g * (thm(10)) + (cosmo / 3.))  !Expansion rate.
    rhob = thm(9)             !Baryon mass density.
    !..........COMPUTE REACTION RATE COEFFICIENTS.
    CALL rate1(t9)
    GO TO (100, 110, 120), irun    !Run network selection.
    100  CONTINUE
    CALL rate4                 !Forward rate for all of reactions.
    110  CONTINUE
    CALL rate3                 !Forward rate for reactions with A < 19.
    120  CONTINUE
    CALL rate2                 !Forward rate for reactions with A < 10.
    !..........SOLVE COUPLED DIFFERENTIAL EQUATIONS.
    CALL sol(loop)
    IF (mbad.gt.0) RETURN        !Abort in case matrix not invertible.

    !20--------COMPUTE DERIVATIVES FOR TEMPERATURE, hv, AND CHEMICAL POTENTIAL------

    !..........INITIALIZE SUMS TO ZERO.
    sumy = 0.
    sumzy = 0.
    sumdy = 0.
    summdy = 0.
    sumzdy = 0.
    !..........ACCUMULATE TO GET SUM.
    DO i = 1, isize
        sumy = sumy + y(i)           !Sum of abundance.
        sumzy = sumzy + zm(i) * y(i)     !Sum of charge*abundance.
        sumdy = sumdy + dydt(i)        !Sum of abundance flow.
        summdy = summdy + dm(i) * dydt(i)  !Sum of (mass excess)*(abundance flow).
        sumzdy = sumzdy + zm(i) * dydt(i)  !Sum of (charge)*(abundance flow).
    END DO
    !..........CHANGES IN TEMPERATURE, hv, AND CHEMICAL POTENTIAL.
    dphdt9 = thm(12) * (-1.070e-4 * hv * sumzy / t9 - thm(11))
    dphdln = -thm(12) * 3.568e-5 * hv * sumzy
    dphdzy = thm(12) * 3.568e-5 * hv
    bar = 9.25e-5 * t9 * sumy + 1.388e-4 * t9 * sumdy / (3. * hubcst)&
            + summdy / (3. * hubcst)
    dlndt9 = -(thm(2) + thm(5) + thm(6) * dphdt9 + thm(9) * 1.388e-4 * &
            sumy) / (thm(1) + thm(3) + thm(4) + thm(7) + thm(9) * bar&
            + thm(6) * (dphdln + dphdzy * sumzdy / (3. * hubcst)))   !(Ref 1).
    dt9 = (3. * hubcst) / dlndt9
    dlt9dt = dt9 / t9
    dhv = -hv * ((3. * hubcst) + 3. * dlt9dt)                    !(Ref 2).
    dphie = dphdt9 * dt9 + dphdln * (3. * hubcst) + dphdzy * sumzdy  !(Ref 3).

    RETURN

    !----------REFERENCES-----------------------------------------------------------
    !     1)  Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
    !          Kellogg Radiation Lab preprint OAP-714,
    !          equation D.35.
    !     2)  Kawano, L., 1992, preprint, equation D.19.
    !     3)  Kawano, L., 1992, preprint, equation D.20.

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE accum

    !----------LINKAGES.
    !     CALLED BY - [subroutine] driver
    !     CALLS     - none

    !----------REMARKS.
    !     Output accumulator.

    !----------PARAMETERS.
    PARAMETER (nvar = 29)          !Number of variables to be evolved.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.
    PARAMETER (itmax = 40)         !Maximum # of lines to be printed.

    !----------COMMON AREAS.
    COMMON /evolp1/ t9, hv, phie, y                   !Evolution parameters.
    COMMON /compr/  cy, ct, t9i, t9f, ytmin, inc        !Computation parameters.
    COMMON /ttime/   t, dt, dlt9dt                   !Time variables.
    COMMON /thermcb/  thm, hubcst                   !Dynamic variables.
    COMMON /nucdat/ am, zm(nnuc), dm(nnuc)           !Nuclide data.
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags,counters.
    COMMON /outdat/ xout, thmout, t9out, tout, dtout, !Output data.&
    etaout,hubout
    COMMON /runopt/ irun, isize, jsize               !Run options.


    !==========================DECLARATION DIVISION=================================

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature (in units of 10**9 K).
    REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
    REAL    phie                 !Chemical potential for electron.
    REAL    y(nnuc)              !Relative number abundances.

    !----------COMPUTATION PARAMETERS.
    INTEGER inc                  !Accumulation increment.

    !----------TIME PARAMETERS.
    REAL    t                    !Time.
    REAL    dt                   !Time step.

    !----------DYNAMIC VARIABLES.
    REAL    thm(14)              !Thermodynamic variables.
    REAL    hubcst               !Expansion rate.

    !----------NUCLIDE DATA.
    REAL    am(nnuc)             !Atomic number of nuclide.

    !----------COUNTERS AND FLAGS.
    INTEGER ltime                !Indicates if output buffer printed.
    INTEGER it                   !# times accumulated in output buffer.
    INTEGER ip                   !# time steps after outputting a line.

    !----------OUTPUT ARRAYS.
    REAL    xout(itmax, nnuc)     !Nuclide mass fractions.
    REAL    thmout(itmax, 6)      !Thermodynamic variables.
    REAL    t9out(itmax)         !Temperature (in units of 10**9 K).
    REAL    tout(itmax)          !Time.
    REAL    dtout(itmax)         !Time step.
    REAL    etaout(itmax)        !Baryon-to-photon ratio.
    REAL    hubout(itmax)        !Expansion rate.

    !----------RUN OPTION.
    INTEGER isize                !Number of nuclides in computation.


    !===========================PROCEDURE DIVISION==================================

    it = it + 1                  !Set up accumulation counter.

    !10--------SET UP OUTPUT VARIABLES----------------------------------------------

    !..........DIVIDE NUMBER FRACTION BY THAT OF PROTON.
    DO i = 1, isize
        xout(it, i) = y(i) / y(2)
    END DO
    xout(it, 2) = y(2) * am(2)      !Exception for proton.
    xout(it, 6) = y(6) * am(6)      !Exception for helium.
    !..........SUM UP ABUNDANCES OF HEAVY NUCLIDES.
    xout(it, 10) = xout(it, 10) + xout(it, 11) + xout(it, 12) + xout(it, 13)&
            + xout(it, 14) + xout(it, 15) + xout(it, 16) + xout(it, 17)&
            + xout(it, 18) + xout(it, 19) + xout(it, 20) + xout(it, 21)&
            + xout(it, 22) + xout(it, 23) + xout(it, 24) + xout(it, 25)&
            + xout(it, 26)   !Li8 to O16.
    !..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC.
    t9out(it) = t9            !Temperature.
    tout(it) = t             !Time.
    thmout(it, 1) = thm(1)        !rho photon.
    thmout(it, 2) = thm(4)        !rho electron.
    thmout(it, 3) = thm(8)        !rho neutrino.
    thmout(it, 4) = thm(9)        !rho baryon.
    thmout(it, 5) = phie          !Chemical potential.
    thmout(it, 6) = thm(10)       !rho total.
    dtout(it) = dt            !Time step.
    etaout(it) = hv / (3.3683e+4)!Baryon to photon ratio.
    hubout(it) = hubcst        !Expansion rate.

    !20--------INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------------

    IF ((it.eq.itmax).or.(ip.lt.inc)) ltime = 1
    RETURN

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE therm

    !----------LINKAGES.
    !     CALLED BY - [subroutine] derivs
    !     CALLS     - [subroutine] bessel, nudens
    !               - [function] ex

    !----------REMARKS.
    !     Computes various temperature dependent thermodynamic quantities.

    !----------PARAMETER.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.
    PARAMETER (q = 2.531)          !(mass(neutron)-mass(proton))/m(electron)

    !----------COMMON AREAS.
    COMMON /evolp1/ t9, hv, phie, y(nnuc)             !Evolution parameters.
    COMMON /compr/  cy, ct, t9i, t9f, ytmin, inc        !Computation parameters.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /thermcb/  thm, hubcst                   !Dynamic variables.
    COMMON /endens/ rhone0, rhob0, rhob, rnb          !Energy densities.
    COMMON /besselcb/ bl1, bl2, bl3, bl4, bl5, !Eval of function bl(z).&
    bm1,bm2, bm3, bm4, bm5, !Eval of function bm(z).&
    bn1,bn2, bn3, bn4, bn5            !Eval of function bn(z).
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.


    !==========================DECLARATION DIVISION=================================

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature (in units of 10**9 K).
    REAL    phie                 !Chemical potential for electron.

    !----------COMPUTATION PARAMETERS.
    REAL    t9i                  !Initial temperature (in 10**9 K).

    !----------EARLY UNIVERSE MODEL PARAMETERS.
    REAL    xnu                  !Number of neutrino species.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------DYNAMIC VARIABLES.
    REAL    thm(14)              !Thermodynamic variables.

    !----------ENERGY DENSITIES.
    REAL    rhone0               !Initial electron neutrino mass density.
    REAL    rhob0                !Initial baryon mass density.
    REAL    rnb                  !Baryon mass density (ratio to init value).

    !----------EVALUATION OF FUNCTIONS bl,bm,bn.
    REAL    bl1, bl2, bl3, bl4, bl5  !Evaluation of function bl(z).
    REAL    bm1, bm2, bm3, bm4, bm5  !Evaluation of function bm(z).
    REAL    bn1, bn2, bn3, bn4, bn5  !Evaluation of function bn(z).

    !----------NEUTRINO PARAMETERS.
    REAL    tnu                  !Neutrino temperature.
    REAL    rhonu                !Neutrino energy density.
    INTEGER nu                   !Type of neutrino.

    !----------LOCAL VARIABLE.
    REAL    z                    !Defined by z = m(electron)*c**2/k*t9.


    !===========================PROCEDURE DIVISION==================================

    !10--------COMPUTE FACTORS------------------------------------------------------

    z = 5.930 / t9                 !z = m(electron)c**2/k(t9).
    tnu = ((rnb)**(1. / 3.)) * t9i   !Neutrino temperature.
    !..........FACTORS OF z.
    z1 = z
    z2 = z * z
    z3 = z * z * z
    z4 = z * z * z * z
    z5 = z * z * z * z * z
    !..........TRIGNOMETRIC FUNCTION VALUES.
    IF (phie.le.17.) THEN        !No chance of overflow.
        cosh1 = cosh(phie)
        cosh2 = cosh(2. * phie)
        cosh3 = cosh(3. * phie)
        cosh4 = cosh(4. * phie)
        cosh5 = cosh(5. * phie)
        sinh1 = sinh(phie)
        sinh2 = sinh(2. * phie)
        sinh3 = sinh(3. * phie)
        sinh4 = sinh(4. * phie)
        sinh5 = sinh(5. * phie)
    ELSE
        cosh1 = 0.
        cosh2 = 0.
        cosh3 = 0.
        cosh4 = 0.
        cosh5 = 0.
        sinh1 = 0.
        sinh2 = 0.
        sinh3 = 0.
        sinh4 = 0.
        sinh5 = 0.
    END IF
    CALL bessel(z)

    !20--------COMPUTE THERMODYNAMIC VARIABLES--------------------------------------

    thm(1) = 8.418 * t9 * t9 * t9 * t9                               !(Ref 1).
    thm(2) = 4. * thm(1) / t9                                    !(Ref 2).
    thm(3) = thm(1) / 3.                                       !(Ref 3).
    thm(4) = 3206. * (bm1 * cosh1 - bm2 * cosh2 + bm3 * cosh3        !(Ref 4).&
    - bm4*cosh4 + bm5*cosh5)
    thm(5) = 3206. * (z / t9) * (bn1 * cosh1 - 2. * bn2 * cosh2          !(Ref 5).&
    + 3.*bn3*cosh3 - 4.*bn4*cosh4 + 5.*bn5*cosh5)
    thm(6) = 3206. * (bm1 * sinh1 - 2. * bm2 * sinh2 + 3. * bm3 * sinh3  !(Ref 6).&
    - 4.*bm4*sinh4 + 5.*bm5*sinh5)
    thm(7) = 3206. * (bl1 * cosh1 / z - bl2 * cosh2 / (2. * z)           !(Ref 7).&
    + bl3*cosh3/(3.*z) - bl4*cosh4/(4.*z)&
            + bl5*cosh5/(5.*z))
    IF ((xi(1).eq.0.).and.(xi(2).eq.0.).and.(xi(3).eq.0)) THEN  !Nondegen.
        thm(8) = xnu * rhone0 * (rnb**(4. / 3.))                      !(Ref 8).
    ELSE                         !Include effects of neutrino degeneracy.
        thm(8) = 0.
        DO nu = 1, xnu              !For every neutrino family.
            CALL nudens              !Compute neutrino energy density.
            thm(8) = thm(8) + 12.79264 * rhonu  !Have 12.79264 from units change.
        END DO
    END IF
    thm(9) = rhob0 * rnb                                       !(Ref 9).
    thm(10) = thm(1) + thm(4) + thm(8) + thm(9)               !(Ref 10).
    thm(11) = -(z**3 / t9) * (sinh1 * (3. * bl1 - z * bm1) - sinh2 * (3. * bl2  !(Ref 11).&
    -2.*z*bm2) + sinh3*(3.*bl3-3.*z*bm3) - sinh4&
            *(3.*bl4-4.*z*bm4) + sinh5*(3.*bl5-5.*z*bm5))
    thm(12) = z**3 * (cosh1 * bl1 - 2. * cosh2 * bl2                   !(Ref 12).&
    + 3.*cosh3*bl3 - 4.*cosh4*bl4 + 5.*cosh5*bl5)
    IF (thm(12).ne.0.) thm(12) = 1. / thm(12)
    thm(13) = 1.000 + 0.565 / z1 - 6.382 / z2 + 11.108 / z3         !(Ref 13).&
    + 36.492/z4 + 27.512/z5
    thm(14) = (5.252 / z1 - 16.229 / z2 + 18.059 / z3 + 34.181 / z4   !(Ref 14).&
    + 27.617/z5)*ex(-q*z)

    RETURN

    !----------REFERENCES AND NOTES-------------------------------------------------
    !     1)  thm(1)  = rho photon
    !         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
    !          page 43, equation A2.)
    !     2)  thm(2)  = d(rho photon)/d(t9)
    !     3)  thm(3)  = (p photon)/c**2
    !         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967,
    !          page 43, equation A3.)
    !     4)  thm(4)  = rho electron+positron
    !         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
    !          page 281, equation B44.)
    !     5)  thm(5)  = d(rho electron+positron)/d(t9)
    !     6)  thm(6)  = d(rho electron+positron)/d(phi e)
    !     7)  thm(7)  = (p electron+positron)/c**2
    !         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
    !          page 279, equation B27.)
    !     8)  thm(8)  = rho neutrino
    !                 = # neutrino species x rho electron neutrino (nondegenerate)
    !                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate)
    !     9)  thm(9)  = rho baryon
    !     10) thm(10) = rho total
    !                 = rho photon + rho electron+positron + rho neutrino
    !                              + rho baryon
    !     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\
    !                   d(t9) \  2  (mc**2)**3                 /
    !     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\
    !                   d(phi e) \  2  (mc**2)**3                 /
    !     13) thm(13) = rate for n->p
    !     14) thm(14) = rate for p->n

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE bessel(z)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] start, therm
    !     CALLS     - [subroutine] knux

    !----------REMARKS.
    !     Evaluates functions bl(z), bm(z), and bn(z) using solutions to
    !     modified Bessel functions.

    !----------COMMON AREAS.
    COMMON /besselcb/ bl1, bl2, bl3, bl4, bl5, !Eval function bl(z).&
    bm1,bm2, bm3, bm4, bm5, !Eval function bm(z).&
    bn1,bn2, bn3, bn4, bn5            !Eval function bn(z).
    COMMON /kays/   bk0, bk1, bk2, bk3, bk4            !Coefficients K.


    !==========================DECLARATION DIVISION=================================

    !----------EVALUATION OF FUNCTIONS bl,bm,bn.
    REAL    bl1, bl2, bl3, bl4, bl5  !Single variables equivalenced to array blz.
    REAL    bm1, bm2, bm3, bm4, bm5  !Single variables equivalenced to array bmz.
    REAL    bn1, bn2, bn3, bn4, bn5  !Single variables equivalenced to array bnz.

    !----------EVALUATIION OF MODIFIED BESSEL FUNCTIONS.
    REAL    bk0, bk1, bk2, bk3, bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).

    !----------LOCAL VARIABLES.
    REAL    blz(5)               !Array containing values from function bl.
    REAL    bmz(5)               !Array containing values from function bm.
    REAL    bnz(5)               !Array containing values from function bn.
    REAL    z                    !Defined by z = m(electron)*c**2/k*t9.
    REAL    r                    !Multiples of z.

    !----------EQUIVALENCE STATEMENTS.
    EQUIVALENCE (blz(1), bl1), (blz(2), bl2), (blz(3), bl3), (blz(4), bl4), &
            (blz(5), bl5)
    EQUIVALENCE (bmz(1), bm1), (bmz(2), bm2), (bmz(3), bm3), (bmz(4), bm4), &
            (bmz(5), bm5)
    EQUIVALENCE (bnz(1), bn1), (bnz(2), bn2), (bnz(3), bn3), (bnz(4), bn4), &
            (bnz(5), bn5)


    !===========================PROCEDURE DIVISION==================================

    !10--------LOCALLY DEFINED FUNCTIONS--------------------------------------------

    bl(z) = bk2 / z                !Function bl.
    bm(z) = 0.25 * (3. * bk3 + bk1) / z  !Function bm.
    bn(z) = 0.5 * (bk4 + bk2) / z      !Function bn.

    !20--------CALCULATE FOR 1 THRU 5 Z---------------------------------------------

    DO i = 1, 5
        r = i * z                      !Multiples of z.
        CALL knux(r)               !Get k0(r),k1(r),k2(r),k3(r),k4(r),k(5).
        blz(i) = bl(r)             !Put value from function bl into array.
        bmz(i) = bm(r)             !Put value from function bm into array.
        bnz(i) = bn(r)             !Put value from function bn into array.
    END DO
    RETURN

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE knux(z)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] bessel
    !     CALLS     - [function] exp

    !----------REMARKS.
    !     A subroutine for modified bessel functions of the second kind
    !     k-nu(z).

    !----------COMMON AREAS.
    COMMON /kays/   bk0, bk1, bk2, bk3, bk4            !Coefficients K.


    !===========================DECLARATION DIVISION================================

    !-----------MODIFIED BESSEL FUNCTION VALUES.
    REAL    bk0, bk1              !Values k0(z),k1(z)
    REAL    bi0, bi1              !Values i0(z),i1(z).
    REAL    bk2, bk3, bk4          !Values k2(z),k3(z),k4(z).

    !-----------EXPANSION COEFFICIENTS.
    REAL    ci0(7)               !Expansion coefficients for i0 (z.le.2).
    REAL    ci1(7)               !Expansion coefficients for i1 (z.le.2).
    REAL    ck0(7)               !Expansion coefficients for k0 (z.le.2).
    REAL    ck1(7)               !Expansion coefficients for k1 (z.le.2).
    REAL    c0(7)                !Expansion coefficients for k0 (z.gt.2).
    REAL    c1(7)                !Expansion coefficients for k1 (z.gt.2).

    !-----------VARIABLES TO BE EVALUATED.
    REAL    z                    !Input variable.
    REAL    y                    !Expansion variable = z/2.
    REAL    t                    !Expansion variable = z/3.75.
    REAL    coeff                !Logrithmic or exponential coefficient.


    !==============================DATA DIVISION====================================

    !----------EXPANSION COEFFICIENTS.
    DATA ci0 / 1., &
            3.5156229, 3.0899424, 1.2067492, &
            0.2659732, 0.0360768, 0.0045813/
    DATA ci1 / 0.5, &
            0.87890594, 0.51498869, 0.15084934, &
            0.02658733, 0.00301532, 0.00032411/
    DATA ck0 /-0.57721566, &
            0.42278420, 0.23069756, 0.03488590, &
            0.00262698, 0.00010750, 0.00000740/
    DATA ck1 / 1., &
            0.15443144, -0.67278579, -0.18156897, &
            -0.01919402, -0.00110404, -0.00004686/
    DATA c0  / 1.25331414, &
            -0.07832358, 0.02189568, -0.01062446, &
            0.00587872, -0.00251540, 0.00053208/
    DATA c1  / 1.25331414, &
            0.23498619, -0.03655620, 0.01504268, &
            -0.00780353, 0.00325614, -0.00068245/


    !===========================PROCEDURE DIVISION==================================

    !10--------COMPUTE K0 AND K1----------------------------------------------------

    IF (z.le.2.) THEN            !(Ref. 1).
        !..........COMPUTE FACTORS.
        t = (z / 3.75)
        y = (z / 2)
        coeff = alog(y)
        !..........VALUES FOR i0(z) and i1(z).
        bi0 = ci0(1)
        bi1 = ci1(1)
        bk0 = ck0(1)
        bk1 = ck1(1)
        DO i = 2, 7
            bi0 = bi0 + ci0(i) * t**(2 * (i - 1))
            bi1 = bi1 + ci1(i) * t**(2 * (i - 1))
            bk0 = bk0 + ck0(i) * y**(2 * (i - 1))
            bk1 = bk1 + ck1(i) * y**(2 * (i - 1))
        END DO
        !..........VALUES FOR k0(z) and k1(z).
        bk0 = -coeff * bi0 + bk0
        bk1 = coeff * bi1 * z + bk1 / z
    ELSE !(z.le.2.)               !(Ref. 2).
        !..........COMPUTE FACTORS.
        y = (2.0 / z)
        coeff = (ex(-z) / sqrt(z))
        !..........VALUES FOR k0(z) and k1(z).
        bk0 = c0(1)
        bk1 = c1(1)
        DO i = 2, 7
            bk0 = bk0 + c0(i) * y**(i - 1)
            bk1 = bk1 + c1(i) * y**(i - 1)
        END DO
        bk0 = coeff * bk0
        bk1 = coeff * bk1
    END IF !(z.le.2.)

    !20--------FIND K2, K3, AND K4 BY ITERATION (Ref. 3)----------------------------

    bk2 = 2. * (bk1 / z) + bk0       !k2(z).
    bk3 = 4. * (bk2 / z) + bk1       !k3(z).
    bk4 = 6. * (bk3 / z) + bk2       !k4(z).

    RETURN

    !----------REFERENCES-----------------------------------------------------------
    !     Handbook of Mathematical Functions (Abramowitz and Stegun),
    !       Dover Publications, Inc., New York
    !       1) Polynomial approximations for z.le.2
    !         page 378, equations 9.8.1 and 9.8.3.
    !         page 379, equations 9.8.5 and 9.8.7.
    !       2) Polynomial approximations for z > 2
    !         page 379, equations 9.8.6 and 9.8.8.
    !       3) Recursion relation from 1st line of 9.6.26, page 376.

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE nudens

    !----------LINKAGES.
    !     CALLED BY - [subroutine] therm
    !     CALLS     - [function] xintd, eval

    !----------REMARKS.
    !     Computes energy density contribution from neutrinos.

    !----------PARAMTER.
    PARAMETER (iter = 50)          !Number of gaussian quads.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !----------EXTERNAL FUNCTIONS.
    EXTERNAL func5               !Integral for neutrinos.
    EXTERNAL func6               !Integral for antineutrinos.


    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    rhonu                !Neutrino energy density.
    INTEGER nu                   !Which neutrino type.

    !----------LOCAL VARIABLES.
    REAL    uplim1               !Upper limit for neutrino energy integral.
    REAL    uplim2               !Upper limit for antineu energy integral.


    !===========================PROCEDURE DIVISION==================================

    !10--------COMPUTE NEUTRINO ENERGY DENSITIES------------------------------------

    IF (abs(xi(nu)).le.0.03) THEN
        !..........SMALL xi APPROXIMATION.
        rhonu = 2. * (3.14159**2 / 30.) * (tnu)**4&
                * (7. / 8. + (15. / (4 * 3.14159**2)) * xi(nu)**2&
                        + (15. / (8. * 3.14159**4)) * xi(nu)**4)
    ELSE
        IF (abs(xi(nu)).ge.30.) THEN
            !..........LARGE xi APPROXIMATION.
            rhonu = ((tnu)**4) / (8. * 3.14159**2) * xi(nu)**4&
                    * (1 + 12. * 1.645 / xi(nu)**2)
        ELSE
            !..........DO INTEGRATION
            uplim1 = (88.029 + xi(nu)) * tnu
            uplim2 = (88.029 - xi(nu)) * tnu
            IF (uplim2.le.0.) THEN
                rhonu = xintd(0., uplim1, func5, iter)
            ELSE
                rhonu = xintd(0., uplim1, func5, iter)&
                        + xintd(0., uplim2, func6, iter)
            END IF
        END IF !(abs(xi(nu)).ge.30.)
    END IF !(abs(xi(nu)).le.0.03)
    RETURN

    !----------REFERENCES-----------------------------------------------------------
    !     Forms of the integrals involved can be found in
    !       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
    !       Freese, K., Kolb, E.W., Turner, M.S., 1983, Phys. Rev. D, 27, 1689.

END


!========================IDENTIFICATION DIVISION================================


!===========================PROCEDURE DIVISION==================================

!10--------1ST PART OF INTEGRAL FOR n->p RATE-----------------------------------
!***************************************************************************

real function func1(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - [function] ex

    !----------REMARKS.
    !     Contains integrands to be integrated.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    cnorm                !Normalizing constant.

    !----------LOCAL VARIABLES.
    REAL    x                    !Value at which function is evaluated.
    REAL    part1                !Exponential expression with photon temp.
    REAL    part2                !Exponential expression with neutrino temp.

    IF (x.le.0.) THEN
        func1 = 0.
    ELSE
        part1 = 1. / (1. + ex(-.511 * x / t9mev))
        part2 = 1. / (1. + ex(+(x - 2.531) * (.511 / tnmev) - xi(1)))
        func1 = cnorm * x * (x - 2.531)**2 * (x**2 - 1)**.5 * part1 * part2
    END IF
    RETURN
end

!20--------2ND PART OF INTEGRAL FOR n->p RATE-----------------------------------

real function func2(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - [function] ex

    !----------REMARKS.
    !     Contains integrands to be integrated.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    cnorm                !Normalizing constant.

    !----------LOCAL VARIABLES.
    REAL    x                    !Value at which function is evaluated.
    REAL    part1                !Exponential expression with photon temp.
    REAL    part2                !Exponential expression with neutrino temp.

    IF (x.le.1.) THEN
        func2 = 0.
    ELSE
        part1 = 1. / (1. + ex(+.511 * x / t9mev))
        part2 = 1. / (1. + ex(-(x + 2.531) * (.511 / tnmev) - xi(1)))
        func2 = cnorm * x * (x + 2.531)**2 * (x**2 - 1)**.5 * part1 * part2
    END IF
    RETURN
end

!30--------1ST PART OF INTEGRAL FOR p->n RATE-----------------------------------

real function func3(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - [function] ex

    !----------REMARKS.
    !     Contains integrands to be integrated.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    cnorm                !Normalizing constant.

    !----------LOCAL VARIABLES.
    REAL    x                    !Value at which function is evaluated.
    REAL    part1                !Exponential expression with photon temp.
    REAL    part2                !Exponential expression with neutrino temp.

    IF (x.le.1.) THEN
        func3 = 0.
    ELSE
        part1 = 1. / (1. + ex(-.511 * x / t9mev))
        part2 = 1. / (1. + ex(+(x + 2.531) * (.511 / tnmev) + xi(1)))
        func3 = cnorm * x * (x + 2.531)**2 * (x**2 - 1)**.5 * part1 * part2
    END IF
    RETURN
end

!40--------2ND PART OF INTEGRAL FOR p->n RATE-----------------------------------

real function  func4(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - [function] ex

    !----------REMARKS.
    !     Contains integrands to be integrated.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    cnorm                !Normalizing constant.

    !----------LOCAL VARIABLES.
    REAL    x                    !Value at which function is evaluated.
    REAL    part1                !Exponential expression with photon temp.
    REAL    part2                !Exponential expression with neutrino temp.

    IF (x.le.1.) THEN
        func4 = 0.
    ELSE
        part1 = 1. / (1. + ex(+.511 * x / t9mev))
        part2 = 1. / (1. + ex(-(x - 2.531) * (.511 / tnmev) + xi(1)))
        func4 = cnorm * x * (x - 2.531)**2 * (x**2 - 1)**.5 * part1 * part2
    END IF
    RETURN
end

!50--------INTEGRAL FOR ENERGY DENSITY OF NEUTRINO------------------------------

real function func5(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - [function] ex

    !----------REMARKS.
    !     Contains integrands to be integrated.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    cnorm                !Normalizing constant.

    !----------LOCAL VARIABLES.
    REAL    x                    !Value at which function is evaluated.
    REAL    part1                !Exponential expression with photon temp.
    REAL    part2                !Exponential expression with neutrino temp.

    func5 = 1. / (2 * 3.14159**2) * x**3 / (1. + exp(x / tnu - xi(nu)))
    RETURN
end

!60--------INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO--------------------------

real function func6(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - [function] ex

    !----------REMARKS.
    !     Contains integrands to be integrated.

    !----------COMMON AREAS.
    COMMON /modpr/  g, tau, xnu, c(3), cosmo, xi        !Model parameters.
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Integration parameters.

    !==========================DECLARATION DIVISION=================================

    !----------MODEL PARAMETERS.
    REAL    xi(3)                !Neutrino degeneracy parameters.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature (units of 10**9 K).
    REAL    cnorm                !Normalizing constant.

    !----------LOCAL VARIABLES.
    REAL    x                    !Value at which function is evaluated.
    REAL    part1                !Exponential expression with photon temp.
    REAL    part2                !Exponential expression with neutrino temp.

    func6 = 1. / (2 * 3.14159**2) * x**3 / (1. + exp(x / tnu + xi(nu)))
    RETURN
end

!----------REFERENCES-----------------------------------------------------------
!     Forms of the integrals involved can be found in
!       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
!       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.




!========================IDENTIFICATION DIVISION================================

FUNCTION xintd (xlow, xhi, func, nq)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] rate1, nudens
    !     CALLS     - none

    !----------REMARKS.
    !     Computes the integral of the function "func".


    !==========================DECLARATION DIVISION=================================

    !----------INPUT VARIABLES.
    REAL    xlow                 !Array of low limits.
    REAL    xhi                  !Array of high limits.
    INTEGER nq                   !Number of six point gaussian quads.

    !----------COMPUTATION VARIABLES.
    REAL    dist                 !Size of quad interval.
    REAL    cent                 !Center of quad interval.
    REAL    x                    !Variables of integration.
    REAL    sum                  !Summation of terms.

    !----------COUNTERS.
    INTEGER nint                 !Interval number.
    INTEGER npnt                 !Point number.
    INTEGER np                   !Total number of points in interval.

    !----------ABSCISSAS AND WEIGHT FACTORS.
    REAL    u(6)                 !Abscissas.
    REAL    w(6)                 !Weight factor.


    !==============================DATA DIVISION====================================

    !----------ABSCISSAS AND WEIGHT FACTORS.
    DATA u/-.93246951420315, -.66120938646627, -.23861918608320, &
            .23861918608320, .66120938646627, .93246951420315/
    DATA w/.17132449237917, .36076157304814, .46791393457269, &
            .46791393457269, .36076157304814, .17132449237917/
    DATA np/6/              !6 point Gaussian integration.


    !===========================PROCEDURE DIVISION==================================

    !10--------DO INTEGRATION-------------------------------------------------------

    sum = 0.
    dist = (xhi - xlow) / float(nq) !Size of quad interval.
    DO nint = 1, nq
        cent = xlow + (float(nint) - 0.5) * dist  !Center of interval.
        DO npnt = 1, np
            x = cent + 0.5 * dist * u(npnt) !Integration point.
            f = func(x)            !Evaluate function x(1).
            sum = sum + f * w(npnt)      !Add up sum.
        END DO
    END DO

    !20--------GET INTEGRAL VALUE---------------------------------------------------

    xintd = sum * dist * 0.5         !Do integral.
    RETURN

END


!========================IDENTIFICATION DIVISION================================

FUNCTION ex(x)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol
    !               - [function] eval
    !     CALLS     - none

    !----------REMARKS.
    !     Exponential function with underflow precaution.


    !===========================PROCEDURE DIVISION==================================

    IF (x.gt.88.029) THEN        !In danger of overflow.
        ex = exp(88.029)
    ELSE
        IF (x.lt.-88.722) THEN     !In danger of underflow.
            ex = 0.
        ELSE                       !Value of x in allowable range.
            ex = exp(x)
        END IF
    END IF
    RETURN

    !----------NOTE-----------------------------------------------------------------
    !     The overflow limit for the VAX/VMS system is exp(88.029).
    !     The underflow limit for the VAX/VMS system is exp(-88.722).

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE sol(loop)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] derivs
    !     CALLS     - [subroutine] eqslin
    !               - [function] ex

    !----------REMARKS.
    !     Computes reverse strong and electromagnetic reaction rates.
    !     Fills and solves matrix equation for dydt(i).

    !----------PARAMETERS.
    PARAMETER (ir = 5)             !Input unit number.
    PARAMETER (iw = 6)             !Output unit number.
    PARAMETER (nrec = 88)          !Number of nuclear reactions.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.

    !-----------COMMON AREAS.
    COMMON /recpr/  iform, ii, jj, kk, ll, rev, q9       !Reaction parameters names.
    COMMON /rates/  f, r                            !Reaction rates.
    COMMON /evolp1/ t9, hv, phie, y                   !Evolution parameters.
    COMMON /evolp2/ dt9, dhv, dphie, dydt             !Evolution parameters.
    COMMON /evolp3/ t90, hv0, phie0, y0               !Evolution parameters.
    COMMON /ttime/   t, dt, dlt9dt                   !Time varying parameters.
    COMMON /thermcb/  thm(14), hubcst                 !Dynamic variables.
    COMMON /endens/ rhone0, rhob0, rhob, rnb          !Energy densities.
    COMMON /lncoef/ a, b, yx                         !Linear eqn coefficients.
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags,counters.
    COMMON /runopt/ irun, isize, jsize               !Run option.


    !==========================DECLARATION DIVISION=================================

    !----------REACTION PARAMETERS.
    INTEGER iform(nrec)          !Reaction code number (1-88).
    INTEGER ii(nrec)             !Incoming nuclide type (1-26).
    INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
    INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
    INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
    REAL    rev(nrec)            !Reverse reaction coefficient.
    REAL    q9(nrec)             !Energy released in reaction.

    !----------REACTION RATES.
    REAL    f(nrec)              !Forward reaction rate coefficients.
    REAL    r(nrec)              !Reverse reaction rate coefficients.

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature (in units of 10**9 K).
    REAL    y(nnuc)              !Relative number abundances.

    !----------EVOLUTION PARAMETERS (DERIVATIVES).
    REAL    dydt(nnuc)           !Change in rel number abundances.

    !----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
    REAL    y0(nnuc)             !Rel # abund at start of iteration.

    !----------TIME VARIABLES.
    REAL    dt                   !Time step.

    !----------DYNAMIC VARIABLES.
    REAL    hubcst               !Expansion rate.

    !----------ENERGY DENSITIES.
    REAL    rhob                 !Baryon mass density.

    !----------COMPONENTS OF MATRIX EQUATION.
    DOUBLE PRECISION a(nnuc, nnuc)!Relates y(t-dt) to y(t).
    REAL    b(nnuc)              !Contains y0 in inverse order.
    REAL    yx(nnuc)             !yy in reverse order.

    !----------COUNTERS AND FLAGS.
    INTEGER loop                 !Counts which Runge-Kutta loop.
    INTEGER ip                   !# time steps after outputting a line.
    INTEGER mbad                 !Indicates if gaussian elimination fails.

    !----------RUN OPTIONS.
    INTEGER isize                !Number of nuclides in computation.
    INTEGER isize1               !Equals isize + 1.
    INTEGER jsize                !Number of reactions in computation.

    !----------EVOLUTION EQUATION COEFFICIENTS.
    INTEGER i, j, k, l              !Equate to ii,jj,kk,ll.
    REAL    ri, rj, rk, rl          !Equate to si,sj,sk,sl.
    REAL    ci, cj, ck, cl          !Coefficients of rate equation.

    !----------LOCAL VARIABLES.
    REAL    yy(nnuc)             !Abundances at end of iteration.
    REAL    si(11), sj(11), sk(11), sl(11)  !# of nuclide i,j,k,l
    REAL    bdln                 !(10**(-5))*volume expansion rate.
    INTEGER ind                  !Equate to iform.
    INTEGER ierror               !Element which does not converge.


    !==============================DATA DIVISION====================================

    !----------NUMBER OF NUCLIDES IN REACTION TYPES 1-11.
    DATA si /1., 1., 1., 1., 1., 2., 3., 2., 1., 1., 2./
    DATA sj /0., 1., 1., 0., 1., 0., 0., 1., 1., 1., 0./
    DATA sk /0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 2./
    DATA sl /1., 1., 1., 2., 2., 1., 1., 1., 2., 3., 1./


    !===========================PROCEDURE DIVISION==================================

    !10--------TEMPERATURE FACTORS AND INITIAL VALUES-------------------------------

    !..........TEMPERATURE FACTORS.
    t932 = t9**1.5              !t9**(3/2).
    t9m32 = 1. / t932              !t9**(-3/2).
    !..........MATRIX SIZE.
    isize1 = isize + 1
    !..........INITIALIZE A-MATRIX.
    DO i = 1, isize
        DO j = 1, isize
            a(j, i) = 0.d0            !Set a-matrix to zero.
        END DO
    END DO

    !20--------COMPUTE FACTORS FOR THE A-MATRIX-------------------------------------

    DO n = 1, jsize
        !..........EQUATE VARIABLES TO ARRAYS.
        ind = iform(n)             !Type of reaction.
        i = ii(n)                  !ID # of incoming nuclide i.
        j = jj(n)                  !ID # of incoming nuclide j.
        k = kk(n)                  !ID # of outgoing nuclide k.
        l = ll(n)                  !ID # of outgoing nuclide l.
        IF ((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) THEN  !Reaction okay.
            ri = si(ind)             !# of incoming nuclide i.
            rj = sj(ind)             !# of incoming nuclide j.
            rk = sk(ind)             !# of outgoing nuclide k.
            rl = sl(ind)             !# of outgoing nuclide l.
            !..........COMPUTE DIFFERENT REACTION RATES.
            GO TO (201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211), ind
            201      CONTINUE                 !1-0-0-1 configuration.
            ci = f(n)              !(Ref 1).
            cj = 0.
            ck = 0.
            cl = r(n)
            GO TO 212
            202      CONTINUE                 !1-1-0-1 configuration.
            r(n) = rev(n) * 1.e+10 * t932 * ex(-q9(n) / t9) * f(n)  !(Ref 2).
            f(n) = rhob * f(n)
            ci = y(j) * f(n) / 2.
            cj = y(i) * f(n) / 2.
            ck = 0.
            cl = r(n)
            GO TO 212
            203      CONTINUE                 !1-1-1-1 configuration.
            f(n) = rhob * f(n)
            r(n) = rev(n) * ex(-q9(n) / t9) * f(n)  !(Ref 3).
            ci = y(j) * f(n) / 2.
            cj = y(i) * f(n) / 2.
            ck = y(l) * r(n) / 2.
            cl = y(k) * r(n) / 2.
            GO TO 212
            204      CONTINUE                 !1-0-0-2 configuration.
            ci = f(n)
            cj = 0.
            ck = 0.
            cl = y(l) * r(n) / 2.
            GO TO 212
            205      CONTINUE                 !1-1-0-2 configuration.
            f(n) = rhob * f(n)
            r(n) = rev(n) * ex(-q9(n) / t9) * f(n)  !(Ref 3).
            ci = y(j) * f(n) / 2.
            cj = y(i) * f(n) / 2.
            ck = 0.
            cl = y(l) * r(n) / 2.
            GO TO 212
            206      CONTINUE                 !2-0-1-1 configuration.
            f(n) = rhob * f(n)
            r(n) = rev(n) * ex(-q9(n) / t9) * f(n)  !(Ref 3).
            ci = y(i) * f(n) / 2.
            cj = 0.
            ck = y(l) * r(n) / 2.
            cl = y(k) * r(n) / 2.
            GO TO 212
            207      CONTINUE                 !3-0-0-1 configuration.
            r(n) = rev(n) * 1.e+20 * t932 * t932 * ex(-q9(n) / t9) * f(n)  !(Ref 4).
            f(n) = rhob * rhob * f(n)
            ci = y(i) * y(i) * f(n) / 6.
            cj = 0.
            ck = 0.
            cl = r(n)
            GO TO 212
            208      CONTINUE                 !2-1-0-1 configuration.
            r(n) = rev(n) * 1.e+20 * t932 * t932 * ex(-q9(n) / t9) * f(n)  !(Ref 4).
            f(n) = rhob * rhob * f(n)
            ci = y(j) * y(i) * f(n) / 3.
            cj = y(i) * y(i) * f(n) / 6.
            ck = 0.
            cl = r(n)
            GO TO 212
            209      CONTINUE                 !1-1-1-2 configuration.
            f(n) = rhob * f(n)
            r(n) = rev(n) * 1.e-10 * t9m32 * rhob * ex(-q9(n) / t9) * f(n)  !(Ref 5).
            ci = y(j) * f(n) / 2.
            cj = y(i) * f(n) / 2.
            ck = y(l) * y(l) * r(n) / 6.
            cl = y(k) * y(l) * r(n) / 3.
            GO TO 212
            210      CONTINUE                 !1-1-0-3 configuration.
            f(n) = rhob * f(n)
            r(n) = rev(n) * 1.e-10 * t9m32 * rhob * ex(-q9(n) / t9) * f(n)  !(Ref 5).
            ci = y(j) * f(n) / 2.
            cj = y(i) * f(n) / 2.
            ck = 0.
            cl = y(l) * y(l) * r(n) / 6.
            GO TO 212
            211      CONTINUE                 !2-0-2-1 configuration.
            f(n) = rhob * f(n)
            r(n) = rev(n) * 1.e-10 * t9m32 * rhob * ex(-q9(n) / t9) * f(n)  !(Ref 5).
            ci = y(i) * f(n) / 2.
            cj = 0.
            ck = y(l) * y(k) * r(n) / 3.
            cl = y(k) * y(k) * r(n) / 6.
            212      CONTINUE

            !30--------CONSTRUCT THE A-MATRIX-----------------------------------------------

            i = isize1 - i           !Invert i index.
            j = isize1 - j           !Invert j index.
            k = isize1 - k           !Invert k index.
            l = isize1 - l           !Invert l index.
            !..........FILL I NUCLIDE COLUMN.
            IF (j.le.isize) a(j, i) = a(j, i) + rj * ci
            IF (k.le.isize) a(k, i) = a(k, i) - rk * ci
            a(i, i) = a(i, i) + ri * ci
            a(l, i) = a(l, i) - rl * ci
            !..........FILL J NUCLIDE COLUMN.
            IF (j.le.isize) THEN
                a(j, j) = a(j, j) + rj * cj
                IF (k.le.isize) a(k, j) = a(k, j) - rk * cj
                a(i, j) = a(i, j) + ri * cj
                a(l, j) = a(l, j) - rl * cj
            END IF
            !..........FILL K NUCLIDE COLUMN.
            IF (k.le.isize) THEN
                IF (j.le.isize) a(j, k) = a(j, k) - rj * ck
                a(k, k) = a(k, k) + rk * ck
                a(i, k) = a(i, k) - ri * ck
                a(l, k) = a(l, k) + rl * ck
            END IF
            !..........FILL L NUCLIDE COLUMN.
            IF (j.le.isize) a(j, l) = a(j, l) - rj * cl
            IF (k.le.isize) a(k, l) = a(k, l) + rk * cl
            a(i, l) = a(i, l) - ri * cl
            a(l, l) = a(l, l) + rl * cl
        END IF !((ind.ne.0).and.(i.le.isize).and.(l.le.isize))
    END DO !n = 1,jsize

    !40--------PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION-----------

    bdln = 1.e-5 * (3. * hubcst)   !(10**(-5))*(Expansion rate).
    DO i = 1, isize
        i1 = isize1 - i            !Invert the rows.
        DO j = 1, isize
            j1 = isize1 - j          !Invert the columns.
            IF (dabs(a(j, i)).lt.bdln * y0(j1) / y0(i1)) THEN
                a(j, i) = 0.d0          !Set 0 if tiny.
            ELSE
                a(j, i) = a(j, i) * dt     !Bring dt over to other side.
            END IF
        END DO
        a(i, i) = 1.d0 + a(i, i)     !Add identity matrix to a-matrix.
        b(i1) = y0(i)             !Initial abundances.
    END DO

    !50--------SOLVE EQUATIONS TO GET DERIVATIVE------------------------------------

    !..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION.
    IF (loop.eq.1) THEN
        CALL eqslin(ip, ierror)
    ELSE
        CALL eqslin(0, ierror)
    END IF
    !..........OBTAIN DERIVATIVE.
    DO i = 1, isize
        yy(i) = yx(isize1 - i)     !Abundance at t+dt.
        dydt(i) = (yy(i) - y0(i)) / dt         !Take derivative.
    END DO

    !60--------POSSIBLE ERROR MESSAGES AND EXIT-------------------------------------

    IF (mbad.ne.0) THEN          !Problem in gaussian elimination.
        IF (mbad.eq.-1) WRITE (iw, 6000) ierror !Error message.
        IF (mbad.ge. 1) WRITE (iw, 6002) mbad   !Error message.
        6000   FORMAT (' ', '** y(', i2, ') fails to converge **')
        6002   FORMAT (' ', '** ', i2, ' th diagonal term equals zero **')
    END IF
    RETURN

    !----------REFERENCES-----------------------------------------------------------
    !     1) The coefficients are given in general as:
    !             ci = ri*(y(j)**rj)*(y(i)**(ri-1)*f(n)/
    !                  ((ri+rj)*fac(ri)*fac(rj))
    !             cj = rj*(y(i)**ri)*(y(j)**(rj-1)*f(n)/
    !                  ((ri+rj)*fac(ri)*fac(rj))
    !             ck = rk*(y(l)**rl)*(y(k)**(rk-1)*f(n)/
    !                  ((rk+rl)*fac(rk)*fac(rl))
    !             cl = rl*(y(k)**rk)*(y(l)**(rl-1)*f(n)/
    !                  ((rk+rl)*fac(rk)*fac(rl))
    !        in which fac(x) is the factorial of x.
    !     2) Form of reverse rate given in
    !        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
    !          tables 1B, 4B, 7B.
    !     3) Form of reverse rate given in
    !        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
    !          tables 2B, 3B, 5B, 6B, 8B, 9B, 10B.
    !     4) Form of reverse rate given in
    !        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
    !          table 11B.
    !     5) Form of reverse rate given in
    !        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
    !          tables 12B, 13B.

END


!========================IDENTIFICATION DIVISION================================

SUBROUTINE eqslin(icnvm, ierror)

    !----------LINKAGES.
    !     CALLED BY - [subroutine] sol
    !     CALLS     - none

    !----------REMARKS.
    !     Solves for new abundances using gaussian elimination
    !     with back substitution, no pivoting.

    !----------PARAMETERS.
    PARAMETER (nnuc = 26)          !Rank of matrix.
    PARAMETER (mord = 1)           !Higher order in correction.
    PARAMETER (eps = 2.e-4)        !Tolerance for convergence (.ge. 1.e-7).

    !----------COMMON AREAS.
    COMMON /compr/  cy, ct, t9i, t9f, ytmin, inc        !Computation parameters.
    COMMON /lncoef/ a, b, y                          !Lin eqn coefficients.
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags, counters.
    COMMON /runopt/ irun, isize, jsize               !Run options.


    !==========================DECLARATION DIVISION=================================

    !----------COMPUTATION PARAMETER.
    INTEGER inc                  !Accumulation increment.

    !----------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
    DOUBLE PRECISION a(nnuc, nnuc)!Coefficient array.
    REAL    b(nnuc)              !Right-hand vector w/o manipulation.
    REAL    y(nnuc)              !Solution vector.

    !----------COUNTERS AND FLAGS.
    INTEGER mbad                 !Indicates type of error.

    !----------RUN OPTION.
    INTEGER isize                !Number of nuclides in computation.

    !----------LOCAL MATRICES AND VECTORS.
    DOUBLE PRECISION a0(nnuc, nnuc)!Coefficient array w/o manipulation.
    DOUBLE PRECISION x(nnuc)     !Right-hand vector.

    !----------LOCAL COMPUTATION VARIABLES.
    DOUBLE PRECISION cx          !Scaling factor in triangularization.
    DOUBLE PRECISION sum         !Sum for backsubstitution.
    REAL   xdy                   !Relative error.

    !----------LOCAL COUNTERS.
    INTEGER nord                 !Order of correction.
    INTEGER icnvm                !Convergence monitor.
    INTEGER ierror               !ith nuclide fails to converge.


    !===========================PROCEDURE DIVISION==================================

    !10--------INITIALIZE VECTOR----------------------------------------------------

    !..........SET COUNTERS TO ZERO.
    nord = 0                     !No corrections yet.
    mbad = 0                     !No errors yet.
    !..........SET RIGHT-HAND AND SOLUTION VECTORS TO INITIAL VALUES.
    DO i = 1, isize
        x(i) = b(i)                !Right-hand vector.
        y(i) = 0.                  !Solution vector.
    END DO
    !..........SAVE MATRIX.
    IF (icnvm.eq.inc) THEN       !Monitor convergence.
        DO i = 1, isize
            DO j = 1, isize
                a0(j, i) = a(j, i)       !Initial value of coefficient array.
            END DO
        END DO
    END IF

    !20--------TRIANGULARIZE MATRIX AND SAVE OPERATOR-------------------------------

    !..........CHECK TO SEE THAT THERE ARE NO ZEROES AT PIVOT POINTS.
    DO i = 1, isize - 1
        IF (a(i, i).eq.0.d0) THEN   !Don't want to divide by zero.
            mbad = i                 !Position of zero coefficient.
            RETURN                   !Terminate matrix evaluation.
        END IF
        !..........TRIANGULARIZE MATRIX.
        DO j = i + 1, isize
            IF (a(j, i).ne.0.d0) THEN !Progress diagonally down the column.
                cx = a(j, i) / a(i, i)     !Scaling factor down the column.
                DO k = i + 1, isize       !Progress diagonally along row.
                    a(j, k) = a(j, k) - cx * a(i, k)  !Subtract scaled coeff along row.
                END DO
                a(j, i) = cx            !Scaled coefficient.
                !..........OPERATE ON RIGHT-HAND VECTOR.
                x(j) = x(j) - cx * x(i)  !Subtract off scaled coefficient.
            END IF
        END DO
    END DO

    !30--------DO BACK SUBSTITUTION-------------------------------------------------

    300  CONTINUE
    x(isize) = x(isize) / a(isize, isize)   !Solution for ultimate position.
    y(isize) = y(isize) + x(isize)
    DO i = isize - 1, 1, -1          !From i = penultimate to i = 1.
        sum = 0.d0
        DO j = i + 1, isize
            sum = sum + a(i, j) * x(j)  !Sum up all previous terms.
        END DO
        x(i) = (x(i) - sum) / a(i, i)
        y(i) = y(i) + x(i)         !Add difference to initial value.
    END DO

    !40--------TESTS AND EXITS------------------------------------------------------

    IF (icnvm.eq.inc) THEN
        DO i = 1, isize
            IF (y(i).ne.0.) THEN
                xdy = dabs(x(i) / y(i))  !Relative error.
                IF (xdy.gt.eps) THEN
                    IF (nord.lt.mord) THEN !Continue to higher orders.
                        nord = nord + 1
                        !..........FIND ERROR IN RIGHT-HAND VECTOR.
                        DO j = 1, isize
                            r = 0.d0         !Initialize r.
                            DO k = 1, isize
                                r = r + a0(j, k) * y(k) !Left side with approximate solution.
                            END DO
                            x(j) = b(j) - r  !Subtract difference from right side.
                        END DO
                        !..........OPERATE ON RIGHT-HAND VECTOR.
                        DO j = 1, isize - 1
                            DO k = j + 1, isize
                                x(k) = x(k) - a(k, j) * x(j)  !Subtract off scaled coefficient.
                            END DO
                        END DO
                        GO TO 300          !Go for another iteratiion.
                    ELSE
                        !..........NOT ENOUGH CONVERGENCE.
                        mbad = -1          !Signal error problem.
                        ierror = i         !ith nuclide for which x/y checked.
                        RETURN
                    END IF !(nord.lt.mord)
                END IF !(xdy.gt.eps)
            END IF !(y(i).ne.0)
        END DO !i = 1,isize
    END IF !(icnvm.eq.inc)
    RETURN                       !No more iterations & relative error small.

END