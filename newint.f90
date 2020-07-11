!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     Changes (to run inder unix f77): 
!     -------------------------------
!     Programme renamed from nucint.for to newint.f
!     COMMON /therm/  -> COMMON /thermcb/
!     COMMON /bessel/ -> COMMON /besselcb/
!     COMMON /check/  -> COMMON /checkcb/
!     COMMON /time/   -> COMMON /ttime/
!
!     Updated Correction to helium-4 abundance (using fitted rates and  
!     *smallest* step-size) - Sarkar, Rep. Prog Phys. 59, 1493 (1996):
!
!     Y_p = Y_p - 0.0025 -> Y_p = Y_p - 0.0003
!
!
!     output nucint.dat  -> newint.dat
!
!========================IDENTIFICATION DIVISION================================

SUBROUTINE check

    !----------REMARKS.
    !     This is an interface subroutine,
    !     a flexible module which allows user to manipulate physical quantities
    !     of interest at certain key points during the computer run.
    !     Included within this subroutine is a roster of all global variables
    !     and their respective COMMON areas.
    !     Current instructions accumulate abundances of deuterium, helium-3,
    !     helium-4, and lithium-7 for eventual plotting, taking into account
    !     the contribution of beryllium-7 to lithium-7 and tritium to helium-3.

    !----------PARAMETERS.
    PARAMETER (nrec = 88)          !Number of nuclear reactions.
    PARAMETER (nvar = 29)          !Number of variables to be evolved.
    PARAMETER (nnuc = 26)          !Number of nuclides in calculation.
    PARAMETER (itmax = 40)         !Maximum # of lines to be printed.

    !----------COMMON AREAS.
    COMMON /recpr0/ reacpr                         !Reaction parameter values.
    COMMON /recpr/  iform, ii, jj, kk, ll, rev, q9       !Reaction parameter names.
    COMMON /rates/  f, r                            !Reaction rates.
    COMMON /evolp1/ t9, hv, phie, y                   !Evolution parameters.
    COMMON /evolp2/ dt9, dhv, dphie, dydt             !Evolution parameters.
    COMMON /evolp3/ t90, hv0, phie0, y0               !Evolution parameters.
    COMMON /compr0/ cy0, ct0, t9i0, t9f0, ytmin0, inc0  !Default comp parameters.
    COMMON /compr/  cy, ct, t9i, t9f, ytmin, inc        !Computation parameters.
    COMMON /modpr0/ c0, cosmo0, xi0                  !Default model parameters.
    COMMON /modpr/  g, tau, xnu, c, cosmo, xi           !Model parameters.
    COMMON /varpr0/ dt0, eta0                       !Default variationl params.
    COMMON /varpr/  dt1, eta1                       !Variational parameters.
    COMMON /ttime/   t, dt, dlt9dt                   !Time variables.
    COMMON /thermcb/  thm, hubcst                   !Dynamic variables.
    COMMON /endens/ rhone0, rhob0, rhob, rnb          !Energy densities.
    COMMON /lncoef/ a, b, yx                         !Linear eqn coefficients.
    COMMON /nucdat/ am, zm, dm                       !Nuclide data.
    COMMON /besselcb/ bl1, bl2, bl3, bl4, bl5, !Eval function bl(z).&
    bm1,bm2, bm3, bm4, bm5, !Eval function bm(z).&
    bn1,bn2, bn3, bn4, bn5            !Eval function bn(z).
    COMMON /kays/   bk0, bk1, bk2, bk3, bk4            !Coefficients K.
    COMMON /flags/  ltime, is, ip, it, mbad            !Flags,counters.
    COMMON /checkcb/  itime                        !Computation location.
    COMMON /outdat/ xout, thmout, t9out, tout, dtout, !Output data.&
    etaout,hubout
    COMMON /nupar/  t9mev, tnmev, tnu, cnorm, nu, rhonu !Neutrino parameters.
    COMMON /runopt/ irun, isize, jsize               !Run options.
    COMMON /outopt/ nout, outfile                   !Output option.


    !==========================DECLARATION DIVISION=================================

    !----------REACTION PARAMETER VALUES.
    REAL    reacpr(nrec, 8)       !Reaction parameters.

    !----------REACTION PARAMETER NAMES.
    INTEGER iform(nrec)          !Reaction type code (1-11).
    INTEGER ii(nrec)             !Incoming nuclide type (1-26).
    INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
    INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
    INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
    REAL    rev(nrec)            !Reverse reaction coefficient.
    REAL    q9(nrec)             !Energy released in reaction (in 10**9 K).

    !----------REACTION RATES.
    REAL    f(nrec)              !Forward reaction rate coefficients.
    REAL    r(nrec)              !Reverse reaction rate coefficients.

    !----------EVOLUTION PARAMETERS.
    REAL    t9                   !Temperature of photons (units of 10**9 K).
    REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
    REAL    phie                 !Chemical potential of electron.
    REAL    y(nnuc)              !Relative number abundances.

    !----------EVOLUTION PARAMETERS (DERIVATIVES).
    REAL    dt9                  !Change in temperature.
    REAL    dhv                  !Change in hv.
    REAL    dphie                !Change in chemical potential.
    REAL    dydt(nnuc)           !Change in relative number abundances.

    !----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
    REAL    y0(nnuc)             !Rel # abundances at end of 1st R-K loop.

    !----------DEFAULT COMPUTATION PARAMETERS.
    REAL    cy0                  !Default cy.
    REAL    ct0                  !Default ct.
    REAL    t9i0                 !Default t9i.
    REAL    t9f0                 !Default t9f.
    REAL    ytmin0               !Default ytmin.
    INTEGER inc0                 !Default accumulation increment.

    !----------COMPUTATION PARAMETERS.
    REAL    cy                   !Time step limiting constant on abundances.
    REAL    ct                   !Time step limiting constant on temperature.
    REAL    t9i                  !Initial temperature (in 10**9 K).
    REAL    t9f                  !Final temperature (in 10**9 K).
    REAL    ytmin                !Smallest abundances allowed.
    INTEGER inc                  !Accumulation increment.

    !----------DEFAULT MODEL PARAMETERS.
    REAL    c0(3)                !Default c.
    REAL    cosmo0               !Default cosmological constant.
    REAL    xi0(3)               !Default neutrino degeneracy parameters.

    !----------EARLY UNIVERSE MODEL PARAMETERS.
    REAL    g                    !Gravitational constant.
    REAL    tau                  !Neutron lifetime (sec).
    REAL    xnu                  !Number of neutrino species.
    REAL    c(3)                 !c(1) is variation of gravitational constant.&
    !c(2) is neutron half-life (min).&
    !c(3) is number of neutrino species.
    REAL    cosmo                !Cosmological constant.
    REAL    xi(3)                !Neutrino degeneracy parameters.&
    !xi(1) is e neutrino degeneracy parameter.&
    !xi(2) is m neutrino degeneracy parameter.&
    !xi(3) is t neutrino degeneracy parameter.

    !----------DEFAULT VARIATIONAL PARAMETERS.
    REAL    dt0                  !Default initial time step.
    REAL    eta0                 !Default baryon-to-photon ratio.

    !----------VARIATIONAL PARAMETERS.
    REAL    dt1                  !Initial time step.
    REAL    eta1                 !Baryon-to-photon ratio.

    !----------TIME VARIABLES.
    REAL    t                    !Time.
    REAL    dt                   !Time step.
    REAL    dlt9dt               !(1/t9)*d(t9)/d(t).

    !----------DYNAMIC VARIABLES.
    REAL    thm(14)              !Thermodynamic variables (energy densities).
    REAL    hubcst               !Expansion rate of the universe.

    !----------ENERGY DENSITIES.
    REAL    rhone0               !Initial electron neutrino energy density.
    REAL    rhob0                !Initial baryon energy density.
    REAL    rhob                 !Baryon energy density.
    REAL    rnb                  !Baryon energy density (ratio to init value).

    !----------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
    DOUBLE PRECISION a(nnuc, nnuc)!Relates y(t+dt) to y(t).
    REAL    b(nnuc)              !Contains y0 in inverse order.
    REAL    yx(nnuc)             !yy in reverse order.

    !----------NUCLIDE DATA.
    REAL    am(nnuc)             !Atomic number of nuclide.
    REAL    zm(nnuc)             !Charge of nuclide.
    REAL    dm(nnuc)             !Mass excess of nuclide.

    !----------EVALUATION OF FUNCTIONS bl,bm,bn.
    REAL    bl1, bl2, bl3, bl4, bl5  !Evaluation of function bl(z).
    REAL    bm1, bm2, bm3, bm4, bm5  !Evaluation of function bm(z).
    REAL    bn1, bn2, bn3, bn4, bn5  !Evaluation of function bn(z).

    !----------EVALUATION OF MODIFIED BESSEL FUNCTIONS.
    REAL    bk0, bk1, bk2, bk3, bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).

    !----------FLAGS AND COUNTERS.
    INTEGER ltime                !Indicates if output buffer printed.
    INTEGER is                   !# total iterations for particular model.
    INTEGER ip                   !# iterations after outputing a line.
    INTEGER it                   !# times accumulated in output buffer.
    INTEGER mbad                 !Indicates if gaussian elimination failed.

    !----------COMPUTATION LOCATION.
    INTEGER itime                !Time check.

    !----------OUTPUT ARRAYS.
    REAL    xout(itmax, nnuc)     !Nuclide mass fractions.
    REAL    thmout(itmax, 6)      !Thermodynamic variables.
    REAL    t9out(itmax)         !Temperature (in units of 10**9 K).
    REAL    tout(itmax)          !Time.
    REAL    dtout(itmax)         !Time step.
    REAL    etaout(itmax)        !Baryon to photon ratio.
    REAL    hubout(itmax)        !Expansion rate.

    !----------NEUTRINO PARAMETERS.
    REAL    t9mev                !Temperature (in units of MeV).
    REAL    tnmev                !Neutrino temperature (in units of MeV).
    REAL    tnu                  !Neutrino temperature.
    REAL    cnorm                !Normalizing constant.
    REAL    rhonu                !Neutrino energy density.
    INTEGER nu                   !Type of neutrino.

    !----------RUN OPTION.
    INTEGER irun                 !Run network size.
    INTEGER isize                !Number of nuclides in computation.
    INTEGER jsize                !Number of reactions in computation.

    !----------OUTPUT FILE STATUS.
    INTEGER nout                 !Number of output requests.
    LOGICAL outfile              !Indicates if output file used.


    !===========================PROCEDURE DIVISION==================================

    !10--------OPEN FILE------------------------------------------------------------

    IF (itime.eq.1) THEN         !Beginning of program.
        OPEN (unit = 3, file = 'newint.dat', status = 'new')
    END IF

    !20--------WRITE INTO FILE------------------------------------------------------

    IF (itime.eq.8) THEN         !Right after a run.
        xout(it, 8) = xout(it, 8) + xout(it, 9)  !Add beryllium to lithium.
        xout(it, 5) = xout(it, 5) + xout(it, 4)  !Add tritium to helium-3.
        xout(it, 6) = xout(it, 6) - 0.0003
        !                    ! my correction for fitted rates+coarse steps
        WRITE (3, 200) c(3), c(2), etaout(it), xout(it, 3), &
                xout(it, 5), xout(it, 6), xout(it, 8)&
        !Output N_nu, tau_n, eta, H2, He3, He4, and Li7.
        200    FORMAT (7(e13.5, ' '))
    END IF

    !30--------CLOSE FILE-----------------------------------------------------------

    IF (itime.eq.10) THEN        !End of program.
        CLOSE (unit = 3)
    END IF
    RETURN

    !----------REFERENCES-----------------------------------------------------------
    !     1) D.A. Dicus, E.W. Kolb, A.M. Gleeson, E.C.G. Sudarshan, V.L. Teplitz,
    !        M.S. Turner, Phys. Rev. D26 (1982) 2694. (Rad corr, Coulomb Corr)
    !     2) D. Seckel, Bartol preprint BA-93-16; G. Guyk and M.S. Turner,
    !         FERMILAB preprint FERMILAB-Pub-93/181-A. (Nucleon mass)
    !     4) S. Dodelson and M.S. Turner, Phys. Rev. D46 (1992) 3372; B. Fields,
    !         S. Dodelson and M.S. Turner, Phys. Rev. D47 (1993) 4309. (Nu heating)
    !     5) S. Sarkar, Rep. Prog Phys. 59 (1996) 1493 (review)

END