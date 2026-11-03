*       Chain-inert pairs: a bound KS pair with a < CHINRT*RMIN that is also
*       below 0.01 of the chain scale enters a chain as a dormant c.m. body
*       instead of being resolved by the chain integrator (the standard
*       inert-binary path, no longer restricted to KZ(26) < 2 and KZ(27) > 0).
*       The chain is terminated and the pair renewed when a member comes
*       within CHINRR semi-major axes of it (the frozen pair is no longer
*       adiabatically isolated) or when its perturbation exceeds CHINGM.
*       The registry keeps the component masses and the semi-major axis of
*       the dormant pairs by the name of their c.m. body (at most two).
      PARAMETER  (CHINRT = 1.0D-03, CHINRR = 30.0D0, CHINGM = 1.0D-02)
      COMMON/CHDORM/  BDORM(2,2),ADORM(2),NMDORM(2),NDORM
