
//==========================================================================
//=== Physical functions.
//==========================================================================

/** Calculates the blackbody radiation (the Planck function).
    Each row of the returned matrix corresponds to a frequency, while each
    column corresponds to a temperature.
    A vector version, taking a single temperature, exists also.

    @param B       output: the blackbody radiation
    @param f       a frequency grid
    @param t       a temperature profile

    @author Patrick Eriksson 08.04.2000 */
void planck (
              MATRIX&     B, 
        const VECTOR&     f,
        const VECTOR&     t );
void planck (
              VECTOR&     B, 
        const VECTOR&     f,
        const Numeric&    t );



//==========================================================================
//=== Tangent altitudes.
//==========================================================================

/** Calculates the geometrical tangent altitude.
    That is, refraction is neglected.

    @return        the tangent altitude
    @param view    the angle between zenith and the LOS
    @param z_plat  the platform altitude

    @author Patrick Eriksson 08.04.2000 */
Numeric ztan_geom(
        const Numeric&     view,
        const Numeric&     z_plat );



//==========================================================================
//=== Conversion and interpolation of pressure and altitude grids.
//==========================================================================

/** Converts an altitude vector to pressures.
    The log. of the pressures are interpolated linearly, 
    i.e. (Matlab notation):

      p = exp(interp1(z0,log(p0),z,'linear'))

    @param p       output: the pressures at z
    @param z0      original altitude grid
    @param p0      original pressure grid
    @param z       new altitude grid

    @author Patrick Eriksson 10.04.2000 */
void z2p(
              VECTOR&     p,
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z );

/** Interpolates a vertical profile at a new set of pressures.
    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      p = interp1(log(p0),x,log(p)'linear')

    @param x       output: the interpolated values at p
    @param p0      original pressure grid
    @param x0      the profile to be interpolated
    @param p       new pressure grid

    @author Patrick Eriksson 12.04.2000 */
void interpp(
              VECTOR&     x,
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p );

