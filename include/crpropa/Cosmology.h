#ifndef CRPROPA_COSMOLOGY_H
#define CRPROPA_COSMOLOGY_H

namespace crpropa {
/**
 * \addtogroup PhysicsDefinitions
 * @{
 */

/**
 @file
 @brief Cosmology functions
 */

/**
 Set the cosmological parameters for a flat universe. To ensure flatness omegaL
 is set to 1 - omegaMatter (f(r) does not have the omegaL parameter).
 @param hubbleParameter	dimensionless Hubble parameter, default = 0.673
 @param omegaMatter		matter parameter, default = 0.315
 @param omegaRadiation			radiation parameter, default = 5.373*1e-5
 @param n			model parameter, default = 1.4
 */
void setCosmologyParameters(double hubbleParameter, double omegaMatter, double omegaRadiation, double n);

/**
 Hubble rate at given redshift
 H(z) = sqrt((-2 * n * R0 / (3 * pow(3 - n, 2) * omegaM)) * ( (n - 3) * omegaM *
 pow(1 + z, 3.0 / n) + 2 * (n - 2) * omegaR * pow(1 + z, (n + 3) / n ) ))
 */
double hubbleRate(double redshift = 0);

// Returns the radiation density parameter
double omegaR0();

// Returns the matter density parameter
double omegaM();

// Returns the hubble parameter
double H0();

double R0();

/**
 Redshift of a comoving object at a given comoving distance to an observer at z = 0.
 d_comoving(z) = c/H0 * int_0^z dz' / E(z')
 */
double comovingDistance2Redshift(double distance);

/**
 Comoving distance between an observer at z = 0 and a comoving object at z.
 d_comoving(z) = c/H0 * int_0^z dz' / E(z')
 */
double redshift2ComovingDistance(double redshift);

/**
 Redshift of a comoving object at a given luminosity distance to an observer at z = 0.
 d_luminosity(z) = (1 + z) * d_comoving(z)
 */
double luminosityDistance2Redshift(double distance);

/**
 Luminosity distance between an observer at z = 0 and a comoving object at z.
 d_luminosity(z) = (1 + z) * d_comoving(z)
 */
double redshift2LuminosityDistance(double redshift);

/**
 Redshift of a comoving object at a given light travel distance to an observer at z = 0.
 d_lighttravel(z) = c/H0 * int_0^z dz' / ((1 + z')  *  E(z'))
 */
double lightTravelDistance2Redshift(double distance);

/**
 Light travel distance between an observer at z = 0 and a comoving object at z.
 d_lighttravel(z) = c/H0 * int_0^z dz' / ((1 + z')  *  E(z'))
 */
double redshift2LightTravelDistance(double redshift);

// Conversion from comoving distance to light travel distance.
double comoving2LightTravelDistance(double distance);

// Conversion from light travel distance to comoving distance.
double lightTravel2ComovingDistance(double distance);

/** @}*/
} // namespace crpropa

#endif // CRPROPA_COSMOLOGY_H
