#include "crpropa/Cosmology.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include <vector>
#include <cmath>
#include <stdexcept>

namespace crpropa {

/**
 @class Cosmology
 @brief Cosmology calculations
 */
struct Cosmology {
	double H0; // Hubble parameter at z=0
	double omegaM; // matter density parameter
    double omegaR0; // radiation density parameter
    double N0; // model parameter, default = 1.4
    double R0;

	static const int n;
	static const double zmin;
	static const double zmax;

	std::vector<double> Z;  // redshift
	std::vector<double> Dc; // comoving distance [m]
	std::vector<double> Dl; // luminosity distance [m]
	std::vector<double> Dt; // light travel distance [m]

	void update() {
		double dH = c_light; // Hubble distance

		std::vector<double> E(n);
		E[0] = 1;

		// Relation between comoving distance r and redshift z (cf. J.A. Peacock, Cosmological physics, p. 89 eq. 3.76)
		// dr = c / H(z) dz, integration using midpoint rule
		double dlz = log10(zmax) - log10(zmin); 
		for (int i = 1; i < n; i++) {
			Z[i] = zmin * pow(10, i * dlz / (n - 1)); // logarithmic even spacing
			double dz = (Z[i] - Z[i - 1]); // redshift step
			E[i] = sqrt((-2 * N0 * R0 / (3 * pow(3 - N0, 2) * omegaM)) *
			            ((N0 - 3) * omegaM * pow(1 + Z[i], 3.0 / N0) +
			             2 * (N0 - 2) * omegaR0 * pow(1 + Z[i], (N0 + 3) / N0)));
			Dc[i] = Dc[i - 1] + dH * dz * (1 / E[i] + 1 / E[i - 1]) / 2;
			Dl[i] = (1 + Z[i]) * Dc[i];
			Dt[i] = Dt[i - 1]
					+ dH * dz
							* (1 / ((1 + Z[i]) * E[i])
									+ 1 / ((1 + Z[i - 1]) * E[i - 1])) / 2;
		}
	}

	Cosmology() {
        // Cosmological parameters (K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014))
		H0 = 67.4 * 1000 * meter / second / Mpc; // default values
		omegaM = 0.315;
        omegaR0 = 5.373*1e-5;
        N0 = 1.4;
		R0 = -(3 * pow(3 - N0, 2) * pow(H0, 2) * omegaM / 
        ( 2 * N0 * ((N0 - 3) * omegaM + 2 * (N0 - 2) * omegaR0))) ;

		Z.resize(n);
		Dc.resize(n);
		Dl.resize(n);
		Dt.resize(n);

		Z[0] = 0;
		Dc[0] = 0;
		Dl[0] = 0;
		Dt[0] = 0;

		update();
	}

	void setParameters(double h, double oM, double oR0, double n) {
		H0 = h * 1e5 / Mpc;
		omegaM = oM;
        omegaR0 = oR0;
        N0 = n;
		update();
	}
};

const int Cosmology::n = 1000; 
const double Cosmology::zmin = 0.0001;
const double Cosmology::zmax = 100;

static Cosmology cosmology; // instance is created at runtime

void setCosmologyParameters(double h, double oM, double oR0, double n) {
	cosmology.setParameters(h, oM, oR0, n);
}

double hubbleRate(double z) {
    double omegaM = cosmology.omegaM;
    double omegaR0 = cosmology.omegaR0;
    double N0 = cosmology.N0;
    double R0 = cosmology.R0;

    return sqrt((-2 * N0 * R0 / (3 * pow(3 - N0, 2) * omegaM)) *
	            ((N0 - 3) * omegaM * pow(1 + z, 3.0 / N0) +
	             2 * (N0 - 2) * omegaR0 * pow(1 + z, (N0 + 3) / N0)));
}

double omegaR0() {
	return cosmology.omegaR0;
}

double omegaM() {
	return cosmology.omegaM;
}

double H0() {
	return cosmology.H0;
}

double R0() {
    return cosmology.R0;
}

double comovingDistance2Redshift(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dc.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dc, cosmology.Z);
}

double redshift2ComovingDistance(double z) {
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > cosmology.zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	return interpolate(z, cosmology.Z, cosmology.Dc);
}

double luminosityDistance2Redshift(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dl.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dl, cosmology.Z);
}

double redshift2LuminosityDistance(double z) {
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > cosmology.zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	return interpolate(z, cosmology.Z, cosmology.Dl);
}

double lightTravelDistance2Redshift(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dt.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dt, cosmology.Z);
}

double redshift2LightTravelDistance(double z) {
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > cosmology.zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	return interpolate(z, cosmology.Z, cosmology.Dt);
}

double comoving2LightTravelDistance(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dc.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dc, cosmology.Dt);
}

double lightTravel2ComovingDistance(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dt.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dt, cosmology.Dc);
}

} // namespace crpropa
