#include "transducer.hpp"



Transducer::Transducer(const Eigen::Vector3d pos, const Eigen::Vector3d dir, const double phi /* =0 */):
    pos(pos), director(dir / dir.norm()), phi(phi), k(2 * M_PI / wavelength)
{
    
    double l = director.norm();
    if (l == 0) {
	std::cout << "Cannot use a zero director" << std::endl;
	throw "Cannot use a zero director";
    }

}



Transducer::Transducer(const std::array<double, 3> pos, const std::array<double, 3> dir, const double phi /* =0 */):
    Transducer(Eigen::Vector3d({pos[0], pos[1], pos[2]}), Eigen::Vector3d({dir[0], dir[1], dir[2]}), phi)
{}



double Transducer::PktoPkA(const double d) const {
    return d * 18;
}



std::complex<double> Transducer::pressure(const Eigen::Vector3d point, const double shift /* =0 */) const {

    // checks
    const Eigen::Vector3d rs = point - pos;
    const double mag_rs = rs.norm();    
    if (mag_rs < 0.001) {
	return 0+0j;
    }    
    const double theta = std::acos(rs.dot(director) / mag_rs);
    if (theta > M_PI/2) {
	return 0+0j;
    }
    
    const std::complex<double> exponential = std::exp((std::complex<double>)(1j * (phi + k * mag_rs + shift))) / mag_rs;

    const double x = k * piston_radius * std::sin(theta);

    double frac;
    if (x < 1e-6) {
    	frac = PktoPkA(p0);
    }
    else {
    	frac = PktoPkA(p0 * std::sin(x) / x);
    }

    return exponential * frac;

}



Eigen::Vector3cd Transducer::nablap(const Eigen::Vector3d point, const double shift /* =0 */) const {

    // checks
    const Eigen::Vector3d rs = point - pos;
    const double mag_rs = rs.norm();
    if (mag_rs < 0.001) {
	return {0j, 0j, 0j};
    }

    const double theta = std::acos(rs.dot(director) / mag_rs);    
    if (theta > M_PI/2) {
	return {0j, 0j, 0j};
    }

    const Eigen::Vector3d rs_hat = rs / mag_rs;

    const std::complex<double> exponent_term = std::exp((std::complex<double>)1j * (phi + k * mag_rs + shift)) / mag_rs;

    const Eigen::Vector3cd d_exponential_term_dr = rs_hat * (std::complex<double>)1j * k * std::exp((std::complex<double>)1j * (phi + k * mag_rs + shift)) / mag_rs - std::exp((std::complex<double>)1j * (phi + k * mag_rs + shift)) * rs_hat / std::pow(mag_rs, 2);
    if (theta < 1e-4) {
    	const double fraction_term = PktoPkA(p0);
    	return fraction_term * d_exponential_term_dr;
    }
    else {
    	const double fraction_term = PktoPkA(p0 * std::sin(k * piston_radius * std::sin(theta)) / (k * piston_radius * std::sin(theta)));
    	const double numerator = PktoPkA(p0 * std::sin(k * piston_radius * std::sin(theta)));
    	const double denominator = k * piston_radius * std::sin(theta);
        const Eigen::Vector3cd d_denominator_dr = k * piston_radius * (-rs_hat.dot(director)) / std::sqrt(1 - rs_hat.dot(director)) * (director - rs_hat * director.dot(rs_hat)) / mag_rs;
    	const Eigen::Vector3cd d_numerator_dr = PktoPkA(p0) * std::cos(k * piston_radius * std::sin(theta)) * d_denominator_dr;
    	const Eigen::Vector3cd d_fraction_dr = (d_numerator_dr * denominator - numerator * d_denominator_dr) / std::pow(denominator, 2);
    	return d_fraction_dr * exponent_term + fraction_term * d_exponential_term_dr;
    }
}

    
// auto Transducer::pressure_equation() const {

//     using namespace sym;
    
//     Var<vidx<'x'> > x;
//     Var<vidx<'y'> > y;
//     Var<vidx<'z'> > z;
//     Var<vidx<'i'> > i;

//     auto rsx = x - pos(0);
//     auto rsy = y - pos(1);
//     auto rsz = z - pos(2);
    
//     auto dr = pow(rsx * rsx + rsy * rsy + rsz * rsz, C<1>()/C<2>());

//     auto sin_theta = pow(C<1>() - pow((rsx*director[0]+rsy*director[1]+rsz*director[2]) / dr, C<2>()), C<1>()/C<2>());
    
//     auto h = k * piston_radius * sin_theta;
//     auto Df = sin(h) / h;

//     auto P = PktoPkA(p0) * Df / dr * exp(i * (phi + k * dr));

//     return P;
    
// }


    
// std::complex<double> Transducer::stator_pressure(const Eigen::Vector3d point, const double shift /* =0 */) const {

//     // checks
//     const Eigen::Vector3d rs = point - pos;
//     const double mag_rs = rs.norm();    
//     if (mag_rs < 0.001) {
// 	return 0+0j;
//     }    
//     const double theta = std::acos(rs.dot(director) / mag_rs);
//     if (theta > M_PI/2) {
// 	return 0+0j;
//     }
    
//     using namespace sym;
//     auto P = simplify(pressure_equation());

//     Var<vidx<'x'> > x;
//     Var<vidx<'y'> > y;
//     Var<vidx<'z'> > z;
//     Var<vidx<'i'> > i;

//     std::complex<double> imag_unit = 1j;
    
//     return sub(sub(sub(sub(P, x=point[0]), y=point[1]), z=point[2]), i=imag_unit);

// }



// Eigen::Vector3cd Transducer::stator_nablap(const Eigen::Vector3d point, const double shift /* =0 */) const {

//     // checks
//     const Eigen::Vector3d rs = point - pos;
//     const double mag_rs = rs.norm();
//     if (mag_rs < 0.001) {
// 	return {0j, 0j, 0j};
//     }

//     const double theta = std::acos(rs.dot(director) / mag_rs);    
//     if (theta > M_PI/2) {
// 	return {0j, 0j, 0j};
//     }

//     using namespace sym;
//     auto P = pressure_equation();

//     Var<vidx<'x'> > x;
//     Var<vidx<'y'> > y;
//     Var<vidx<'z'> > z;
//     Var<vidx<'i'> > i;

//     auto Px = derivative(P, x);
//     auto Py = derivative(P, y);
//     auto Pz = derivative(P, z);
    
//     std::complex<double> imag_unit = 1j;
    
//     return {sub(sub(sub(sub(Px, x=point[0]), y=point[1]), z=point[2]), i=imag_unit),
// 	    sub(sub(sub(sub(Py, x=point[0]), y=point[1]), z=point[2]), i=imag_unit),
// 	    sub(sub(sub(sub(Pz, x=point[0]), y=point[1]), z=point[2]), i=imag_unit)};
    	    
// }
