#include "ParticleSystem.hpp"



ParticleSystem::ParticleSystem() {}



void ParticleSystem::add_transducer(const Transducer t) {

    transducers.push_back(t);

}



void ParticleSystem::add_transducers(const std::vector<Transducer> transducer_list) {

    for (Transducer t : transducer_list) {

	add_transducer(t);
	
    }

}



void ParticleSystem::clear() {

    transducers.clear();
    
}



std::complex<double> ParticleSystem::pressure_sum(const Eigen::Vector3d point, const double shift /* =0 */) const {
    
    std::complex<double> sum = 0+0j;

    for (const Transducer& t : transducers) {
	sum += t.pressure(point, shift);
    }

    return sum;

}



Eigen::Vector3cd ParticleSystem::nablap_sum(const Eigen::Vector3d point, const double shift /* =0 */) const {
    
    Eigen::Vector3cd sum = {0+0j, 0+0j, 0+0j};

    for (const Transducer& t : transducers) {
	sum += t.nablap(point, shift);
    }

    return sum;

}



double ParticleSystem::Gorkov_potential(const Particle& p) const {

    std::complex<double> pressure = pressure_sum(p.pos);
    Eigen::Vector3cd nablap = nablap_sum(p.pos);
    
    double omega = 2.0 * M_PI * frequency;
    Eigen::Vector3cd v = nablap / (1j * air_density * omega);
    Eigen::Vector3cd conjv = v.conjugate();
    std::complex<double> vdotconjv = v[0]*conjv[0] + v[1]*conjv[1] + v[2]*conjv[2];    

    double avg_psq = (pressure * std::conj(pressure)).real() / 2.0;
    double avg_vsq = vdotconjv.real() / 2.0;
    
    std::complex<double> test = v[0]*conjv[0] + v[1]*conjv[1] + v[2]*conjv[2];

    double k_0 = 1.0 / (air_density * std::pow(sound_speed_air, 2.0));
    double k_p = 1.0 / (p.density() * std::pow(sound_speed_particle, 2.0));
    double k_tilde = k_p / k_0;
    double f_1 = 1.0 - k_tilde;

    double rho_tilde = p.density() / air_density;
    double f_2 = 2.0 * (rho_tilde - 1.0) / (2.0 * rho_tilde + 1.0);

    //std::cout << p.volume() * (f_1 * 0.5 * k_0 * avg_psq - f_2 * 0.75 * air_density * avg_vsq) << std::endl;

    // std::cout << "pvolume " << p.volume() << "\nf_1 " << f_1 << "\nk_0 " << k_0 << "\navg_psq " << avg_psq << "\nf_2 " << f_2 << "\nair_density " << air_density << "\navg_vsq " << avg_vsq << "\npmass " << p.mass << "\ng " << gravity << "\npos2 " << p.pos[2] << std::endl;
    
    return p.volume() * (f_1 * 0.5 * k_0 * avg_psq - f_2 * 0.75 * air_density * avg_vsq) - p.mass * gravity * p.pos[2];
    
}



void ParticleSystem::focus(const Eigen::Vector3d point) {

    for (Transducer& t : transducers) {
	
        Eigen::Vector3d d = point - t.pos;
	double dmag = d.norm();

	t.phi = (1 - std::fmod(dmag / t.wavelength, 1.0)) * 2 * M_PI;

	std::cout << "phi: " << t.phi << "\n";
	
    }

}
