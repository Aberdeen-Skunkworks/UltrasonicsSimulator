#include "simulation.hpp"



Simulation::Simulation() {}



void Simulation::add_transducer(const Transducer t) {

    transducers.push_back(t);

}



void Simulation::add_transducers(const std::vector<Transducer> transducer_list) {

    for (Transducer t : transducer_list) {

	add_transducer(t);
	
    }

}



void Simulation::clear() {

    transducers.clear();
    
}



std::complex<double> Simulation::pressure_sum(const Eigen::Vector3d point, const double shift /* =0 */) const {
    
    std::complex<double> sum = 0+0j;

    for (const Transducer& t : transducers) {
	sum += t.pressure(point, shift);
    }

    return sum;

}



Eigen::Vector3cd Simulation::nablap_sum(const Eigen::Vector3d point, const double shift /* =0 */) const {
    
    Eigen::Vector3cd sum = {0+0j, 0+0j, 0+0j};

    for (const Transducer& t : transducers) {
	sum += t.nablap(point, shift);
    }

    return sum;

}



double Simulation::Gorkov_potential(const Particle& p) const {

    std::complex<double> pressure = pressure_sum(p.pos);
    Eigen::Vector3cd nablap = nablap_sum(p.pos);
    
    double omega = 2.0 * M_PI * frequency;
    Eigen::Vector3cd v = nablap / (1j * air_density * omega);
    Eigen::Vector3cd conjv = v.conjugate();
    std::complex<double> vdotconjv = v[0]*conjv[0] + v[1]*conjv[1] + v[2]*conjv[2];    

    double avg_psq = (pressure * std::conj(pressure)).real() / 2.0;
    double avg_vsq = vdotconjv.real() / 2.0;

    double k_0 = 1.0 / (air_density * std::pow(sound_speed_air, 2.0));
    double k_p = 1.0 / (p.density() * std::pow(sound_speed_particle, 2.0));
    double k_tilde = k_p / k_0;
    double f_1 = 1.0 - k_tilde;

    double rho_tilde = p.density() / air_density;
    double f_2 = 2.0 * (rho_tilde - 1.0) / (2.0 * rho_tilde + 1.0);

    // std::cout << p.volume() << ", " << p.mass << ", " << p.diameter << "\n";
    
    return p.volume() * (f_1 * 0.5 * k_0 * avg_psq - f_2 * 0.75 * air_density * avg_vsq) - p.mass * gravity * p.pos[2];
    
}



void Simulation::focus(const Eigen::Vector3d point) {

    for (Transducer& t : transducers) {
	
        Eigen::Vector3d d = point - t.pos;
	double dmag = d.norm();

	t.phi = (1 - std::fmod(dmag / t.wavelength, 1.0)) * 2 * M_PI;

	std::cout << "phi: " << t.phi << "\n";
	
    }

}



void Simulation::focus(const std::array<double, 3> array_point) {

    Eigen::Vector3d point({array_point[0], array_point[1], array_point[2]});
    
    for (Transducer& t : transducers) {
	
        Eigen::Vector3d d = point - t.pos;
	double dmag = d.norm();

	t.phi = (1 - std::fmod(dmag / t.wavelength, 1.0)) * 2 * M_PI;

	std::cout << "phi: " << t.phi << "\n";
	
    }

}



Field<double> Simulation::Gorkov_potential_field(const std::array<size_t, 3>& N, const std::array<double, 3>& L, const std::array<double, 3>& origin, const double particle_mass, const double particle_diameter) const {

    Field<double> gorkov(N, L, origin);

    for (size_t i = 0; i < N[0]; i++) {
	for (size_t j = 0; j < N[1]; j++) {
	    for (size_t k = 0; k < N[2]; k++) {
		double x = origin[0] + L[0] * (i+0.5) / N[0];
		double y = origin[1] + L[1] * (j+0.5) / N[1];
		double z = origin[2] + L[2] * (k+0.5) / N[2];

		const std::array<double, 3> particle_pos({x, y, z});
		Particle p(particle_pos, particle_mass, particle_diameter);

		gorkov({i, j, k}) = Gorkov_potential(p);
		// std::cout << x << ", " << y << ", " << z << " -> " << gorkov({i, j, k}) << "\n";
	    }
	}
    }

    return gorkov;

}



double Simulation::laplacian_sum(const std::vector<std::array<double, 3> >& points, const double width, const double particle_mass, const double particle_diameter, const bool x, const bool y, const bool z) const {

    double sum = 0;

    for (const std::array<double, 3>& p : points) {

	double x = p[0];
	double y = p[1];
	double z = p[2];
	
	const std::array<double, 3> p1_pos({x+width/2, y, z});
	Particle p1(p1_pos, particle_mass, particle_diameter);
	double f1 = Gorkov_potential(p1);
	
	const std::array<double, 3> p2_pos({x, y, z});
	Particle p2(p2_pos, particle_mass, particle_diameter);
	double f2 = Gorkov_potential(p2);
	
	const std::array<double, 3> p3_pos({x-width/2, y, z});
	Particle p3(p3_pos, particle_mass, particle_diameter);
	double f3 = Gorkov_potential(p3);

	sum += (f1 - 2*f2 + f3) / std::pow(width, 2);
	
    }
    
    return sum;

}



void Simulation::optimise_Gorkov_laplacian(const std::vector<std::array<double, 3> >& optimisation_points,
					   const double laplacian_width,
					   const double particle_mass,
					   const double particle_diameter) {
    
    nlopt::opt opt(nlopt::GN_ESCH, transducers.size());

    const std::vector<double> lb(transducers.size(), 0);
    const std::vector<double> ub(transducers.size(), 2*M_PI);

    LaplacianOptimisationData data(optimisation_points, laplacian_width, particle_mass, particle_diameter);
    
    Test tdata = {this, &data};

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    //opt.set_max_objective(optimise_laplacian_function, &data);
    opt.set_max_objective(owrapper, &tdata);

    opt.set_maxtime(1);

    opt.set_xtol_rel(1e-4);
    std::vector<double> x(transducers.size(), M_PI);
    double maxf;

    try{
    	nlopt::result result = opt.optimize(x, maxf);
    	// std::cout << "found minimum of " << std::setprecision(10) << minf << "\n";
	std::cout << "found max of " << maxf << "\n";
	// std::cout << "phases:" << "\n";
	// for (size_t i = 0; i < x.size(); i++) {
	//     std::cout << i << " - " << x[i] << "\n";
	// }
    }
    catch(std::exception &e) {
    	std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    // setting phi values to optimised values
    for (size_t i = 0; i < x.size(); i++) {
	transducers[i].phi = x[i];
    }
    
}



double Simulation::optimise_laplacian_function(const std::vector<double> &x, std::vector<double> &grad, void* opt_func_data) {
    LaplacianOptimisationData* d = reinterpret_cast<LaplacianOptimisationData*>(opt_func_data);
    std::vector<std::array<double, 3> > optimisation_points = d->optimisation_points;
    double laplacian_width = d->laplacian_width;
    double particle_mass = d->particle_mass;
    double particle_diameter = d->particle_diameter;

    for (size_t i = 0; i < x.size(); i++) {
	transducers[i].phi = x[i];
    }

    double lap = laplacian_sum(optimisation_points, laplacian_width, particle_mass, particle_diameter);
    
    if (optimisation_counter % 100 == 0) {
	std::cout << optimisation_counter << ", " << lap << "\n";
    }
    optimisation_counter++;
    
    return lap;
    
}



double owrapper(const std::vector<double>& x, std::vector<double>& grad, void* wrapper_func_data) {
    Test* t = reinterpret_cast<Test*>(wrapper_func_data);
    Simulation* sim = t->sim;
    return sim->optimise_laplacian_function(x, grad, t->data);
}
