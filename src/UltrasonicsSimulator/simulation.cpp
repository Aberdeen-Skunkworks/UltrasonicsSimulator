#include "simulation.hpp"



Simulation::Simulation():
    gravity({0, -9.81, 0})
{}



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

    const std::complex<double> pressure = pressure_sum(p.pos);
    const Eigen::Vector3cd nablap = nablap_sum(p.pos);

    const double omega = 2.0 * M_PI * frequency;
    const Eigen::Vector3cd v = nablap / (1j * air_density * omega);
    const Eigen::Vector3cd conjv = v.conjugate();
    const std::complex<double> vdotconjv = v[0]*conjv[0] + v[1]*conjv[1] + v[2]*conjv[2];

    const double avg_psq = (pressure * std::conj(pressure)).real() / 2.0;
    const double avg_vsq = vdotconjv.real() / 2.0;

    const double k_0 = 1.0 / (air_density * std::pow(sound_speed_air, 2.0));
    const double k_p = 1.0 / (p.density() * std::pow(sound_speed_particle, 2.0));
    const double k_tilde = k_p / k_0;
    const double f_1 = 1.0 - k_tilde;

    const double rho_tilde = p.density() / air_density;
    const double f_2 = 2.0 * (rho_tilde - 1.0) / (2.0 * rho_tilde + 1.0);

    // std::cout << p.volume() << ", " << p.mass << ", " << p.diameter << "\n";

    return p.volume() * (f_1 * 0.5 * k_0 * avg_psq - f_2 * 0.75 * air_density * avg_vsq) - p.mass * gravity.dot(p.pos);

}



void Simulation::focus(const Eigen::Vector3d point) {

    for (Transducer& t : transducers) {

        Eigen::Vector3d d = point - t.pos;
	double dmag = d.norm();

	t.phi = (1 - std::fmod(dmag / t.wavelength, 1.0)) * 2 * M_PI;

    }

}



void Simulation::focus(const std::array<double, 3> array_point) {

    Eigen::Vector3d point({array_point[0], array_point[1], array_point[2]});

    for (Transducer& t : transducers) {

        Eigen::Vector3d d = point - t.pos;
	double dmag = d.norm();

	t.phi = (1 - std::fmod(dmag / t.wavelength, 1.0)) * 2 * M_PI;

	// std::cout << "phi: " << t.phi << "\n";

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



double Simulation::laplacian_sum(const std::vector<std::array<double, 3> >& points,
				 const double width,
				 const double particle_mass,
				 const double particle_diameter,
				 const bool dx,
				 const bool dy,
				 const bool dz) const {

    double sum = 0;

    std::array<bool, 3> d({dx, dy, dz});

    for (const std::array<double, 3>& p : points) {

	for (size_t n = 0; n < d.size(); n++) {

	    if (d[n]) {
		// std::cout << "getting laplacian in " << n << "dir" << "\n";

		std::array<double, 3> p_hi = p;
		p_hi[n] += width/2;
		std::array<double, 3> p_lo = p;
		p_lo[n] -= width/2;

		// double x = p[0];
		// double y = p[1];
		// double z = p[2];

		// const std::array<double, 3> p1_pos({x-width/2, y, z});
		const std::array<double, 3> p1_pos(p_lo);
		Particle p1(p1_pos, particle_mass, particle_diameter);
		double f1 = Gorkov_potential(p1);

		// const std::array<double, 3> p2_pos({x, y, z});
		const std::array<double, 3> p2_pos(p);
		Particle p2(p2_pos, particle_mass, particle_diameter);
		double f2 = Gorkov_potential(p2);

		// const std::array<double, 3> p3_pos({x+width/2, y, z});
		const std::array<double, 3> p3_pos(p_hi);
		Particle p3(p3_pos, particle_mass, particle_diameter);
		double f3 = Gorkov_potential(p3);

		// std::cout << "1 " << f1 << "\t2 " << f2 << "\t3 " << f3 << "\n";

		sum += (f1 - 2*f2 + f3) / std::pow(width, 2);
		// std::cout << n << "\n" << p_lo[0] << " " << p_lo[1] << " " << p_lo[2];
		// std::cout << "\n" << p[0] << " " << p[1] << " " << p[2];
		// std::cout << "\n" << p_hi[0] << " " << p_hi[1] << " " << p_hi[2] << "\n";
		// std::cout << "value: " << f1 << ", " << f2 << ", " << f3 << ", " << (f1 - 2*f2 + f3) / std::pow(width, 2) << "\n";

	    }
	}
    }

    // std::cout << "sum " << sum << "\n";

    return sum;

}



void Simulation::optimise_Gorkov_laplacian(const std::vector<std::array<double, 3> >& optimisation_points,
					   const double laplacian_width,
					   const double particle_mass,
					   const double particle_diameter,
					   const bool dx,
					   const bool dy,
					   const bool dz,
					   const double max_optimisation_time,
					   const size_t max_iterations,
					   const double xtol,
					   const std::string algo) {

    optimisation_counter = 0;

    if (algo == "LBFGS"){
	std::cout << "Using LBFGS algorithm" << "\n";

	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-6;
	param.max_iterations = max_iterations;

	LBFGSpp::LBFGSSolver<double> solver(param);
	LBFGS_laplacian_optimiser fun(this, optimisation_points, laplacian_width, particle_mass, particle_diameter, dx, dy, dz);

	// Initial guess
	Eigen::VectorXd x = Eigen::VectorXd::Zero(transducers.size());
	for (size_t i = 0; i < transducers.size(); i++) {
	    x[i] = transducers[i].phi;
	    //x[i] = 3.14;
	}
	// std::cout << "x init\n" << x << "\n";

	// x will be overwritten to be the best point found
	double fx;
	int niter = solver.minimize(fun, x, fx);

	std::cout << niter << " iterations" << std::endl;
	std::cout << "x = \n" << x.transpose() << std::endl;
	std::cout << "f(x) = " << -fx << std::endl;
    }

    else if (algo == "IPOPT"){
	
	// Create a new instance of your nlp (use a SmartPtr, not raw)
	Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Ipopt::LaplacianOptimiser(this, optimisation_points, laplacian_width, particle_mass, particle_diameter, dx, dy, dz);

	// Create a new instance of IpoptApplication (use a SmartPtr,
	//  not raw) We are using the factory, since this allows us to
	//  compile this example with an Ipopt Windows DLL
	Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

	// Change some options
	// Note: The following choices are only examples, they might not be
	//       suitable for your optimization problem.
	app->Options()->SetNumericValue("tol", xtol);
	app->Options()->SetIntegerValue("max_iter", max_iterations);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");

	// Intialize the IpoptApplication and process the options
	Ipopt::ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Ipopt::Solve_Succeeded) {
	    printf("\n\n*** Error during initialization!\n");
	    // return (void) status;
	}

	// Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

	if (status == Ipopt::Solve_Succeeded) {
	    printf("\n\n*** The problem solved!\n");
	}
	else {
	    printf("\n\n*** The problem FAILED!\n");
	}

	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically 
	// be deleted.

	// return (int) status;
	
    }

    else if (algo == "LN_BOBYQA" or algo == "GN_ESCH"){

	nlopt::algorithm nlopt_enum;
	if (algo == "LN_BOBYQA") {
	    std::cout << "Using LN_BOBYQA (Local Optimisation)" << "\n";
	    nlopt_enum = nlopt::LN_BOBYQA;
	}
	else {
	    std::cout << "Using GN_ESCH (Genetic Algorithm)" << "\n";
	    nlopt_enum = nlopt::GN_ESCH;
	}

	nlopt::opt opt(nlopt_enum, transducers.size());

	const std::vector<double> lb(transducers.size(), 0);
	const std::vector<double> ub(transducers.size(), 2*M_PI);

	LaplacianOptimisationData data(optimisation_points, laplacian_width, particle_mass, particle_diameter, dx, dy, dz);

	LaplacianSimCombo tdata = {this, &data};

	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.set_max_objective(optimise_laplacian_function_wrapper, &tdata);

	opt.set_maxtime(max_optimisation_time);

	opt.set_xtol_rel(xtol);

	std::vector<double> x;
	for (const Transducer& t : transducers) {
	    x.push_back(t.phi);
	}

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

    else {
	std::cout << "Error: Unrecognised algorithm: " << algo << "\n";
    }

}



double Simulation::optimise_laplacian_function(const std::vector<double> &x, std::vector<double> &grad, void* opt_func_data) {
    LaplacianOptimisationData* d = reinterpret_cast<LaplacianOptimisationData*>(opt_func_data);
    std::vector<std::array<double, 3> > optimisation_points = d->optimisation_points;
    double laplacian_width = d->laplacian_width;
    double particle_mass = d->particle_mass;
    double particle_diameter = d->particle_diameter;
    bool dx = d->dx;
    bool dy = d->dy;
    bool dz = d->dz;

    for (size_t i = 0; i < x.size(); i++) {
	transducers[i].phi = x[i];
    }

    double lap = laplacian_sum(optimisation_points, laplacian_width, particle_mass, particle_diameter, dx, dy, dz);

    if (optimisation_counter % 100 == 0) {
	std::cout << optimisation_counter << ", " << lap << "\n";
    }
    optimisation_counter++;

    return lap;

}



double optimise_laplacian_function_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* wrapper_func_data) {

    LaplacianSimCombo* c = reinterpret_cast<LaplacianSimCombo*>(wrapper_func_data);

    Simulation* sim = c->sim;

    return sim->optimise_laplacian_function(x, grad, c->data);

}



LBFGS_laplacian_optimiser::LBFGS_laplacian_optimiser(Simulation* sim, const std::vector<std::array<double, 3> > opt_points, const double width, const double mass, const double diameter, const bool dx, const bool dy, const bool dz):
    sim(sim),
    opt_points(opt_points),
    width(width),
    mass(mass),
    diameter(diameter),
    dx(dx),
    dy(dy),
    dz(dz)
{}



double LBFGS_laplacian_optimiser::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

    // std::cout << "x's\n";
    for (size_t i = 0; i < x.size(); i++) {
	// std::cout << x[i] << ", ";
	sim->transducers[i].phi = x[i];
    }
    // std::cout << "\n";

    // This optimiser can only minimise, so all of the laplacian sums
    // have a minus before them as we want to maximise this

    const double lap = - sim->laplacian_sum(opt_points, width, mass, diameter, dx, dy, dz);

    for (size_t i = 0; i < x.size(); i++) {

	const double incr = 2e-6;

	sim->transducers[i].phi += incr;
	const double lap_hi = - sim->laplacian_sum(opt_points, width, mass, diameter, dx, dy, dz);

	sim->transducers[i].phi -= incr;

	grad[i] = (lap_hi - lap) / incr;

    }

    if (sim->optimisation_counter % 100 == 0) {
	std::cout << sim->optimisation_counter << ", " << lap << "\n";
    }
    sim->optimisation_counter++;

    std::cout << "lap: " << -lap << "\n";

    return lap;

}



namespace Ipopt {

    bool LaplacianOptimiser::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
					  Index& nnz_h_lag, IndexStyleEnum& index_style)
    {
	// we have a number of variables equal to the number of
	// transducers in the simulation
	n = sim->transducers.size();

	// we have no constraints so we set m to 0, and jacobian and
	// hessian to have no non-zero entries
	m = 0;

	nnz_jac_g = 0;

	nnz_h_lag = 0;

	index_style = TNLP::C_STYLE;

	return true;

    }



    bool LaplacianOptimiser::get_bounds_info(Index n, Number* x_l, Number* x_u,
					     Index m, Number* g_l, Number* g_u)
    {
	// here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
	// If desired, we could assert to make sure they are what we think they are.
	assert(n == sim->transducers.size());
	assert(m == 0);

	// the variables have lower bounds of 1
	for (Index i = 0; i < sim->transducers.size(); i++) {
	    x_l[i] = 0.0;
	}

	// the variables have upper bounds of 5
	for (Index i = 0; i < sim->transducers.size(); i++) {
	    x_u[i] = 2*M_PI;
	}

	return true;
    }


    bool LaplacianOptimiser::get_starting_point(Index n, bool init_x, Number* x,
						bool init_z, Number* z_L, Number* z_U,
						Index m, bool init_lambda,
						Number* lambda)
    {
	// Here, we assume we only have starting values for x, if you code
	// your own NLP, you can provide starting values for the dual variables
	// if you wish to use a warmstart option
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	// initialize to the given starting point
	for (Index i = 0; i < sim->transducers.size(); i++) {
	    x[i] = fmod(sim->transducers[i].phi, 2*M_PI);
	}

	return true;
    }



    bool LaplacianOptimiser::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
    {
	assert(n == sim->transducers.size());

	// std::cout << "x's\n";
	for (size_t i = 0; i < sim->transducers.size(); i++) {
	    // std::cout << x[i] << ", ";
	    sim->transducers[i].phi = x[i];
	}
	// std::cout << "\n";

	obj_value = - sim->laplacian_sum(opt_points, width, mass, diameter, dx, dy, dz);

	return true;
    }



    bool LaplacianOptimiser::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
    {
	assert(n == sim->transducers.size());

	// std::cout << "x's\n";
	for (size_t i = 0; i < sim->transducers.size(); i++) {
	    // std::cout << x[i] << ", ";
	    sim->transducers[i].phi = x[i];
	}
	// std::cout << "\n";

	const double f = - sim->laplacian_sum(opt_points, width, mass, diameter, dx, dy, dz);

	for (size_t i = 0; i < sim->transducers.size(); i++) {

	    const double incr = 2e-6;

	    sim->transducers[i].phi += incr;
	    const double lap_hi = - sim->laplacian_sum(opt_points, width, mass, diameter, dx, dy, dz);

	    sim->transducers[i].phi -= incr;

	    grad_f[i] = (lap_hi - f) / incr;

	}

	return true;
    }



    /** Method to return the constraint residuals */
    bool LaplacianOptimiser::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
    {return true;}

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    bool LaplacianOptimiser::eval_jac_g(Index n, const Number* x, bool new_x,
					Index m, Index nele_jac, Index* iRow, Index *jCol,
					Number* values)
    {return true;}

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    bool LaplacianOptimiser::eval_h(Index n, const Number* x, bool new_x,
				    Number obj_factor, Index m, const Number* lambda,
				    bool new_lambda, Index nele_hess, Index* iRow,
				    Index* jCol, Number* values)
    {return true;}



    void LaplacianOptimiser::finalize_solution(SolverReturn status,
					       Index n, const Number* x, const Number* z_L,
					       const Number* z_U, Index m, const Number* g,
					       const Number* lambda, Number obj_value,
					       const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
    {
	// here is where we would store the solution to variables, or write to a file, etc
	// so we could use the solution.

	// std::cout << "x's\n";
	for (size_t i = 0; i < sim->transducers.size(); i++) {
	    // std::cout << x[i] << ", ";
	    sim->transducers[i].phi = x[i];
	}
	// std::cout << "\n";

	// For this example, we write the solution to the console
	printf("\n\nSolution of the primal variables, x\n");
	for (Index i=0; i<n; i++) {
	    printf("x[%d] = %e\n", i, x[i]);
	}

	printf("\n\nObjective value\n");
	printf("f(x*) = %e\n", obj_value);
    }

    

    LaplacianOptimiser::LaplacianOptimiser(Simulation* sim, const std::vector<std::array<double, 3> > opt_points, const double width, const double mass, const double diameter, const bool dx, const bool dy, const bool dz):
	sim(sim),
	opt_points(opt_points),
	width(width),
	mass(mass),
	diameter(diameter),
	dx(dx),
	dy(dy),
	dz(dz)
    {}


    
    LaplacianOptimiser::~LaplacianOptimiser() {}
    
}
