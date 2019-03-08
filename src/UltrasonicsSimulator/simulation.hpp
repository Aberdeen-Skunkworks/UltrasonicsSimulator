#pragma once
#include <vector>
#include <complex>

#include "nlopt.hpp"
#include "LBFGS.h"

#include "transducer.hpp"
#include "particle.hpp"
#include "field.hpp"
#include "optimisation_data.hpp"

#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"


class Simulation {
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    const Eigen::Vector3d gravity;

    
    const double frequency = 40000;

    
    const double sound_speed_air = 346; // m/s

    
    const double sound_speed_particle = 2600; // m/s

    
    const double air_density = 1.2; // kg/m^3


    size_t optimisation_counter = 0;

    
    Simulation();

    
    std::vector<Transducer> transducers;    

    
    void add_transducer(const Transducer);

    
    void add_transducers(const std::vector<Transducer>);

    
    void clear();

    
    std::complex<double> pressure_sum(const Eigen::Vector3d, const double shift=0) const;

    
    Eigen::Vector3cd nablap_sum(const Eigen::Vector3d, const double shift=0) const;
    

    double Gorkov_potential(const Particle&) const;

    
    void focus(const Eigen::Vector3d);

    
    void focus(const std::array<double, 3>);

    
    double laplacian_sum(const std::vector<std::array<double, 3> >&, const double, const double, const double, const bool x, const bool y, const bool z) const;

    
    Field<double> Gorkov_potential_field(const std::array<size_t, 3>&, const std::array<double, 3>&, const std::array<double, 3>&, const double, const double) const;
    
    
    void optimise_Gorkov_laplacian(const std::vector<std::array<double, 3> >&, const double, const double, const double, const bool, const bool, const bool, const double, const size_t, const double, const std::string);

    
    double optimise_laplacian_function(const std::vector<double>& , std::vector<double>&, void*);
    
};



double optimise_laplacian_function_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* wrapper_func_data);



class LBFGS_laplacian_optimiser {

    Simulation* sim;
    const std::vector<std::array<double, 3> > opt_points;
    const double width;
    const double mass;
    const double diameter;
    const bool dx;
    const bool dy;
    const bool dz;
    
public:

    LBFGS_laplacian_optimiser(Simulation*, const std::vector<std::array<double, 3> >, const double, const double, const double, const bool, const bool, const bool);
    
    double operator()(const Eigen::VectorXd&, Eigen::VectorXd&);
    
};



namespace Ipopt {

    class LaplacianOptimiser : public TNLP {

	Simulation* sim;
	const std::vector<std::array<double, 3> > opt_points;
	const double width;
	const double mass;
	const double diameter;
	const bool dx;
	const bool dy;
	const bool dz;

    public:
	/** default constructor */
	LaplacianOptimiser(Simulation*, const std::vector<std::array<double, 3> >, const double,
			   const double, const double, const bool, const bool, const bool);

	/** default destructor */
	virtual ~LaplacianOptimiser();

	/**@name Overloaded from TNLP */
	//@{
	/** Method to return some info about the nlp */
	virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				  Index& nnz_h_lag, IndexStyleEnum& index_style);

	/** Method to return the bounds for my problem */
	virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
				     Index m, Number* g_l, Number* g_u);

	/** Method to return the starting point for the algorithm */
	virtual bool get_starting_point(Index n, bool init_x, Number* x,
					bool init_z, Number* z_L, Number* z_U,
					Index m, bool init_lambda,
					Number* lambda);

	/** Method to return the objective value */
	virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	/** Method to return the gradient of the objective */
	virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	/** Method to return the constraint residuals */
	virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	/** Method to return:
	 *   1) The structure of the jacobian (if "values" is NULL)
	 *   2) The values of the jacobian (if "values" is not NULL)
	 */
	virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
				Index m, Index nele_jac, Index* iRow, Index *jCol,
				Number* values);

	/** Method to return:
	 *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
	 *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
	 */
	virtual bool eval_h(Index n, const Number* x, bool new_x,
			    Number obj_factor, Index m, const Number* lambda,
			    bool new_lambda, Index nele_hess, Index* iRow,
			    Index* jCol, Number* values);

	//@}

	/** @name Solution Methods */
	//@{
	/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
	virtual void finalize_solution(SolverReturn status,
				       Index n, const Number* x, const Number* z_L, const Number* z_U,
				       Index m, const Number* g, const Number* lambda,
				       Number obj_value,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq);
	//@}

    private:
	/**@name Methods to block default compiler methods.
	 * The compiler automatically generates the following three methods.
	 *  Since the default compiler implementation is generally not what
	 *  you want (for all but the most simple classes), we usually 
	 *  put the declarations of these methods in the private section
	 *  and never implement them. This prevents the compiler from
	 *  implementing an incorrect "default" behavior without us
	 *  knowing. (See Scott Meyers book, "Effective C++")
	 *  
	 */
	//@{
	//  LaplacianOptimiser();
	// LaplacianOptimiser(const LaplacianOptimiser&);
	// LaplacianOptimiser& operator=(const LaplacianOptimiser&);
	//@}
    };
	  
}
