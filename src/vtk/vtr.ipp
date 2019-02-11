


template<class T>
void vtr::dump(const Field<T>& f, const std::string vtr_filename) {

    std::ofstream vtr_file;
    std::stringstream vtrs;
    vtrs << vtr_filename;
    vtr_file.open(vtrs.str());
  
    vtr_file << header(f.N());
    vtr_file << point_data(f);
    vtr_file << cell_data();
    vtr_file << coordinates(f.N(), f.L(), f.orig);
    vtr_file << bottom();  
    
    vtr_file.close();
    
}


    
template<class T>
void vtr::dump_series(const std::vector<Field<T> >& farray, const std::string out_name) {

    const std::string pvd_filename = out_name + ".pvd";
    const int dir_err = mkdir(out_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    std::ofstream pvd_file;
    pvd_file.open(pvd_filename);
    pvd_file << pvd::header();
    
    for (size_t i = 0; i < farray.size(); i++) {
	
	std::stringstream vtr_filename_ss;
	vtr_filename_ss << out_name << "/" << out_name << "_" << i << ".vtr";

	pvd_file << pvd::entry(vtr_filename_ss.str(), i);

	dump(farray[i], vtr_filename_ss.str());
	
    }

    pvd_file << pvd::bottom();
    pvd_file.close();

}


    
template<class T, typename std::enable_if<std::is_base_of<Eigen::EigenBase<T>, T>::value>::type>
std::string vtr::point_data(const Field<T>& f) {
    size_t number_of_components = f({0, 0, 0}).size();
    std::stringstream vtr_str;
    vtr_str << "      <PointData>\n";
    vtr_str << "        <DataArray type=\"Float64\" NumberOfComponents=\"" << number_of_components << "\" Name=\"Field\" format=\"ascii\">\n          ";
    for (size_t z = 0; z < f.N(2); z++) {
	for (size_t y = 0; y < f.N(1); y++) {
	    for (size_t x = 0; x < f.N(0); x++) {
		std::array<size_t, 3> c{x, y, z};
		for (size_t i = 0; i < number_of_components; i++) {
		    vtr_str << f(c)[i] << " ";
		}
	    }
	}
    }
    vtr_str << "\n        </DataArray>\n";
    vtr_str << "      </PointData>\n";
    return vtr_str.str();
}
