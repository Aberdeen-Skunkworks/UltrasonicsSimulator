#include "vtr.hpp"


    
std::string vtr::header(const std::array<size_t, 3> N) {
    std::stringstream vtr_hdr;
    vtr_hdr << "<?xml version=\"1.0\"?>\n";
    vtr_hdr << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\">\n";
    vtr_hdr << "  <RectilinearGrid WholeExtent=\"0 " << N[0]-1 << " 0 " << N[1]-1 << " 0 " << N[2]-1 << "\">\n";
    vtr_hdr << "    <Piece Extent=\"0 " << N[0]-1 << " 0 " << N[1]-1 << " 0 " << N[2]-1 << "\">\n";
    return vtr_hdr.str();    
}



std::string vtr::bottom() {
    std::stringstream vtr_btm;
    vtr_btm << "    </Piece>\n";
    vtr_btm << "  </RectilinearGrid>\n";
    vtr_btm << "</VTKFile>";
    return vtr_btm.str();
}


    
std::string vtr::cell_data() {
    std::stringstream vtr_str;
    vtr_str << "      <CellData>\n";
    vtr_str << "      </CellData>\n";
    return vtr_str.str();
}



std::string vtr::coordinates(const std::array<size_t, 3> N, const std::array<double, 3> L, const std::array<double, 3> origin) {
    std::stringstream vtr_str;
    vtr_str << "      <Coordinates>\n";
    for (size_t n = 0; n < 3; n++) {
	vtr_str << "        <DataArray type=\"Float64\" format=\"ascii\">\n          ";
	for (size_t r = 0; r < N[n]; r++) {
	    vtr_str << origin[n] + (r + 0.5) * L[n] / N[n] << " ";
	}
	vtr_str << "\n        </DataArray>\n";		
    }
    vtr_str << "      </Coordinates>\n";
    return vtr_str.str();
}



std::string vtr::point_data(const Field<double>& f) {
    std::stringstream vtr_str;
    vtr_str << "      <PointData>\n";
    vtr_str << "        <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"Field\" format=\"ascii\">\n          ";
    for (size_t z = 0; z < f.N(2); z++) {
	for (size_t y = 0; y < f.N(1); y++) {
	    for (size_t x = 0; x < f.N(0); x++) {
		std::array<size_t, 3> c{x, y, z};
		vtr_str << f(c) << " ";
	    }
	}
    }
    vtr_str << "\n        </DataArray>\n";
    vtr_str << "      </PointData>\n";
    return vtr_str.str();
}



std::string vtr::point_data(const Field<Eigen::Vector3d>& f) {
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


