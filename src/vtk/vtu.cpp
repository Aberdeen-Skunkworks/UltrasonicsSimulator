#include "vtu.hpp"



// void vtu::dump(const std::vector<Transducer>& transducers, const std::string vtu_filename) {
    
//     std::ofstream vtu_file;
//     std::stringstream vtus;
//     vtus << vtu_filename;
//     vtu_file.open(vtus.str());
  
//     vtu_file << header(transducers.size());
//     vtu_file << point_data(transducers);
//     vtu_file << cell_data();
//     vtu_file << points(transducers);
//     vtu_file << cells();
//     vtu_file << bottom();  
    
//     vtu_file.close();

// }



// std::string vtu::header(const size_t number_of_transducers) {
//     std::stringstream vtu_hdr;
//     vtu_hdr << "<?xml version=\"1.0\"?>\n";
//     vtu_hdr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
//     vtu_hdr << "  <UnstructuredGrid>\n";
//     vtu_hdr << "    <Piece NumberOfPoints=\"" << number_of_transducers << "\" NumberOfCells=\"0\">\n";
//     return vtu_hdr.str();    
// }



// std::string vtu::bottom() {
//     std::stringstream vtu_btm;
//     vtu_btm << "    </Piece>\n";
//     vtu_btm << "  </UnstructuredGrid>\n";
//     vtu_btm << "</VTKFile>";
//     return vtu_btm.str();
// }



// std::string vtu::point_data(const std::vector<Transducer>& transducers) {
//     std::stringstream vtu_str;
//     vtu_str << "      <PointData>\n";
//     vtu_str << "        <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"Field\" format=\"ascii\">\n          ";
//     for (size_t z = 0; z < f.N(2); z++) {
// 	for (size_t y = 0; y < f.N(1); y++) {
// 	    for (size_t x = 0; x < f.N(0); x++) {
// 		std::array<size_t, 3> c{x, y, z};
// 		vtu_str << f(c) << " ";
// 	    }
// 	}
//     }
//     vtu_str << "\n        </DataArray>\n";
//     vtu_str << "      </PointData>\n";
//     return vtu_str.str();
// }
    


// std::string vtu::cell_data() {
//     std::stringstream vtu_str;
//     vtu_str << "      <CellData>\n";
//     vtu_str << "      </CellData>\n";
//     return vtu_str.str();
// }



// std::string vtu::points(const std::vector<Transducer>& transducers) {
//     std::stringstream vtu_str;
//     vtu_str << "      <Points>\n";
//     vtu_str << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//     vtu_str << "        </DataArray>\n";
//     vtu_str << "      </Points>\n";
//     return vtu_str.str();
// }



// std::string vtu::cells() {
//     std::stringstream vtu_str;
//     vtu_str << "      <Cells>\n";
//     vtu_str << "      </Cells>\n";
//     return vtu_str.str();
// }
