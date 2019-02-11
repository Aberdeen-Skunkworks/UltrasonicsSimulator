#include "pvd.hpp"



std::string pvd::header() {
    std::stringstream pvd_hdr;
    pvd_hdr << "<?xml version=\"1.0\"?>\n";
    pvd_hdr << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd_hdr << "  <Collection>\n";
    return pvd_hdr.str();
}



std::string pvd::entry(const std::string filename, const size_t entry) {
    std::stringstream pvd_entry;
    pvd_entry << "    <DataSet timestep=\"" << entry << "\" part=\"0\" file=\"" << filename << "\" />\n";
    return pvd_entry.str();
}



std::string pvd::bottom() {
    std::stringstream pvd_btm;
    pvd_btm << "  </Collection>\n";
    pvd_btm << "</VTKFile>";
    return pvd_btm.str();
}
