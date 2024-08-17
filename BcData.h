#ifndef BCDATA_H
#define BCDATA_H

#include <libmesh/mesh.h>
#include <tbox/Database.h>

class BcData
{
public:
    BcData(const libMesh::Mesh& mesh, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    const double p_vein;
    const double ux_proximal_artery;
    const double uy_proximal_artery;
    const double uz_proximal_artery;
    const double ux_distal_artery;
    const double uy_distal_artery;
    const double uz_distal_artery;
    const libMesh::Mesh& d_mesh;

    std::vector<std::array<double, 3>> d_vein_node_coordinates;
    std::vector<std::array<double, 3>> d_proximal_artery_node_coordinates;
    std::vector<std::array<double, 3>> d_distal_artery_node_coordinates;  

    bool is_vein(const std::array<double, 3>& pt) const;
    bool is_distal_artery(const std::array<double, 3>& pt) const;
    bool is_proximal_artery(const std::array<double, 3>& pt) const;

};

#endif // BCDATA_H
