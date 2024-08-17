#ifndef BCDATA_H
#define BCDATA_H

#include <libmesh/mesh.h>
#include <tbox/Database.h>

class BcData
{
public:
    BcData(const libMesh::Mesh& mesh, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    const double U1;
    const double U2;
    const double t_load;
    const double tg_load;
    const double wall;
    const double d_in;
    const double d_out;
    const double h1_in;
    const double h2_in;
    const double h_out;
    const double z_min;
    const double z_max;
    const libMesh::Mesh& d_mesh;

    std::vector<std::array<double, 3>> d_vein_node_coordinates;
    std::vector<std::array<double, 3>> d_proximal_artery_node_coordinates;
    std::vector<std::array<double, 3>> d_distal_artery_node_coordinates;  

    bool is_vein(const std::array<double, 3>& pt) const;
    bool is_distal_artery(const std::array<double, 3>& pt) const;
    bool is_proximal_artery(const std::array<double, 3>& pt) const;

};

#endif // BCDATA_H
