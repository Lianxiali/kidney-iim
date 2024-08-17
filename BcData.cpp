

#include <ibamr/config.h>

#include <ibamr/app_namespaces.h>

#include "VelocityBcCoefs.h"
#include "BcData.h"
#include "utility.h"  

#include <CartesianPatchGeometry.h>
#include <SAMRAI_config.h>


BcData::BcData(const libMesh::Mesh& mesh, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : 
    d_mesh(mesh),
    U1(input_db->getDouble("U1")),
    U2(input_db->getDouble("U2")),
    t_load(input_db->getDouble("T_LOAD")),
    tg_load(input_db->getDouble("TG_LOAD")),
    wall(input_db->getDouble("WALL")),
    d_in(input_db->getDouble("D_IN")),
    d_out(input_db->getDouble("D_OUT")),
    h1_in(input_db->getDouble("H1_IN")),
    h2_in(input_db->getDouble("H2_IN")),
    h_out(input_db->getDouble("H_OUT")),
    z_min(input_db->getDouble("Z_MIN")),
    z_max(input_db->getDouble("Z_MAX"))
{
        // Set boundary nodes based on the mesh
    pout << "d_vein_node_coordinates.size() = " << d_vein_node_coordinates.size() <<std::endl;
    if(d_vein_node_coordinates.size()==0)
        get_boundary_nodes(d_mesh, d_vein_node_coordinates, d_proximal_artery_node_coordinates, d_distal_artery_node_coordinates);
}

bool 
BcData::is_vein(const std::array<double, 3>& pt) const
{
    return is_point_inside_polygon(d_vein_node_coordinates, pt);
}

bool 
BcData::is_distal_artery(const std::array<double, 3>& pt) const
{
    return is_point_inside_polygon(d_distal_artery_node_coordinates, pt);
}
bool 
BcData::is_proximal_artery(const std::array<double, 3>& pt) const
{
    return is_point_inside_polygon(d_proximal_artery_node_coordinates, pt);
}