

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
    p_vein(input_db->getDouble("p_vein")),
    ux_proximal_artery(input_db->getDouble("ux_proximal_artery")),
    uy_proximal_artery(input_db->getDouble("uy_proximal_artery")),
    uz_proximal_artery(input_db->getDouble("uz_proximal_artery")),
    ux_distal_artery(input_db->getDouble("ux_distal_artery")),
    uy_distal_artery(input_db->getDouble("uy_distal_artery")),
    uz_distal_artery(input_db->getDouble("uz_distal_artery"))
{
    // Set boundary nodes based on the mesh
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