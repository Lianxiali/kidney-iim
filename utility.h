#ifndef UTILITY_H
#define UTILITY_H

// #include "libmesh/libmesh.h"
// #include "libmesh/mesh.h"
// #include "libmesh/mesh_base.h"
// #include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/boundary_info.h"
// #include "libmesh/dof_map.h"
// #include <ibtk/ibtk_utilities.h>
// #include <ibamr/app_namespaces.h>


#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

double get_area_of_triangulation_from_point(const std::vector<std::array<double, 3>>& points, const std::array<double, 3>& center);

// Function to sort nodes in counterclockwise order (assuming 2D polygon)
void sort_nodes_ccw(std::vector<std::array<double, 3>> &points)
{
    // Calculate centroid
    double cx = 0.0, cy = 0.0;
    for (const auto &point : points)
    {
        cx += point[0];
        cy += point[1];
    }
    cx /= points.size();
    cy /= points.size();

    // Sort points by angle with respect to the centroid
    std::sort(points.begin(), points.end(), [cx, cy](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return std::atan2(a[1] - cy, a[0] - cx) < std::atan2(b[1] - cy, b[0] - cx);
    });
}

double get_area_of_polygon(const std::vector<std::array<double, 3>>& points)
{
    std::array<double, 3> center = {0.0, 0.0, 0.0};
    int n = points.size();
    
    // Calculate the centroid of the polygon
    for (const auto& point : points)
    {
        center[0] += point[0];
        center[1] += point[1];
        center[2] += point[2];
    }
    center[0] /= n;
    center[1] /= n;
    center[2] /= n;
    
    return get_area_of_triangulation_from_point(points, center);
}

double get_area_of_triangulation_from_point(const std::vector<std::array<double, 3>>& points, const std::array<double, 3>& center)
{
    double area = 0.0;
    int n = points.size();    
   
    // Calculate the area using the triangles formed by the centroid and the polygon edges
    for (int i = 0; i < n; ++i)
    {
        // Get the current and next point, wrapping around at the end
        const auto& p1 = points[i];
        const auto& p2 = points[(i + 1) % n];
        
        // Calculate the sides of the triangle
        double a = std::sqrt(std::pow(center[0] - p1[0], 2) +
                             std::pow(center[1] - p1[1], 2) +
                             std::pow(center[2] - p1[2], 2));
                             
        double b = std::sqrt(std::pow(center[0] - p2[0], 2) +
                             std::pow(center[1] - p2[1], 2) +
                             std::pow(center[2] - p2[2], 2));
                             
        double c = std::sqrt(std::pow(p1[0] - p2[0], 2) +
                             std::pow(p1[1] - p2[1], 2) +
                             std::pow(p1[2] - p2[2], 2));
                             
        // Calculate the semi-perimeter
        double s = (a + b + c) / 2.0;
        
        // Calculate the area of the triangle using Heron's formula
        area += std::sqrt(s * (s - a) * (s - b) * (s - c));
    }
    
    return area;
}

// Function to calculate the area of a polygon using the Shoelace formula
double calculate_area(const std::vector<std::array<double, 3>> &points)
{
    double area = 0.0;
    int n = points.size();
    for (int i = 0; i < n; ++i)
    {
        const auto &p1 = points[i];
        const auto &p2 = points[(i + 1) % n];
        area += p1[0] * p2[1] - p2[0] * p1[1];
    }
    return std::abs(area) / 2.0;
}

void get_boundary_nodes(
    const libMesh::Mesh& mesh0, 
    std::vector<std::array<double, 3>>& vein_node_coordinates,
    std::vector<std::array<double, 3>>& proximal_artery_node_coordinates,
    std::vector<std::array<double, 3>>& distal_artery_node_coordinates
    )
{
    libMesh::Mesh mesh = mesh0;
    libMesh::BoundaryInfo& boundary_info = mesh.get_boundary_info();

    // boundary_info.regenerate_id_sets();
    boundary_info.print_info();    
    // Added: Boundary condition
    libMesh::boundary_id_type vein_id = boundary_info.get_id_by_name("vein");
    libMesh::boundary_id_type proximal_artery_id = boundary_info.get_id_by_name("proximal_artery");
    libMesh::boundary_id_type distal_artery_id = boundary_info.get_id_by_name("distal_artery");
    boundary_info.print_summary();
    pout << "vein id: " << vein_id << std::endl;
    pout << "proximal_artery id: " << proximal_artery_id << std::endl;
    pout << "distal_artery id: " << distal_artery_id << std::endl;

    for (int i = 0; i < 6; i++)
    {
        // libmesh has default names for each side
        pout << "set " << i << " name: " << boundary_info.get_sideset_name(i) << std::endl;
    }

    // Prepare to collect boundary nodes
    std::set<libMesh::dof_id_type> vein_node_ids;
    std::set<libMesh::dof_id_type> proximal_artery_node_ids;
    std::set<libMesh::dof_id_type> distal_artery_node_ids;   
         
    // Iterate through all the elements
    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        for (auto s : elem->side_index_range())
        {
            unsigned side = s + 2;
            if (elem->neighbor_ptr(s) == nullptr) // It's a boundary side
            {
                // Select the specific boundary nodes with boundary_id = vein_id
                if (boundary_info.has_boundary_id(elem, side, vein_id)) 
                {
                    // Get the node indices on this side
                    std::vector<unsigned int> side_node_indices = elem->nodes_on_side(s);
                    for (auto node_index : side_node_indices)
                    {
                        vein_node_ids.insert(elem->node_id(node_index));
                        pout << "vein node: " << elem->node_id(node_index) << std::endl;
                    }
                }
                else if (boundary_info.has_boundary_id(elem, side, proximal_artery_id))
                {
                    std::vector<unsigned int> side_node_indices = elem->nodes_on_side(s);
                    for (auto node_index : side_node_indices)
                    {
                        proximal_artery_node_ids.insert(elem->node_id(node_index));
                        pout << "proximal artery node: " << elem->node_id(node_index) << std::endl;
                    }
                }   
                else if (boundary_info.has_boundary_id(elem, side, distal_artery_id))
                {
                    std::vector<unsigned int> side_node_indices = elem->nodes_on_side(s);
                    for (auto node_index : side_node_indices)
                    {
                        distal_artery_node_ids.insert(elem->node_id(node_index));
                        pout << "distal artery node: " << elem->node_id(node_index) << std::endl;
                    }
                }
            }
        }
    }

    for (const auto& node_id : vein_node_ids)
    {
        const libMesh::Node& node = mesh.node_ref(node_id);
        std::array<double, 3> coords = {0.0, 0.0, 0.0}; // Initialize with zeroes
        coords[0] = node(0);
        coords[1] = node(1); 
        coords[2] = node(2);    
        vein_node_coordinates.push_back(coords);        
    }

    for (const auto& node_id : proximal_artery_node_ids)
    {
        const libMesh::Node& node = mesh.node_ref(node_id);
        std::array<double, 3> coords = {0.0, 0.0, 0.0}; // Initialize with zeroes
        coords[0] = node(0);
        coords[1] = node(1); 
        coords[2] = node(2);    
        proximal_artery_node_coordinates.push_back(coords);                
    }

    for (const auto& node_id : distal_artery_node_ids)
    {
        const libMesh::Node& node = mesh.node_ref(node_id);
        std::array<double, 3> coords = {0.0, 0.0, 0.0}; // Initialize with zeroes
        coords[0] = node(0);
        coords[1] = node(1); 
        coords[2] = node(2);    
        distal_artery_node_coordinates.push_back(coords);                
    }

    std::ofstream outfile("boundary_nodes.txt");
    sort_nodes_ccw(vein_node_coordinates);
    outfile << "vein nodes " << vein_node_ids.size() << std::endl; 
    for (const auto& node : vein_node_coordinates)
    {
        outfile << node[0] << " " << node[1] << " " << node[2] << std::endl;              
    }

    outfile << "proximal artery nodes " << proximal_artery_node_ids.size() << std::endl; 
    sort_nodes_ccw(proximal_artery_node_coordinates);
    for (const auto& node : proximal_artery_node_coordinates)
    {
        outfile << node[0] << " " << node[1] << " " << node[2] << std::endl;              
    }

    outfile << "distal artery nodes " << distal_artery_node_ids.size() << std::endl; 
    sort_nodes_ccw(distal_artery_node_coordinates);
    for (const auto& node : distal_artery_node_coordinates)
    {
        outfile << node[0] << " " << node[1] << " " << node[2] << std::endl;              
    }

    pout << "distal artery boundary area: " << get_area_of_polygon(distal_artery_node_coordinates) << std::endl;      
    pout << "proximal artery boundary area: " << get_area_of_polygon(proximal_artery_node_coordinates) << std::endl;      
    pout << "vein boundary area: " << get_area_of_polygon(vein_node_coordinates) << std::endl;      

    outfile.close();
}

bool is_point_inside_polygon(const std::vector<std::array<double, 3>>& points, const std::array<double, 3>& pt)
{
    double area0 = get_area_of_polygon(points);
    double area1 = get_area_of_triangulation_from_point(points, pt);

    double rel_err = std::fabs(area0 - area1);
    return (rel_err < 1e-5);
}

#endif // UTILITY_H
