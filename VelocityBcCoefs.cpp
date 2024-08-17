// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "VelocityBcCoefs.h"

#include <CartesianPatchGeometry.h>

bool is_point_inside_polygon(const std::vector<std::array<double, 3>>& points, const std::array<double, 3>& pt);

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace libMesh;

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs::VelocityBcCoefs
(
    const INSHierarchyIntegrator* fluid_solver, 
    const BcData& bc_data, 
    int comp_idx
)
: 
d_fluid_solver(fluid_solver), 
d_comp_idx(comp_idx), 
d_bc_data(bc_data)
{

} // VelocityBcCoefs

VelocityBcCoefs::~VelocityBcCoefs() = default; // Use default destructor

void
VelocityBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                            Pointer<ArrayData<NDIM, double>>& bcoef_data,
                            Pointer<ArrayData<NDIM, double>>& gcoef_data,
                            const Pointer<Variable<NDIM>>&  variable,
                            const Patch<NDIM>& patch,
                            const BoundaryBox<NDIM>& bdry_box,
                            double fill_time) const
{
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
pout << " I'm doing setBCCoefs " << std::endl;
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data.isNull());
#endif

    const Box<NDIM>& bc_coef_box = acoef_data->getBox();

#if !defined(NDEBUG)
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif

    const Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();

    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const hier::Index<NDIM>& i = bc();
        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);

        {
            if (d_comp_idx != axis)
            {
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }

            {
                double X[NDIM];
                std::array<double, 3> pt = {0.0, 0.0, 0.0};

                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                    pt[d] = X[d];
                }
                pout << "pt = " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
                if (d_bc_data.is_vein(pt))
                {
                    pout << "X is inside the vein boundary. " << std::endl;
                }
                else if (d_bc_data.is_proximal_artery(pt))
                {
                    pout << "X is inside the proxmial artery boundary. " << std::endl;
                }
                else if (d_bc_data.is_distal_artery(pt))
                {
                    pout << "X is inside the distal artery boundary. " << std::endl;
                }
                else
                {
                    a = 1.0;
                    b = 0.0;
                    g = 0.0;                    
                }
            }
        }

        // z_hi face (distal artery)
        if (location_index == 5)
        {
            if (d_comp_idx != axis)
            {
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }

        }
    }
    return;
} // setBcCoefs



IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable

double
VelocityBcCoefs::parabolic_flow(double t, double y, double r) const
{
    double c = (y > 0.0) ? d_bc_data.U1 : d_bc_data.U2;
    double U = c * (1.0 - r * r);
    return time_ramp(t) * U;
}

double
VelocityBcCoefs::time_ramp(double t) const
{
    return (t < d_bc_data.t_load) ? 0.0 : 1.0;
}


/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

////////////////////////////
