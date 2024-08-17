#ifndef INCLUDED_VELOCITY_BC_COEFS
#define INCLUDED_VELOCITY_BC_COEFS

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>
#include "BcData.h"
#include <RobinBcCoefStrategy.h>
#include <ibamr/app_namespaces.h>

#include <libmesh/mesh.h>



/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class VelocityBcCoefs is an implementation of the strategy class
 * RobinBcCoefStrategy that is used to specify velocity boundary conditions that
 * are determined by a circulation model.
 */
class VelocityBcCoefs : public RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs(const INSHierarchyIntegrator* fluid_solver, const BcData& bc_data, int comp_idx, const libMesh::Mesh& mesh);

    /*!
     * \brief Destructor.
     */
    virtual ~VelocityBcCoefs();

    /*!
     * \name Implementation of RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     */
    void setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                    Pointer<ArrayData<NDIM, double>>& bcoef_data,
                    Pointer<ArrayData<NDIM, double>>& gcoef_data,
                    const Pointer<hier::Variable<NDIM>>& variable,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double fill_time) const override;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
    IntVector<NDIM> numberOfExtensionsFillable() const override;

    //\}

    double parabolic_flow(double t, double y, double r) const;
    double time_ramp(double t) const;
    bool is_vein(const std::array<double, 3>& pt) const;
    bool is_distal_artery(const std::array<double, 3>& pt) const;
    bool is_proximal_artery(const std::array<double, 3>& pt) const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VelocityBcCoefs(const VelocityBcCoefs& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs& operator=(const VelocityBcCoefs& that) = delete;

    //    const CirculationModel* const d_circ_model;
    const INSHierarchyIntegrator* const d_fluid_solver;
    const int d_comp_idx;
    const BcData d_bc_data;
    const libMesh::Mesh& d_mesh;
    std::vector<std::array<double, 3>> d_vein_node_coordinates;
    std::vector<std::array<double, 3>> d_proximal_artery_node_coordinates;
    std::vector<std::array<double, 3>> d_distal_artery_node_coordinates;      
};

#endif // INCLUDED_VELOCITY_BC_COEFS
