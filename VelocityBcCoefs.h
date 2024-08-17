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
    VelocityBcCoefs(const INSHierarchyIntegrator* fluid_solver, const BcData& bc_data, int comp_idx);

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
};

#endif // INCLUDED_VELOCITY_BC_COEFS
