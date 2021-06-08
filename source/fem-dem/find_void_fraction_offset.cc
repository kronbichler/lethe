#include <fem-dem/find_void_fraction_offset.h>

// The constructor of FindVoidFractionOffset class
template <int dim>
FindVoidFractionOffset<dim>::FindVoidFractionOffset()
{}

// Void fraction calculation using offset method
template <int dim>
void
FindVoidFractionOffset<dim>::calculate_void_fraction(
  FEValues<dim> &                       fe_values_void_fraction,
  const unsigned int &                  dofs_per_cell,
  const unsigned int &                  n_q_points,
  DoFHandler<dim> &                     void_fraction_dof_handler,
  Particles::ParticleHandler<dim, dim> &particle_handler,
  AffineConstraints<double> &           void_fraction_constraints,
  TrilinosWrappers::SparseMatrix &      system_matrix_void_fraction,
  TrilinosWrappers::MPI::Vector &       system_rhs_void_fraction)
{
  this->initialize_matrices(dofs_per_cell);
}

template class FindVoidFractionOffset<2>;
template class FindVoidFractionOffset<3>;
