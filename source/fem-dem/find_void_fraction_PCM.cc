#include <fem-dem/find_void_fraction_PCM.h>

// The constructor of FindVoidFractionPCM class
template <int dim>
FindVoidFractionPCM<dim>::FindVoidFractionPCM()
{}

// Void fraction calculation using particle's center method
template <int dim>
void
FindVoidFractionPCM<dim>::calculate_void_fraction(
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

  for (const auto &cell : void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          this->local_matrix_void_fraction = 0;
          this->local_rhs_void_fraction    = 0;

          double particles_volume_in_cell = 0;

          // Loop over particles in cell
          // Begin and end iterator for particles in cell
          const auto pic = particle_handler.particles_in_cell(cell);
          for (auto &particle : pic)
            {
              auto particle_properties = particle.get_properties();
              particles_volume_in_cell +=
                M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
                (2 * dim);
            }
          double cell_volume = cell->measure();

          // Calculate cell void fraction
          double cell_void_fraction =
            (cell_volume - particles_volume_in_cell) / cell_volume;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  // fe_values_void_fraction.get_function_values(
                  //  nodal_void_fraction_relevant, phi_vf);
                  this->phi_vf[k] = fe_values_void_fraction.shape_value(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      this->local_matrix_void_fraction(i, j) +=
                        (this->phi_vf[j] * this->phi_vf[i]) *
                        fe_values_void_fraction.JxW(q);
                    }
                  this->local_rhs_void_fraction(i) +=
                    this->phi_vf[i] * cell_void_fraction *
                    fe_values_void_fraction.JxW(q);
                }
            }
          cell->get_dof_indices(this->local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            this->local_matrix_void_fraction,
            this->local_rhs_void_fraction,
            this->local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }
}

template class FindVoidFractionPCM<2>;
template class FindVoidFractionPCM<3>;
