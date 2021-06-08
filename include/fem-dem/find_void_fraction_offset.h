/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <fem-dem/find_void_fraction.h>

#ifndef find_void_fraction_offset_h
#  define find_void_fraction_offset_h

/**
 * Calculation of void fraction using offset method
 *
 * @note For more information read: "Alobaid, F. and Epple, B., 2013.
 * Improvement, validation and application of CFD/DEM model to dense
 * gas–solid flow in a fluidized bed. Particuology, 11(5), pp.514-526."
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class FindVoidFractionOffset : public FindVoidFraction<dim>
{
public:
  FindVoidFractionOffset<dim>();

  /**
   * Carries out the calculation of void fraction using offset
   * method.
   *
   * @param fe_values_void_fraction Void fraction FE values
   * @param dofs_per_cell
   * @param n_q_points
   * @param void_fraction_dof_handler Void fraction DOF handler
   * @param particle_handler
   * @param void_fraction_constraints
   * @param system_matrix_void_fraction
   * @param system_rhs_void_fraction
   */
  virtual void
  calculate_void_fraction(
    FEValues<dim> &                       fe_values_void_fraction,
    const unsigned int &                  dofs_per_cell,
    const unsigned int &                  n_q_points,
    DoFHandler<dim> &                     void_fraction_dof_handler,
    Particles::ParticleHandler<dim, dim> &particle_handler,
    AffineConstraints<double> &           void_fraction_constraints,
    TrilinosWrappers::SparseMatrix &      system_matrix_void_fraction,
    TrilinosWrappers::MPI::Vector &       system_rhs_void_fraction) override;
};

#endif /* find_void_fraction_offset_h */
