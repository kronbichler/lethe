/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/time_integration_utilities.h"
#include <core/grids.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>

#include <dem/dem.h>
#include <dem/dem_properties.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

#include <iostream>

#ifndef find_void_fraction_h
#  define find_void_fraction_h

/**
 * Base interface for classes that carry out the calculation of void fraction
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class FindVoidFraction
{
public:
  /**
   * This is the base function of the calculation of void fraction.
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
    TrilinosWrappers::MPI::Vector &       system_rhs_void_fraction) = 0;

protected:
  /**
   * Carries out initialization of the classes and vectors for the calculation
   * of void fraction
   *
   * @param dofs_per_cell
   */
  void
  initialize_matrices(const unsigned int &dofs_per_cell);

  FullMatrix<double>                   local_matrix_void_fraction;
  Vector<double>                       local_rhs_void_fraction;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<double>                  phi_vf;
};

#endif /* find_void_fraction_h */
