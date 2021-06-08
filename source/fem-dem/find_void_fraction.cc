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

#include <fem-dem/find_void_fraction.h>

template <int dim>
void
FindVoidFraction<dim>::initialize_matrices(const unsigned int &dofs_per_cell)
{
  local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  local_rhs_void_fraction(dofs_per_cell);
  local_dof_indices.reserve(dofs_per_cell);
  phi_vf.reserve(dofs_per_cell);
}


template class FindVoidFraction<2>;
template class FindVoidFraction<3>;
