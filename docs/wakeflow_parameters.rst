
Wakeflow Parameters
===================

Model Parameters
----------------

+-----------------------+--------------------------+-------------------------+------------------------------------------+
| Parameter             | Type                     | Default Value           | Description                              |
+=======================+==========================+=========================+==========================================+
| `name`                | `str`                    | `'default_results_dir'` | Default results directory.               |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `system`              | `str`                    | `'default_parent_dir'`  | Default parent directory of results.     |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `m_star`              | `float`                  | `1.0`                   | Mass of central star in Solar masses.    |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `m_planet`            | `float` or `List[float]` | `0.1`                   | Mass of embedded planet in Jupiter       |
|                       |                          |                         | masses. Giving a list generates a        |
|                       |                          |                         | model per planet mass.                   |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `r_outer`             | `float`                  | `500`                   | Outer radius of disk in AU.              |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Inner radius of disk in AU. Ignored      |
| `r_inner`             | `float`                  | `100`                   | for Cartesian grid.                      |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `r_planet`            | `float`                  | `250`                   | Orbital radius of planet in AU.          |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Reference radius in AU, where `hr`       |
| `r_ref`               | `float`                  | `r_planet`              | and `dens_ref` are defined.              |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Critical radius in AU for exponentially  |
|                       |                          |                         | tapered disk. Exponential taper not      |
| `r_c`                 | `float`                  | `0`                     | used if set to 0.                        |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | The sound speed is proportional to       |
| `q`                   | `float`                  | `0.25`                  | R^-q.                                    |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `p`                   | `float`                  | `1.0`                   | The density is proportional to R^-p.     |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | The disk aspect ratio at the reference   |
| `hr`                  | `float`                  | `0.1`                   | radius.                                  |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | The gas density at the reference         |
| `dens_ref`            | `float`                  | `1.0`                   | radius.                                  |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `cw_rotation`         | `bool`                   | `False`                 | Does the disk rotate clockwise?          |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | `'cartesian'` or `'cylindrical'` or      |
| `grid_type`           | `str`                    | `'cartesian'`           | `'mcfost'`.                              |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Number of grid points in the x           |
| `n_x`                 | `int`                    | `400`                   | direction                                |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Number of grid points in the y           |
| `n_y`                 | `int`                    | `400`                   | direction                                |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `n_r`                 | `int`                    | `200`                   | Number of radial grid points             |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `n_phi`               | `int`                    | `160`                   | Number of azimuthal grid points          |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Number of grid points in the z           |
| `n_z`                 | `int`                    | `50`                    | direction                                |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Use logarithmically spaced radii?        |
| `r_log`               | `int`                    | `False`                 | Ignored for `'mcfost'` grid.             |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Create plots of the perturbations        |
| `make_midplane_plots` | `bool`                   | `True`                  | in the disk mid-plane?                   |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Display plots of the perturbations       |
| `show_midplane_plots` | `bool`                   | `True`                  | in the disk mid-plane?                   |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
| `dimensionless`       | `bool`                   | `False`                 | Return dimensionless results?            |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Incude the perturbations from the        |
|                       |                          |                         | planet? Set to `False` to generate       |
| `include_planet`      | `bool`                   | `True`                  | unperturbed disk.                        |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Include the results from the linear      |
| `include_linear`      | `bool`                   | `True`                  | regime nearby the planet?                |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Save the calculated perturbations to     |
| `save_perturbations`  | `bool`                   | `True`                  | files?                                   |
+-----------------------+--------------------------+-------------------------+------------------------------------------+
|                       |                          |                         | Save the perturbations plus the          |
| `save_total`          | `bool`                   | `True`                  | unperturbed disk together to files?      |
+-----------------------+--------------------------+-------------------------+------------------------------------------+

MCFOST Interface Parameters
---------------------------

.. list-table::
   :widths: 25 15 15 100
   :header-rows: 1

   * - Parameter
     - Type
     - Default Value
     - Description
   * - `write_FITS`
     - `bool`
     - 
     - 
   * - `inclination`
     - `float`
     - 
     - 
   * - `PA`
     - `float`
     - 
     - 
   * - `PAp`
     - `float`
     - 
     - 
   * - `run_mcfost`
     - `bool`
     - 
     - 
   * - `temp_star`
     - `float`
     - 
     - 
   * - `distance`
     - `float`
     - 
     - 
   * - `v_max`
     - `float`
     - 
     - 
   * - `n_v`
     - `int`
     - 
     - 

Developer Parameters
--------------------

.. list-table::
   :widths: 25 15 15 100
   :header-rows: 1

   * - Parameter
     - Type
     - Default Value
     - Description
   * - `adiabatic_index`
     - `float`
     - 
     - 
   * - `damping_malpha`
     - `float`
     - 
     - 
   * - `CFL`
     - `float`
     - 
     - 
   * - `scale_box`
     - `float`
     - 
     - 
   * - `scale_box_ang`
     - `float`
     - 
     - 
   * - `tf_fac`
     - `float`
     - 
     - 
   * - `show_teta_debug_plots`
     - `bool`
     - 
     - 
   * - `box_warp`
     - `bool`
     - 
     - 
   * - `use_box_IC`
     - `bool`
     - 
     - 