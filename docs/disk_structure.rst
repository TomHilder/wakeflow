
Disk Structure
==============

This section of the documentation outlines the background, unperturbed disk model used by ``wakeflow``.

Temperature
-----------

We assume that the sound speed :math:`c_s` obeys a simple radial power law:

.. math::

    c_s \propto R^{-q},

where :math:`R` is a cylindrical radius coordinate and :math:`q` is some real number. Thus temperature :math:`T` must scale

.. math::

    T \propto c_s^2 \propto R^{-2q}.

Density
-------

The density structure is derived asssuming the disk is in vertical hydrostatic equilibrium (eg. Pringle 1981). The density :math:`\rho` is given by

.. math::

    \rho(R,z) = \rho(R_\mathrm{ref}) \left( \frac{R}{R_\mathrm{ref}} \right)^{-p} \exp{\left(\frac{G M_*}{c_s^2} \left[ \frac{1}{\sqrt{R^2 + z^2}} - \frac{1}{R} \right]\right)},

where :math:`z` is the height, :math:`R_\mathrm{ref}` is some reference radius, :math:`p` is some real number, :math:`G` is the gravitational constant and :math:`M_*` is the mass of the central star. 

``wakeflow`` also supports the commonly used "exponentially-tapered" density structure, given by

.. math::

    \rho(R,z) = \rho(R_\mathrm{ref}) \left( \frac{R}{R_\mathrm{ref}} \right)^{-p} \exp{\left( -\left[ \frac{R}{R_\mathrm{c}} \right]^{2-p} \right)} \exp{\left(\frac{G M_*}{c_s^2} \left[ \frac{1}{\sqrt{R^2 + z^2}} - \frac{1}{R} \right]\right)},

where :math:`R_\mathrm{c}` is the "critical radius". This density structure is automatically used in ``wakeflow`` if you specify ``r_c`` to be non-zero.

Velocities
----------

The velocities are derived assuming radial force balance, where the additional factor that modifies the rotation from Keplerian is due to the radial pressure forces acting on the gas in the disk. `eg. Nelson et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.2610N/abstract>`_. The radial and vertical motions are set to zero, while the rotation is given by

.. math::

   \Omega(R,z) = \Omega_\mathrm{K} \left[ -(p+2q) \left( \frac{H}{R} \right)^2 + (1-2q) + \frac{2qR}{\sqrt{R^2 + z^2}}  \right]^{1/2},

where :math:`\Omega_\mathrm{K}=\sqrt{\frac{GM_*}{R^3}}` is Keplerian rotation and :math:`H=c_s/\Omega_\mathrm{K}` is the disk scale height.

Surface Density
---------------

Very commonly the density structure in disk models is specified in terms of the surface density :math:`\Sigma`, defined by

.. math::

    \Sigma = \int_{-\infty}^{\infty} \rho \, dz.

Assuming we care only about regions of the disk where :math:`z \ll R`, one can show that density and surface density are related by (see for example the lecture notes by Armitage, 2022)

.. math::

    \Sigma = \sqrt{2\pi} H \rho.

Thus if you parameterise the disk density by :math:`\Sigma \propto R^{-\gamma}`, then :math:`\gamma` is related to :math:`p` and :math:`q` by

.. math::

    p = \frac{3}{2} - q + \gamma.

Thus :math:`p` and :math:`\gamma` are not in general the same, although it is tempting to think that they could be.