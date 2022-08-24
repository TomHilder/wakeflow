# wakeflow.py
# Written by Thomas Hilder

"""
Contains the WakeflowModel class, intended for use by users to generate, configure and run models of planet wakes.
"""

import subprocess, os
from .model_setup         import _load_config_file, _run_setup, _Parameters
from .grid                import _Grid
from .linear_perts        import _LinearPerts
from .non_linear_perts    import _NonLinearPerts
#from phantom_interface   import PhantomDump

# class for use by the user to interact with Wakeflow
class WakeflowModel():
    """
    Model object allowing you to configure, generate and save planet wake model results.

    Attributes
    ----------
    model_params : dict
        Dictionary containing the parameters to be used in the model as specified by the user.
    default_params : dict
        Dictionary containing the default WakeflowModel parameters.
    """

    # instantiate class
    def __init__(self) -> None:
        """ Instantiate model
        """

        print("Model initialised.")

    # set the parameters for the model
    def configure(
        self,
        name:                 str = "default_results_dir",
        system:               str = "default_parent_dir",
        m_star:             float = 1.0,
        m_planet:           float = 0.1,
        r_outer:            float = 500,
        r_inner:            float = 100,
        r_planet:           float = 250,
        r_ref:              float = None,
        r_c:                float = 0,
        q:                  float = 0.25,
        p:                  float = 1.0,
        hr:                 float = 0.10,
        dens_ref:           float = 1.0,
        cw_rotation:         bool = False,
        grid_type:          float = "cartesian",
        n_x:                  int = 400,
        n_y:                  int = 400,
        n_r:                  int = 200,
        n_phi:                int = 160,
        n_z:                  int = 50,
        r_log:               bool = False,
        make_midplane_plots: bool = True,
        show_midplane_plots: bool = True,
        dimensionless:       bool = False,
        include_planet:      bool = True,
        include_linear:      bool = True,
        save_perturbations:  bool = True,
        save_total:          bool = True,
        write_FITS:          bool = False,
        run_mcfost:          bool = False,
        inclination:        float = -225,
        PA:                 float = 45,
        PAp:                float = 45,
        temp_star:          float = 9250,
        distance:           float = 101.5,
        v_max:              float = 3.2,
        n_v:                float = 40
    ) -> None:
        """
        Configure the model according the specifications by the user.

        Parameters
        ----------
        name : str 
            results are saved a directory called name.
        system : str
            name of parent directory of results.
        m_star : float
            central star mass in solar masses.
        m_planet : float
            planet mass in Jupiter masses. Can also give list(float) to run multiple planet masses.
        r_outer : float
            outer disk radius in au.
        r_inner : float
            inner disk radius in au.
        r_planet : float
            orbital radius of planet in au.
        r_ref : float
            reference radius r_ref in au.
        r_c : float
            critical radius r_c in au, used for exponentially tapered density profile. ignored if set to 0.
        q : float 
            q index for sound speed profile, defined as c_s \propto r^{-q}.
        p : float
            p index for density profile (NOT surface density), defined as rho \propto r^{-p} or \propto r^{-p}*\exp{-r/r_c}^{2-p}.
        hr : float 
            disk aspect ratio at r_ref.
        dens_ref : float
            density at r_ref in g/cm^3.
        cw_rotation : bool
            clockwise disk rotation (True) or anticlockwise disk rotation (False).
        grid_type : float
            grid type, either "cylindrical" or "cartesian" or "mcfost". Must choose "mcfost" to write a .FITS file.
        n_x : int
            number of grid points in x.
        n_y : int
            number of grid points in y.
        n_r : int
            number of grid points in radius.
        n_phi : int
            number of grid points in azimuth.
        n_z : int
            number of grid points in z.
        r_log : bool
            True for logarithmically spaced radii, False for linear. Ignored for "mcfost" grid type.
        make_midplane_plots : bool
            Create midplane plots of density and velocity perturbations?
        show_midplane_plots : bool
            Display midplane plots of density and velocity perturbations to user?
        dimensionless : bool
            True for dimensionless results otherwise False.
        include_planet : bool
            Include the planet-induced perturbations? Set to False to generate unperturbed disk.
        include_linear : bool
            Include the results from the linear regime?
        save_perturbations : bool
            Save the perturbations?
        save_total : bool
            Save the totals (perturbations + background disk)?
        write_FITS : bool
            Generate a .FITS file to run in MCFOST? Requires "mcfost" grid type.
        run_mcfost : bool
            Call MCFOST to generate line emission cubes on results, requires write_FITS=True.
        inclination : float
            Inclination in degrees for MCFOST.
        PA : float
            Position angle in degrees for MCFOST.
        PAp : float
            Planet position in degrees angle for MCFOST.
        temp_star : float
            Star temperature in Kelvin for MCFOST.
        distance : float
            System distance in parsecs for MCFOST.
        v_max : float
            maxmimum magnitude velocity for cube in MCFOST.
        n_v : float
            number of velocity channels for MCFOST.
        """

        # developer parameters, intentionally not easy to change as can very easily invalidate results
        # 
        adiabatic_index       = 1.6666667     # adiabatic index
        damping_malpha        = 0.0           # artificial damping NOT IMPLEMENTED
        CFL                   = 0.5           # Courant stability factor (require <0.5)
        scale_box             = 1.0           # linear box length scale factor in radial direction
        scale_box_ang         = 1.0           # linear box length scale factor in angular direction
        tf_fac                = 1.0           # scale factor for t coordinate where self-similar solution is used
        show_teta_debug_plots = False         # show (t,eta,chi) space developer plots
        box_warp              = True          # interpret y coordinate of linear regime as arc length, or truly vertical? True (default) for former
        use_box_IC            = False         # use only part of linear regime in box as initial condition for non-linear evolution

        # generate dictionary for model parameters by grabbing all local variables
        self.model_params = locals()

        # remove "self" item
        del self.model_params["self"]

        # confirmation message
        print("Model configured.")

    # set the parameters for the model by handing a .yaml file
    def configure_from_file(
        self,
        param_file: str
    ) -> None:
        """
        Configure the model by reading in .yaml file provided by user.

        Parameters
        ----------
        param_file : str 
            path to .yaml file provided by user specifying run parameters.
        """

        # generate default parameters
        self.configure()
        self.default_params = self.model_params
        del self.model_params

        # read in file to dictionary
        self.model_params = _load_config_file(param_file, self.default_params)

        # confirmation message
        print(f"Model configuration read from file: {param_file}")

    # generate the model using the configuration specified by the user
    def run(self, overwrite: bool = False) -> None:
        """
        Generate results for model, requires user to have called either configure or configure_from_file first.

        Parameters
        ----------
        overwrite : bool 
            Overwrite previous results with identical name?
        """
        
        # run setup
        try:
            params = _run_setup(self.model_params, overwrite=overwrite)
        except ArithmeticError:
            raise Exception("Model has not been configured.")

        # grab list of planet masses
        if params.m_planet is not None:
            planet_masses = [params.m_planet]
        else:
            planet_masses = params.m_planet_array

        # run wakeflow for each planet mass
        for mass_p in planet_masses:
            params.m_planet = mass_p
            print(f"\n* Creating {mass_p} Mj model:")
            self._run_wakeflow(params)

        # run mcfost for each model
        if params.run_mcfost == True:
            print("\n* Running MCFOST on each model:")
            for mass_p in planet_masses:
                print(f"{mass_p} Mj...")
                working_dir = os.getcwd()
                os.chdir(f"{params.system}/{params.name}/{mass_p}Mj/")
                subprocess.call(
                    ["mcfost", "mcfost.para", "-df", "wakeflow_model.fits", "-mol", "-freeze-out", "20", "-photodissociation", "-photodesorption"], 
                    stdout=subprocess.DEVNULL
                )
                os.chdir(working_dir)

        print("\n* Done!")

    # internal method that is called by self.run to generate the results for a specific set of parameters
    # may be called more than once if the user has specified multiple planet masses
    def _run_wakeflow(self, params: _Parameters) -> None:
        """
        Internal use method for generating the planet wake by calling other parts of Wakeflow.

        Parameters
        ----------
        params : Parameters 
            Paramaters object specfiying options for the model to be generated.
        """

        print("Generating unperturbed background disk")

        # make empty grid for unperturbed disk
        grid_background = _Grid(params)
        grid_background._make_grid()

        # fill grid with Keplerian, power law disk
        grid_background._make_keplerian_disk()

        if params.use_planet:

            print("Extracting linear perturbations nearby planet")

            # make empty grid for linear perturbations
            grid_lin_perts = _Grid(params)
            grid_lin_perts._make_grid()
            grid_lin_perts._make_empty_disk()

            # extract linear perturbations from file
            lin_perts = _LinearPerts(params)
            lin_perts._cut_box_annulus_segment()

            # add the linear perturbations onto grid
            grid_lin_perts._add_linear_perturbations(lin_perts, grid_background.rho)

            # make empty grid for non-linear perturbations
            grid_nonlin_perts = _Grid(params)
            grid_nonlin_perts._make_grid()
            grid_nonlin_perts._make_empty_disk()

            # initialise non-linear perturbations
            nonlin_perts = _NonLinearPerts(params, grid_nonlin_perts)

            # extract initial condition from the linear perturbations
            nonlin_perts._extract_ICs(lin_perts)
            if params.use_box_IC:
                nonlin_perts._extract_ICs_ann(lin_perts)

            # solve for non-linear perturbations
            nonlin_perts._get_non_linear_perts()

            # add non-linear perturbations to grid
            grid_nonlin_perts._add_non_linear_perturbations(nonlin_perts, grid_background.rho)

            # merge grids for result
            if params.include_linear:
                grid_background._merge_grids(grid_lin_perts)

            # merge grids for results
            grid_background._merge_grids(grid_nonlin_perts)

            # flip results if desired
            if params.user_cw_rotation:
                grid_background._flip_results()

            # merge grids to save or plot perturbations
            if params.make_midplane_plots or params.save_perturbations:

                if params.include_linear:
                    grid_nonlin_perts._merge_grids(grid_lin_perts)

                # flip results if desired
                if params.user_cw_rotation:
                    grid_nonlin_perts._flip_results()

                if params.dimensionless:
                    grid_nonlin_perts._remove_dimensions()

                if params.make_midplane_plots:
                    if params.show_midplane_plots:
                        print('\n* Displaying results:')
                    grid_nonlin_perts._show_disk2D(0, show=params.show_midplane_plots, save=True, dimless=params.dimensionless)

        if params.save_perturbations or params.save_total:
            print("\n* Saving results:")

        # save perturbations
        if params.save_perturbations:
            #print("Saving perturbations to file")
            grid_nonlin_perts._save_results("delta", "Perturbations")

        if params.dimensionless:
            grid_background._remove_dimensions(scale_dens=True)

        # save perts + background
        if params.save_total:
            #print("Saving background + perturbations to file")
            grid_background._save_results("total", "Total        ")

        # write fits file
        if params.write_FITS:
            grid_background._write_fits_file()