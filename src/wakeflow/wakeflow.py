"""
wakeflow.py

Written by Thomas Hilder
Last modified 11.08.2022

Contains the WakeflowModel class, intended for use by users to generate, configure and run models of planet wakes.
"""

import subprocess, os
from .model_setup         import load_config_file, run_setup
from .grid                import Grid
from .linear_perts        import LinearPerts
from .non_linear_perts    import NonLinearPerts
#from phantom_interface   import PhantomDump

class WakeflowModel():

    def __init__(self) -> None:

        ## find default config file
        #default_config_file = pkg_resources.resource_filename('wakeflow', 'data/default_config.yaml')
        #
        ## load default parameters
        #self.default_params = load_config_file(default_config_file)

        print("Model initialised.")

    def configure(
        self,
        name:                 str = "default_results_dir",
        system:               str = "default_parent_dir",
        m_star:             float = 1.0,
        m_planet:           float = 0.1,
        r_outer:            float = 500,
        r_inner:            float = 100,
        r_planet:           float = 250,
        r_ref:              float = 250,
        r_c:                float = 0,
        q:                  float = 0.25,
        p:                  float = 1.0,
        hr:                 float = 0.10,
        dens_ref:           float = 1.0,
        cw_rotation:         bool = False,
        type:               float = "cartesian",
        n_x:                  int = 400,
        n_y:                  int = 400,
        n_r:                  int = 200,
        n_phi:                int = 160,
        n_z:                  int = 50,
        r_log:                int = False,
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

        # developer parameters
        adiabatic_index         = 1.6666667
        damping_malpha          = 0.0
        CFL                     = 0.5
        scale_box               = 1.0
        scale_box_ang           = 1.0
        tf_fac                  = 1.0
        show_teta_debug_plots   = False
        box_warp                = True
        use_box_IC              = False  

        # generate dictionary for model parameters by grabbing all local variables
        self.model_params = locals()

        # remove "self" item
        del self.model_params["self"]

        # confirmation message
        print("Model configured.")

    def configure_from_file(
        self,
        param_file: str
    ) -> None:

        # generate default parameters
        self.configure()
        self.default_params = self.model_params
        del self.model_params

        # read in file to dictionary
        self.model_params = load_config_file(param_file, self.default_params)

        # confirmation message
        print(f"Model configuration read from file: {param_file}")

#    def configure_old(
#        self, 
#        param_dict: dict = None, 
#        param_file: str = None, 
#    ) -> None:
#
#        if param_dict is not None and param_file is not None:
#            raise Exception("Please use either dictionary or file to configure, not both")
#
#        # update parameters from provided dictionary
#        if param_dict is not None:
#            self.model_params = {**self.default_params, **param_dict}
#            print(f"Model configuration updated from dictionary: {param_dict}")
#
#        # read parameters from provided file
#        elif param_file is not None:
#            self.model_params = load_config_file(param_file, self.default_params)
#            print(f"Model configuration read from file: {param_file}")
#
#        # leave parameters unchanged
#        else:
#            print("Model configuration left as default.")

    def run(self, overwrite: bool = False) -> None:
        
        # run setup
        try:
            params = run_setup(self.model_params, overwrite=overwrite)
        except ArithmeticError:
            raise Exception("Model has not been configured.")

        #try:
        #    params = run_setup(self.model_params, default_param_dict=self.default_params, overwrite=overwrite)
        #except AttributeError:
        #    self.model_params = self.default_params
        #    params = run_setup(self.model_params, default_param_dict=self.default_params, overwrite=overwrite)

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

    def _run_wakeflow(self, params):

        print("Generating unperturbed background disk")

        # make empty grid for unperturbed disk
        grid_background = Grid(params)
        grid_background.make_grid()

        # fill grid with Keplerian, power law disk
        grid_background.make_keplerian_disk()

        if params.use_planet:

            print("Extracting linear perturbations nearby planet")

            # make empty grid for linear perturbations
            grid_lin_perts = Grid(params)
            grid_lin_perts.make_grid()
            grid_lin_perts.make_empty_disk()

            # extract linear perturbations from file
            lin_perts = LinearPerts(params)
            lin_perts.cut_box_annulus_segment()

            # add the linear perturbations onto grid
            grid_lin_perts.add_linear_perturbations(lin_perts, grid_background.rho)

            # make empty grid for non-linear perturbations
            grid_nonlin_perts = Grid(params)
            grid_nonlin_perts.make_grid()
            grid_nonlin_perts.make_empty_disk()

            # initialise non-linear perturbations
            nonlin_perts = NonLinearPerts(params, grid_nonlin_perts)

            # extract initial condition from the linear perturbations
            nonlin_perts.extract_ICs(lin_perts)
            if params.use_box_IC:
                nonlin_perts.extract_ICs_ann(lin_perts)

            # solve for non-linear perturbations
            nonlin_perts.get_non_linear_perts()

            # add non-linear perturbations to grid
            grid_nonlin_perts.add_non_linear_perturbations(nonlin_perts, grid_background.rho)

            # merge grids for result
            if params.include_linear:
                grid_background.merge_grids(grid_lin_perts)

            # merge grids for results
            grid_background.merge_grids(grid_nonlin_perts)

            # flip results if desired
            if params.user_cw_rotation:
                grid_background.flip_results()

            # merge grids to save or plot perturbations
            if params.make_midplane_plots or params.save_perturbations:

                if params.include_linear:
                    grid_nonlin_perts.merge_grids(grid_lin_perts)

                # flip results if desired
                if params.user_cw_rotation:
                    grid_nonlin_perts.flip_results()

                if params.dimensionless:
                    grid_nonlin_perts.remove_dimensions()

                if params.make_midplane_plots:
                    if params.show_midplane_plots:
                        print('\n* Displaying results:')
                    grid_nonlin_perts.show_disk2D(0, show=params.show_midplane_plots, save=True, dimless=params.dimensionless)

        if params.save_perturbations or params.save_total:
            print("\n* Saving results:")

        # save perturbations
        if params.save_perturbations:
            #print("Saving perturbations to file")
            grid_nonlin_perts.save_results("delta", "Perturbations")

        if params.dimensionless:
            grid_background.remove_dimensions(scale_dens=True)

        # save perts + background
        if params.save_total:
            #print("Saving background + perturbations to file")
            grid_background.save_results("total", "Total        ")

        # write fits file
        if params.write_FITS:
            grid_background.write_fits_file()