import subprocess, os, pkg_resources
from .model_setup         import load_config_file, run_setup
from .grid                import Grid
from .linear_perts        import LinearPerts
from .non_linear_perts    import NonLinearPerts
#from phantom_interface   import PhantomDump

class WakeflowModel():

    def __init__(self):

        # find default config file
        default_config_file = pkg_resources.resource_filename('wakeflow', 'data/default_config.yaml')

        # load default parameters
        self.default_params = load_config_file(default_config_file)
        print("Model initialised.")

    def configure(self, param_dict: dict = None, param_file: str = None) -> None:

        if param_dict is not None and param_file is not None:
            raise Exception("Please use either dictionary or file to configure, not both")

        # update parameters from provided dictionary
        if param_dict is not None:
            self.model_params = {**self.default_params, **param_dict}
            print(f"Model configuration updated from dictionary: {param_dict}")

        # read parameters from provided file
        elif param_file is not None:
            self.model_params = load_config_file(param_file, self.default_params)
            print(f"Model configuration read from file: {param_file}")

        # leave parameters unchanged
        else:
            print("Model configuration left as default.")

    def run(self, overwrite: bool = False) -> None:
        
        # run setup
        try:
            params = run_setup(self.model_params, default_param_dict=self.default_params, overwrite=overwrite)
        except AttributeError:
            self.model_params = self.default_params
            params = run_setup(self.model_params, default_param_dict=self.default_params, overwrite=overwrite)

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