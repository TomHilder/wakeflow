# wakefit_interface.py
# Written by Thomas Hilder

"""
Contains the WakefitModel class, intended for use in Wakefit to generate wake models
"""

import subprocess, os
from .model_setup         import _load_config_file, _run_setup, _Parameters
from .grid                import _Grid
from .linear_perts        import _LinearPerts
from .non_linear_perts    import _NonLinearPerts
from .wakeflow            import WakeflowModel

# class for wakefit to use to generate velocity perturbations for models
class _WakefitModel(WakeflowModel):

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