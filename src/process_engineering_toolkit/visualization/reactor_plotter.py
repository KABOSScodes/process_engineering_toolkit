import os
import matplotlib.pyplot as plt


class ReactorPlotter:
    """
    General plotting utility for reactor objects and profiles.

    Design principles:
    - Accepts either reactor objects (with .profile()) or profile dicts
    - Centralizes styling and saving
    - Easy to extend with new plot types
    """

    def __init__(self, style="default"):
        self._set_style(style)

    ### Public plotting methods

    def conversion_vs_volume(self, obj, save=False, filename=None, location=None, show=True):
        profile = self._get_profile(obj)

        plt.figure()
        plt.plot(profile["volume"], profile["conversion"], label="Conversion")

        plt.xlabel("Volume [m³]")
        plt.ylabel("Conversion")
        plt.legend()

        self._finalize_plot(save, filename or "conversion_vs_volume", location, show)

    def concentration_vs_volume(self, obj, species=None, save=False, filename=None, location=None, show=True):
        profile = self._get_profile(obj)

        plt.figure()

        conc = profile["concentration"]

        species_to_plot = self._select_species(conc, species)

        for sp in species_to_plot:
            plt.plot(profile["volume"], conc[sp], label=sp)

        plt.xlabel("Volume [m³]")
        plt.ylabel("Concentration")
        plt.legend()

        self._finalize_plot(save, filename or "concentration_vs_volume", location, show)

    def concentration_vs_conversion(self, obj, species=None, save=False, filename=None, location=None, show=True):
        profile = self._get_profile(obj)

        plt.figure()

        conc = profile["concentration"]
        species_to_plot = self._select_species(conc, species)

        for sp in species_to_plot:
            plt.plot(profile["conversion"], conc[sp], label=sp)

        plt.xlabel("Conversion")
        plt.ylabel("Concentration")
        plt.legend()

        self._finalize_plot(save, filename or "concentration_vs_conversion", location, show)

    ### Internal helpers

    def _get_profile(self, obj):
        """
        Accepts:
        - profile dict
        - reactor object with .profile()
        """
        if isinstance(obj, dict):
            return obj

        if hasattr(obj, "profile"):
            return obj.profile()

        raise ValueError("Input must be a profile dict or reactor with .profile()")

    def _select_species(self, concentration_dict, species):
        """
        Handles species selection logic
        """
        if species is None:
            return concentration_dict.keys()

        if isinstance(species, str):
            return [species]

        return species

    def _finalize_plot(self, save, filename, location, show):
        """
        Handles saving + showing
        """
        if save:
            filepath = self._build_path(filename, location)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')

        if show:
            plt.show()

        plt.close()

    def _build_path(self, filename, location):
        """
        Handles directory creation and path building
        """
        if location is None:
            location = "plots"

        os.makedirs(location, exist_ok=True)

        return os.path.join(location, f"{filename}.png")

    def _set_style(self, style):
        """
        Central styling control
        """
        if style == "default":
            plt.rcParams.update({
                "figure.figsize": (6, 4),
                "axes.labelsize": 12,
                "axes.titlesize": 14,
                "xtick.labelsize": 10,
                "ytick.labelsize": 10,
                "legend.fontsize": 10,
                "lines.linewidth": 2,
            })

        elif style == "readme":
            plt.rcParams.update({
                "figure.figsize": (6, 4),
                "axes.labelsize": 11,
                "legend.fontsize": 9,
                "lines.linewidth": 2,
            })


# ### Species handling is flexible ###
# plotter.concentration_vs_volume(pfr)              # all species
# plotter.concentration_vs_volume(pfr, "A")         # one species
# plotter.concentration_vs_volume(pfr, ["A", "B"])  # subset

# ### For README ###
# from pet.visualization import ReactorPlotter

# plotter = ReactorPlotter(style="readme")

# plotter.conversion_vs_volume(
#     pfr,
#     save=True,
#     location="plots/readme_plots"
# )

# plotter.concentration_vs_volume(pfr)

# ### Next upgrade ###
# 
# plotter.compare_conversion([pfr1, pfr2])

