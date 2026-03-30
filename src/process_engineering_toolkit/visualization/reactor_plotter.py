import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


class ReactorPlotter:
    """
    General plotting utility for reactor objects and profiles.
    
    Parameters
    ----------
    default_style : bool, optional
        If True (default), applies a professional dark theme with science plots.
        If False, uses matplotlib defaults, leaving styling to the user.
    """

    def __init__(self, default_style=True):
        self.default_style = default_style
        
        if default_style:
            self._apply_global_style()

        # Color system
        self._color_map = {}  # key → color
        self._colors = plt.get_cmap("tab20").colors  # 20 distinct colors
        self._color_index = 0

    # ---------------------------
    # Styling
    # ---------------------------

    def _apply_global_style(self):
        """Apply global matplotlib style settings (called once in __init__)"""
        try:
            plt.style.use(['science', 'notebook', 'grid'])
        except OSError:
            # Fallback if scienceplots not available
            pass
        
        plt.rcParams.update({
            'figure.facecolor': "#2D2828",
            'savefig.facecolor': '#2E2E2E',
            'axes.titlecolor': 'white',
            'axes.labelcolor': 'white',
            'figure.figsize': (6, 4),
            'axes.labelsize': 12,
            'axes.titlesize': 18,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'lines.linewidth': 2,
            'lines.markersize': 6,
        })

    def _apply_axis_style(self, ax, legend_loc='lower right'):
        """Apply per-axis styling (tick configuration, etc)"""
        if not self.default_style:
            return
        
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        
        ax.tick_params(axis='both', which='major', length=4, width=1.2, colors='black')
        ax.tick_params(axis='both', which='minor', length=2, width=1.2, colors='black')
        ax.tick_params(axis='both', labelcolor='white')
        ax.legend(loc=legend_loc)

    def _apply_vlines(self, ax, vlines):
        """Add vertical lines to the plot using the deterministic color system"""
        if vlines is None:
            return
        
        for i, x in enumerate(vlines):
            color = self._get_color(f"vline_{i}")
            ax.axvline(x=x, linestyle='--', label=f"vline: {x:.3g}", color=color)

    def _apply_hlines(self, ax, hlines):
        """Add horizontal lines to the plot using the deterministic color system"""
        if hlines is None:
            return
        
        for i, y in enumerate(hlines):
            color = self._get_color(f"hline_{i}")
            ax.axhline(y=y, linestyle='--', label=f"hline: {y:.3g}", color=color)

    # ---------------------------
    # Public plotting methods
    # ---------------------------

    def conversion_vs_volume(
        self,
        profile,
        species_for_conversion=None,
        volume_unit=None,
        vlines=None,
        hlines=None,
        legend_loc='lower right',
        save=False,
        filename=None,
        location=None,
        show=True,
    ):

        V, V_unit = self._prepare_data(profile["volume"], volume_unit)
        X, _ = self._prepare_data(profile["conversion"])

        fig, ax = plt.subplots()

        color = self._get_color("Conversion")

        if species_for_conversion is not None:
            ylabel = f"Conversion, {species_for_conversion}"
        else:
            ylabel = "Conversion"
        ax.plot(V, X, label="Conversion", color=color)

        ax.set_xlabel(f"Volume, [{V_unit}]")
        ax.set_ylabel(ylabel)
        ax.set_title("Conversion vs Volume")
        ax.legend()
        self._apply_vlines(ax, vlines)
        self._apply_hlines(ax, hlines)
        self._apply_axis_style(ax, legend_loc=legend_loc)

        self._finalize_plot(save, filename or "conversion_vs_volume", location, show)

    def concentration_vs_volume(
        self,
        profile,
        species=None,
        volume_unit=None,
        concentration_unit=None,
        save=False,
        filename=None,
        location=None,
        show=True,
    ):

        V, V_unit = self._prepare_data(profile["volume"], volume_unit)
        conc = profile["concentration"]

        species_to_plot = sorted(self._select_species(conc, species))

        fig, ax = plt.subplots()

        for sp in species_to_plot:
            self._validate_species(sp, conc)

            C, C_unit = self._prepare_data(conc[sp], concentration_unit)
            color = self._get_color(sp)

            ax.plot(V, C, label=sp, color=color)

        ax.set_xlabel(f"Volume [{V_unit}]")
        ax.set_ylabel(f"Concentration [{C_unit}]")
        ax.set_title("Concentration vs Volume")
        ax.legend()
        
        self._apply_axis_style(ax)

        self._finalize_plot(save, filename or "concentration_vs_volume", location, show)

    def concentration_vs_conversion(
        self,
        obj,
        species=None,
        concentration_unit=None,
        save=False,
        filename=None,
        location=None,
        show=True,
    ):
        profile = self._get_profile(obj)

        X, X_unit = self._prepare_data(profile["conversion"])
        conc = profile["concentration"]

        species_to_plot = sorted(self._select_species(conc, species))

        fig, ax = plt.subplots()

        for sp in species_to_plot:
            self._validate_species(sp, conc)

            C, C_unit = self._prepare_data(conc[sp], concentration_unit)
            color = self._get_color(sp)

            ax.plot(X, C, label=sp, color=color)

        ax.set_xlabel(f"Conversion [{X_unit}]")
        ax.set_ylabel(f"Concentration [{C_unit}]")
        ax.set_title("Concentration vs Conversion")
        ax.legend()
        
        self._apply_axis_style(ax)

        self._finalize_plot(save, filename or "concentration_vs_conversion", location, show)

    # ---------------------------
    # Internal helpers
    # ---------------------------

    def _get_color(self, key):
        """
        Deterministic color assignment using sorted keys.
        Ensures same key always gets same color.
        """
        if key not in self._color_map:
            if self._color_index >= len(self._colors):
                raise ValueError("Exceeded maximum number of unique plot colors (20).")

            self._color_map[key] = self._colors[self._color_index]
            self._color_index += 1

        return self._color_map[key]

    def _get_profile(self, obj, species_for_conversion=None): # I think we should simplify and only 
        if isinstance(obj, dict):
            return obj

        if hasattr(obj, "profile"):
            return obj.profile()

        raise ValueError("Input must be a profile dict or reactor with .profile()")

    def _prepare_data(self, quantity, target_unit=None):
        if not hasattr(quantity, "to"):
            raise TypeError("Expected Pint Quantity with units")

        if target_unit is not None:
            quantity = quantity.to(target_unit)

        magnitude = quantity.magnitude
        unit_str = f"{quantity.units:~}"
        unit_str = unit_str.replace(" ", "").replace("**", "^")

        return magnitude, unit_str

    def _select_species(self, concentration_dict, species):
        if species is None:
            return list(concentration_dict.keys())

        if isinstance(species, str):
            return [species]

        return list(species)

    def _validate_species(self, species, concentration_dict):
        if species not in concentration_dict:
            raise ValueError(f"Species '{species}' not found in profile")

    def _finalize_plot(self, save, filename, location, show):
        if save:
            filepath = self._build_path(filename, location)
            plt.savefig(filepath, dpi=300, bbox_inches="tight")

        if show:
            plt.show()

        plt.close()

    def _build_path(self, filename, location):
        if location is None:
            location = "plots"

        os.makedirs(location, exist_ok=True)

        return os.path.join(location, f"{filename}.png")


# ### Next upgrade ###
# 
# plotter.compare_conversion([pfr1, pfr2])

