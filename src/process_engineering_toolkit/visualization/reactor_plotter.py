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
        """Add vertical lines to the plot"""
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

    def plot(
        self,
        profile,
        X,
        Y,
        x_species=None,
        x_reaction=None,
        y_species=None,
        y_reaction=None,
        # species_for_conversion=None,
        x_unit=None,
        y_unit=None,
        vlines=None,
        hlines=None,
        legend_loc='lower right',
        save=False,
        filename="profile_plot",
        location=None,
        show=True,
        finalize=True,
    ):
        """
        Unified plotting method for reactor profiles.
        
        Parameters
        ----------
        profile : dict
            Profile dictionary with keys: 'volume', 'flow', 'concentration', 
            'mole_fraction', 'rate', 'conversion'
        X : str
            X-axis profile key. Single-value keys: 'volume', 'conversion'.
            Dict-valued keys: 'flow', 'concentration', 'mole_fraction', 'rate'.
        Y : str
            Y-axis profile key. Same valid keys as X.
        x_species : str, optional
            Species name for X-axis (required if X is a dict-valued key).
        x_reaction : str, optional
            Reaction name for X-axis (required if X is 'rate').
        y_species : str or list, optional
            Species name(s) for Y-axis (optional if Y is a dict-valued key; 
            all species plotted if not specified).
        y_reaction : str or list, optional
            Reaction name(s) for Y-axis (optional if Y is 'rate').
        x_unit : str, optional
            Target unit for X-axis data. If None, uses profile units.
        y_unit : str, optional
            Target unit for Y-axis data. If None, uses profile units.
        vlines : list, optional
            X-coordinates for vertical lines.
        hlines : list, optional
            Y-coordinates for horizontal lines.
        legend_loc : str, optional
            Legend location (default 'lower right').
        save : bool, optional
            If True, save plot to file.
        filename : str, optional
            Output filename without extension (default "profile_plot").
        location : str, optional
            Output directory (default "plots").
        show : bool, optional
            If True, display plot (default True).
        finalize : bool, optional
            If True, finalize the plot (default True).

        """
        # Validate profile keys
        self._validate_profile_keys(profile, X, Y)
        
        # Extract X and Y data
        x_data, x_unit_str, x_label = self._extract_data(
            profile, X, species=x_species, reaction=x_reaction, target_unit=x_unit, is_x_axis=True
        )
        y_data, y_unit_str, y_label = self._extract_data(
            profile, Y, species=y_species, reaction=y_reaction, target_unit=y_unit, is_x_axis=False
        )
        
        # Add species/reaction to X-axis label if specified (for clarity)
        if x_species is not None:
            x_label = f"{x_label} ({x_species})"
        elif x_reaction is not None:
            x_label = f"{x_label} ({x_reaction})"

        if Y == "conversion":
            y_label = f"Conversion of {profile['species_for_conversion']}"
        elif X == "conversion":
            x_label = f"Conversion of {profile['species_for_conversion']}"
        
        # Format axis labels: "Label, [unit]" if has units, else just "Label"
        x_full_label = f"{x_label}, [{x_unit_str}]" if x_unit_str and x_unit_str != 'dimensionless' else x_label
        y_full_label = f"{y_label}, [{y_unit_str}]" if y_unit_str and y_unit_str != 'dimensionless' else y_label
        
        # Create figure and axes
        fig, ax = plt.subplots()
        
        # Plot lines
        self._plot_lines(ax, x_data, y_data, x_label, y_label)
        
        # Set labels
        ax.set_xlabel(x_full_label)
        ax.set_ylabel(y_full_label)
        ax.set_title(f"{y_label} vs {x_label}")
        
        # Apply decorations
        self._apply_vlines(ax, vlines)
        self._apply_hlines(ax, hlines)
        self._apply_axis_style(ax, legend_loc=legend_loc)
        
        # Finalize
        self._finalize_plot(save, filename, location, show, finalize)

        if not finalize:
            return fig, ax # Return figure and axes for further manipulation if needed

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

    def _prepare_data(self, quantity, target_unit=None):
        """
        Extract magnitude and formatted unit string from Pint Quantity.
        """
        if not hasattr(quantity, "to"):
            raise TypeError("Expected Pint Quantity with units")

        if target_unit is not None:
            quantity = quantity.to(target_unit)

        magnitude = quantity.magnitude
        unit_str = f"{quantity.units:~}"
        unit_str = unit_str.replace(" ", "").replace("**", "^")

        return magnitude, unit_str

    def _finalize_plot(self, save, filename, location, show, finalize=True):
        """Save plot to file if requested, display, and close."""
        if save:
            filepath = self._build_path(filename, location)
            plt.savefig(filepath, dpi=300, bbox_inches="tight")

        if show and finalize: # Only show if finalize=True, allowing user to further manipulate the plot if needed
            plt.show()

        if finalize: # Only close if finalize=True, allowing user to further manipulate the plot if needed
            plt.close()

    def _build_path(self, filename, location):
        """Build file path for saving plots."""
        if location is None:
            location = "plots"

        os.makedirs(location, exist_ok=True)

        return os.path.join(location, f"{filename}.png")

    def _validate_profile_keys(self, profile, X, Y):
        """Validate that X and Y are valid profile keys."""
        valid_keys = profile.keys()
        if X not in valid_keys:
            raise ValueError(f"Invalid X key: '{X}'. Must be one of {valid_keys}")
        if Y not in valid_keys:
            raise ValueError(f"Invalid Y key: '{Y}'. Must be one of {valid_keys}")

    def _extract_data(self, profile, key, species=None, reaction=None, target_unit=None, is_x_axis=False):
        """
        Extract data from profile for a given key.
        
        Automatically detects whether the key maps to a scalar or dict, and handles
        extraction accordingly. Labels are auto-generated from the key name.
        
        Parameters
        ----------
        profile : dict
            Profile dictionary
        key : str
            Profile key to extract
        species : str or list, optional
            Species name(s) for dict-valued keys. If list, returns dict of selected species.
        reaction : str or list, optional
            Reaction name(s) for 'rate' key. If list, returns dict of selected reactions.
        target_unit : str, optional
            Target unit for conversion
        is_x_axis : bool, optional
            If True, requires specification for dict-valued keys
        
        Returns
        -------
        data : ndarray or dict
            If scalar key: ndarray of magnitudes
            If dict key with selection: dict of selected items or single array
            If dict key without selection (Y-axis): dict of all items
        unit_str : str
            Formatted unit string
        label : str
            Auto-generated label from key name
        """
        data = profile[key]
        label = key.capitalize()
        
        # Handle dict-valued keys
        if isinstance(data, dict):
            # X-axis requires specification
            if is_x_axis and species is None and reaction is None:
                raise ValueError(
                    f"X-axis key '{key}' requires specification. "
                    f"Available: {list(data.keys())}"
                )
            
            # Determine what to extract
            selection = species if species is not None else reaction
            
            # If nothing specified (Y-axis), return all items as dict
            if selection is None:
                result = {}
                unit_str = None
                for item_key in sorted(data.keys()):
                    quantity = data[item_key]
                    if target_unit is not None:
                        quantity = quantity.to(target_unit)
                    magnitude, unit = self._prepare_data(quantity)
                    result[item_key] = magnitude
                    unit_str = unit
                return result, unit_str, label
            
            # Convert single selection to list for uniform handling
            if isinstance(selection, str):
                selection = [selection]
            
            # Validate selections exist
            for item in selection:
                if item not in data:
                    raise ValueError(f"'{item}' not found in {key}: {list(data.keys())}")
            
            # Extract selected items
            result = {}
            unit_str = None
            for item in sorted(selection):
                quantity = data[item]
                if target_unit is not None:
                    quantity = quantity.to(target_unit)
                magnitude, unit = self._prepare_data(quantity)
                result[item] = magnitude
                unit_str = unit
            
            # For X-axis with single item, return as array not dict
            if is_x_axis and len(selection) == 1:
                return next(iter(result.values())), unit_str, label
            
            # For Y-axis or multiple items, return dict
            return result, unit_str, label
        
        else:
            # Scalar-valued key (volume, conversion, etc.)
            if target_unit is not None:
                data = data.to(target_unit)
            magnitude, unit_str = self._prepare_data(data)
            return magnitude, unit_str, label

    def _plot_lines(self, ax, x_data, y_data, x_label, y_label):
        """
        Plot lines on the given axes.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on.
        x_data : ndarray
            X-axis data (must be single array).
        y_data : ndarray or dict
            Y-axis data. If dict, plots multiple lines (one per key).
        x_label : str
            X-axis label.
        y_label : str
            Y-axis label.
        """
        if isinstance(y_data, dict):
            # Plot multiple lines
            for key in sorted(y_data.keys()):
                color = self._get_color(key)
                ax.plot(x_data, y_data[key], label=key, color=color)
        else:
            # Plot single line
            color = self._get_color(y_label)
            ax.plot(x_data, y_data, label=y_label, color=color)


# ### Next upgrade ###
# 
# plotter.compare_conversion([pfr1, pfr2])

