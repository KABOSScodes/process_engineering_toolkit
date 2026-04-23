# CSTR
from scipy.optimize import root
from scipy.interpolate import interp1d
import numpy as np

from .base import Reactor


class CSTR(Reactor):

    def solve(self, V):
        """
        Solve CSTR at steady state.
        
        Parameters
        ----------
        V : float or tuple
            Single volume (float) or volume span (tuple): (V_start, V_end)
        """
        F0 = self.build_inlet_vector()
        
        # Handle both single value and span
        if isinstance(V, (tuple, list)):
            V_array = np.linspace(V[0], V[1], 100)
        else:
            V_array = np.array([V])
        
        # Solve CSTR at each volume point
        F_out_list = []
        for V_val in V_array:
            def mole_balance(F):
                r = self.reaction_rates(F)
                return F0 - F + V_val * (self.SM @ r)
            
            sol = root(mole_balance, F0)
            if not sol.success:
                raise RuntimeError(f"CSTR solver failed at V={V_val}")
            F_out_list.append(sol.x)
        
        F_out_array = np.array(F_out_list).T  # Shape: (n_species, n_volumes)
        
        # Create interpolator for consistency with PFR
        interpolator = interp1d(V_array, F_out_array, kind='cubic', bounds_error=False)
        
        self.solution = {
            "volume": V_array,
            "flows": F_out_array,
            "interpolator": interpolator
        }
        
        return self.solution
    
    def _get_profile_units(self):
        """
        Return unit structure matching the profile dictionary.
        Units are defined at the category level where possible.
        """

        units = {
            "volume": "m^3",
            "flow": "mol/s",
            "concentration": "mol/m^3",
            "mole_fraction": "dimensionless",
            "rate": "mol/(m^3*s)",
            "conversion": "dimensionless",
        }

        return units