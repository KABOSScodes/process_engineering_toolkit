# import numpy as np
from scipy.integrate import solve_ivp

from .base import Reactor


class PFR(Reactor):

    def solve(self, V_span):

        F0 = self.build_inlet_vector()

        def ode(V, F):
            r = self.reaction_rates(F)
            dFdV = self.SM @ r
            return dFdV

        sol = solve_ivp(
            ode,
            V_span,
            F0,
            method="BDF",
            dense_output=True
        )

        self.solution = {
            "volume": sol.t,
            "flows": sol.y,
            "interpolator": sol.sol
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