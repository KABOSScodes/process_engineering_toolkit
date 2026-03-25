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

    # def volume_for_conversion(self, species, X_target, n_points=500):
    #     """
    #     Compute the reactor volume required to reach a given conversion of `species`.
    #     """

    #     if self.solution is None:
    #         raise RuntimeError("Reactor must be solved first")

    #     # Get profile on a dense grid
    #     profile = self.profile(species_for_conversion=species, n_points=n_points)
    #     X = profile["conversion"]
    #     V = profile["volume"]

    #     # Interpolate conversion as a function of volume
    #     f = interp1d(X, V, bounds_error=True)

    #     # Compute required volume
    #     V_required = f(X_target)

    #     return float(V_required)

    # def conversion_at_volume(self, species, V_input, n_points=500):
    #     """
    #     Conversion of `species` achievable at a reactor volume `V`.
    #     """

    #     if self.solution is None:
    #         raise RuntimeError("Reactor must be solved first")

    #     # Get profile on a dense grid
    #     profile = self.profile(species_for_conversion=species, n_points=n_points)
    #     X = profile["conversion"]
    #     V = profile["volume"]

    #     # Interpolate volume as a function of conversion
    #     f = interp1d(V, X, bounds_error=True)

    #     # Compute achievable conversion
    #     X_V = f(V_input)

    #     return float(X_V)
    