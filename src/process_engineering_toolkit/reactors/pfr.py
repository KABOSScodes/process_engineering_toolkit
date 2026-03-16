import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

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
    
    def profile(self, species_for_conversion=None, n_points=200):

        if self.solution is None:
            raise RuntimeError("Reactor must be solved first")

        V_start = self.solution["volume"][0]
        V_end = self.solution["volume"][-1]

        V = np.linspace(V_start, V_end, n_points)

        # Evaluate interpolated solution
        F = self.solution["interpolator"](V)

        # --- concentrations ---
        C = np.zeros_like(F)

        for i in range(n_points):
            C[:, i] = self.concentrations(F[:, i])

        # --- mole fractions ---
        y = F / np.sum(F, axis=0)

        # --- reaction rates ---
        r = np.zeros((len(self.reactions.reactions), n_points))

        for i in range(n_points):
            r[:, i] = self.reaction_rates(F[:, i])

        # --- conversion ---
        conversion = None

        if species_for_conversion is not None:
            idx = self.species.index(species_for_conversion)
            FA0 = F[idx, 0]
            conversion = (FA0 - F[idx]) / FA0

        rate_dict = {
            f"R{i+1}": r[i]
            for i in range(r.shape[0])
        }

        profile = {
            "volume": V,
            "flows": dict(zip(self.species, F)),
            "concentration": dict(zip(self.species, C)),
            "mole_fraction": dict(zip(self.species, y)),
            "rate": rate_dict,
            "conversion": conversion
        }

        return profile

    def volume_for_conversion(self, species, X_target, n_points=500):
        """
        Compute the reactor volume required to reach a given conversion of `species`.
        """

        if self.solution is None:
            raise RuntimeError("Reactor must be solved first")

        # Get profile on a dense grid
        profile = self.profile(species_for_conversion=species, n_points=n_points)
        X = profile["conversion"]
        V = profile["volume"]

        # Interpolate conversion as a function of volume
        f = interp1d(X, V, bounds_error=True)

        # Compute required volume
        V_required = f(X_target)

        return float(V_required)

    def conversion_at_volume(self, species, V_input, n_points=500):
        """
        Conversion of `species` achievable at a reactor volume `V`.
        """

        if self.solution is None:
            raise RuntimeError("Reactor must be solved first")

        # Get profile on a dense grid
        profile = self.profile(species_for_conversion=species, n_points=n_points)
        X = profile["conversion"]
        V = profile["volume"]

        # Interpolate volume as a function of conversion
        f = interp1d(V, X, bounds_error=True)

        # Compute achievable conversion
        X_V = f(V_input)

        return float(X_V)
    