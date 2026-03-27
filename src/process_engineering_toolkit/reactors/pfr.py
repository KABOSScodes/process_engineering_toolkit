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
    