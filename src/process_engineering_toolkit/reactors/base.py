from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import root, root_scalar
from scipy.interpolate import interp1d

from ..units import standardize_units, _attach_units, ureg


# Base class for all reactor models
class Reactor(ABC):

    def __init__(self, reactions, inlet_streams, parameters):
        self.reactions = reactions
        self.inlet_streams = inlet_streams

        # Normalize streams to list
        if isinstance(inlet_streams, (list, tuple)):
            self.inlet_streams = list(inlet_streams)
        else:
            self.inlet_streams = [inlet_streams]

        p = self._validate_params(parameters)
        self.phase = p["phase"]
        # self.operation = p["operation"]
        self.T = standardize_units(p["T"], ureg.K)
        self.P = standardize_units(p["P"], ureg.Pa)
        self.v0 = standardize_units(p["v0"], ureg.m ** 3 / ureg.s)
        # self.density = p["density"]

        self.species = self._build_reactor_species()
        self.SM = self._build_expanded_stoichiometric_matrix()

        species_index = {sp: i for i, sp in enumerate(self.species)}

        for rxn in self.reactions.reactions:
            rxn._compile(species_index)

        self.R = 8.314 # J/(mol*K)
        self.solution = None
    
    def _validate_params(self, params):

        defaults = {
            "phase": None,
            "T": None,
            "P": None,
            "v0": None,
        }

        # Check for invalid keys
        for k in params:
            if k not in defaults:
                raise ValueError(f"Unknown parameter: '{k}'")

        # Merge
        p = {**defaults, **params}

        phase = p["phase"]
        T = p["T"]
        P = p["P"]
        v0 = p["v0"]

        # -- Validate phase --
        if phase not in ("gas", "liquid"):
            raise ValueError(f"Invalid phase '{phase}'. Must be 'gas' or 'liquid'.")

        # Gas phase
        if phase == "gas":

            if T is None:
                raise ValueError("Gas-phase reactor requires temperature (T).")

            if P is None:
                raise ValueError("Gas-phase reactor requires pressure (P).")

            if v0 is not None:
                raise ValueError("Gas-phase should not define v0 (it is computed).")

        # Liquid phase
        elif phase == "liquid":

            if v0 is None:
                raise ValueError("Liquid-phase reactor requires volumetric flow 'v0'.")

            # T is optional

        return p

    def _build_reactor_species(self):
        species_set = set()

        # From reactions
        species_set.update(self.reactions.species)

        # From inlet streams
        for stream in self.inlet_streams:
            species_set.update(stream.composition.keys())

        return sorted(species_set)
    
    def _build_expanded_stoichiometric_matrix(self): 
        n_species = len(self.species)
        n_reactions = self.reactions.SM.shape[1]

        SM_expanded = np.zeros((n_species, n_reactions))

        for i, sp in enumerate(self.species):
            if sp in self.reactions.species:
                original_index = self.reactions.species.index(sp)
                SM_expanded[i, :] = self.reactions.SM[original_index, :]
            else:
                # inert species → zero row
                SM_expanded[i, :] = 0.0

        return SM_expanded

    def build_inlet_vector(self): 
        F = np.zeros(len(self.species))
        for stream in self.inlet_streams:
            flows = stream.flows
            for i, sp in enumerate(self.species):
                F[i] += flows.get(sp, 0.0)
        return F

    @abstractmethod
    def solve(self, **kwargs):
        pass

    def concentrations(self, F):

        v = self.volumetric_flow(F)

        return F / v
    
    def volumetric_flow(self, F):

        if self.phase == "gas":
            F_T = np.sum(F)
            return F_T * self.R * self.T / self.P

        elif self.phase == "liquid":
            return self.v0  # constant

        else:
            raise ValueError("Unknown phase")

    def reaction_rates(self, F):

        C = self.concentrations(F)
        T = self.T

        r = np.zeros(len(self.reactions.reactions))

        for j, rxn in enumerate(self.reactions.reactions):
            r[j] = rxn.rate(C, T)

        return r

    def equilibrium_conversion(self, species):

        F_eq = self.equilibrium_state()

        F0 = self.build_inlet_vector()
        idx = self.species.index(species)

        FA0 = F0[idx]
        if FA0 == 0:
            raise ValueError(f"Inlet of species {species} is 0.")

        return (FA0 - F_eq[idx]) / FA0

    def equilibrium_state(self, initial_guess=None, tol=1e-8):

        F0 = self.build_inlet_vector()
        n_rxn = self.SM.shape[1]

        # --- Scaling factor (based on inlet magnitude) ---
        scale = np.maximum(np.abs(F0).max(), 1.0)

        # --- Residual in scaled variables ---
        def residual_scaled(extents_scaled):
            extents = extents_scaled * scale
            F = F0 + self.SM @ extents

            # Penalize non-physical states
            if np.any(F < -1e-12):
                return np.ones(n_rxn) * 1e6

            return self.reaction_rates(F)

        # --- Initial guesses ---
        if initial_guess is not None:
            guesses = [np.array(initial_guess) / scale]
        else:
            guesses = [
                np.zeros(n_rxn),
                np.ones(n_rxn) * 1e-6,
                np.ones(n_rxn) * 0.1,
            ]

        # --- Solve attempts ---
        for guess in guesses:

            sol = root(residual_scaled, guess)

            if not sol.success:
                continue

            extents = sol.x * scale
            F_eq = F0 + self.SM @ extents

            # --- Physical check ---
            if np.any(F_eq < -1e-8):
                continue

            # --- Residual check ---
            res = residual_scaled(sol.x)
            if np.linalg.norm(res) > tol:
                continue

            return F_eq
        
        # If all attempts fail
        raise RuntimeError(
            "Equilibrium solver failed to converge to a physical solution."
        )
    
    @abstractmethod
    def _get_profile_units(self, **kwargs):
        pass

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

        if species_for_conversion is None:
            raise ValueError("Must specify species_for_conversion to compute conversion profile.")
        idx = self.species.index(species_for_conversion)
        FA0 = F[idx, 0]
        conversion = (FA0 - F[idx]) / FA0

        rate_dict = {
            f"R{i+1}": r[i]
            for i in range(r.shape[0])
        }

        # Gather profiles
        profile = {
            "species_for_conversion": species_for_conversion,
            "volume": V,
            "flow": dict(zip(self.species, F)),
            "concentration": dict(zip(self.species, C)),
            "mole_fraction": dict(zip(self.species, y)),
            "rate": rate_dict,
            "conversion": conversion
        }

        # Attach units to profile:
        units = self._get_profile_units()
        profile = _attach_units(profile, units)

        return profile

    def volume_for_conversion(self, species, X_target, n_points=500):
        """
        Compute the reactor volume required to reach a given conversion of `species`.
        """

        if self.solution is None:
            raise RuntimeError("Reactor must be solved first")

        # Get profile on a dense grid
        profile = self.profile(species_for_conversion=species, n_points=n_points) # Store profile as attribute and use instead
        X = profile["conversion"].magnitude
        V = profile["volume"].magnitude

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
        profile = self.profile(species_for_conversion=species, n_points=n_points) # Store profile as attribute and use instead
        X = profile["conversion"].magnitude
        V = profile["volume"].magnitude

        # Interpolate volume as a function of conversion
        f = interp1d(V, X, bounds_error=True)

        # Compute achievable conversion
        X_V = f(V_input)

        return float(X_V)