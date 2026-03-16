from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import root_scalar


# Base class for all reactor models
class Reactor(ABC):

    def __init__(self, reactions, inlet_streams, parameters):
        self.reactions = reactions
        self.inlet_streams = inlet_streams
        self.parameters = parameters

        self.species = self._build_reactor_species()
        self.SM = self._build_expanded_stoichiometric_matrix()

        species_index = {sp: i for i, sp in enumerate(self.species)}

        for rxn in self.reactions.reactions:
            rxn._compile(species_index)

        self.R = 8.314 # J/(mol*K)

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
                F[i] += flows.get(sp, 0.0) # Will add proper functionality later to handle flows expressed in various units (e.g., mol/s, kg/s, etc.)
        return F

    @abstractmethod
    def solve(self, **kwargs):
        pass

    def concentrations(self, F):

        T = self.parameters["T"]
        operation = self.parameters["operation"]

        if operation == "constant_density":
            v0 = self.parameters["v0"]
            return F / v0

        elif operation == "constant_pressure":
            P = self.parameters["P"]
            F_T = np.sum(F)

            y = F / F_T
            C_T = P / (self.R * T)

            return y * C_T

        else:
            raise ValueError("Unknown operation mode")

    def reaction_rates(self, F):

        C = self.concentrations(F)
        T = self.parameters["T"]

        r = np.zeros(len(self.reactions.reactions))

        for j, rxn in enumerate(self.reactions.reactions):
            r[j] = rxn.rate(C, T)

        return r

    def equilibrium_conversion(self, species):

        species_idx = self.species.index(species)

        F0 = self.build_inlet_vector()
        FA0 = F0[species_idx]

        rxn_index = 0
        nuA = self.SM[species_idx, rxn_index]

        def residual(X):

            extent = FA0 * X / abs(nuA)

            F = F0 + self.SM[:, rxn_index] * extent

            r = self.reaction_rates(F)

            return r[rxn_index]

        sol = root_scalar(residual, bracket=[0,1], method="brentq")

        return sol.root
