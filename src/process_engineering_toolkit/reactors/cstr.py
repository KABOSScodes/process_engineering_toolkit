# CSTR
import numpy as np
from scipy.optimize import fsolve

class CSTR:
    def __init__(self, volume, reactions, inlet_flow):
        """
        volume: reactor volume in m3
        reactions: list of Reaction objects
        inlet_flow: Flow object (aggregated inlet streams)
        """
        self.volume = volume
        self.reactions = reactions
        self.inlet_flow = inlet_flow
        
        # Species dictionary: name -> index
        self.species_names = list(self.inlet_flow.species)
        self.species_index = {name: i for i, name in enumerate(self.species_names)}
        self.num_species = len(self.species_names)
        self.num_reactions = len(reactions)
        
        # Stoichiometric matrix: nu[i,j] = stoichiometry of species i in reaction j
        self.nu = np.zeros((self.num_species, self.num_reactions))
        for j, reaction in enumerate(reactions):
            for sp, coeff in reaction.stoichiometry.items():
                if sp in self.species_index:
                    i = self.species_index[sp]
                    self.nu[i, j] = coeff

    def _reaction_rates(self, concentrations):
        """
        Calculate rate of each reaction given species concentrations.
        concentrations: np.array of species concentrations in mol/m3
        returns: np.array of reaction rates [r1, r2, ...] in mol/(m3*s)
        """
        rates = np.zeros(self.num_reactions)
        for j, reaction in enumerate(self.reactions):
            rates[j] = reaction.rate(concentrations, self.species_index)
        return rates

    def _mass_balance_equations(self, outlet_conc):
        """
        Steady-state mass balance for CSTR: F_in - F_out + V*r = 0
        outlet_conc: np.array of species concentrations (mol/m3)
        returns: np.array of residuals
        """
        F_in = np.array([self.inlet_flow.total_flow(sp) for sp in self.species_names])
        F_out = outlet_conc  # assuming volumetric flow = 1 m3/s for simplicity, can scale later
        r = self._reaction_rates(outlet_conc)
        residuals = F_in - F_out + self.volume * self.nu @ r
        return residuals

    def compute_outlet(self, guess=None):
        """
        Solve steady-state CSTR for outlet concentrations
        guess: optional initial guess for concentrations
        returns: dict {species_name: outlet_concentration}
        """
        if guess is None:
            guess = np.array([self.inlet_flow.total_flow(sp) for sp in self.species_names])
        
        outlet_conc = fsolve(self._mass_balance_equations, guess)
        return dict(zip(self.species_names, outlet_conc))

