from abc import ABC, abstractmethod

import numpy as np

####################### RATE LAWS #######################

class RateLaw(ABC):

    def __init__(self, stoichiometry):
        self.stoichiometry = stoichiometry
        self.reactant_indices = []
        self.reactant_orders = []
        self.product_indices = []
        self.product_orders = []
        self.R = 8.314 # J/mol/K

    @abstractmethod
    def rate(self, concentrations, T: float) -> float:
        """
        Returns reaction rate r_j [mol/(m^3·s)]
        """
        pass

class MassActionRateLaw(RateLaw):

    def __init__(
        self,
        stoichiometry,
        kf=None,
        kb=None,
        k0_f=None,
        Ea_f=None,
        k0_b=None,
        Ea_b=None,
        ):

        super().__init__(stoichiometry)

        # Store raw parameters
        self.kf = kf
        self.kb = kb
        self.k0_f = k0_f
        self.Ea_f = Ea_f
        self.k0_b = k0_b
        self.Ea_b = Ea_b

        # Determine reversibility
        self.reversible = any(x is not None for x in (kb, k0_b, Ea_b))

    def _compute_forward_k(self, T):
        if self.kf is not None:
            return self.kf
        return self.k0_f * np.exp(-self.Ea_f / (self.R * T))


    def _compute_backward_k(self, T):
        if not self.reversible:
            return 0.0

        if self.kb is not None:
            return self.kb

        return self.k0_b * np.exp(-self.Ea_b / (self.R * T))

    def rate(self, C, T):

        kf = self._compute_forward_k(T)

        rf = kf
        for idx, order in zip(self.reactant_indices, self.reactant_orders):
            rf *= C[idx] ** order

        if not self.reversible:
            return rf

        kb = self._compute_backward_k(T)

        rb = kb
        for idx, order in zip(self.product_indices, self.product_orders):
            rb *= C[idx] ** order

        return rf - rb