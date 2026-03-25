from abc import ABC #, abstractmethod
import numpy as np

# rate_laws.py

# Dictionary to hold all registered rate laws
RATE_LAWS = {}

def register_rate_law(name: str):
    """Decorator to register a rate law class with a string name."""
    def wrapper(cls):
        RATE_LAWS[name.lower()] = cls
        return cls
    return wrapper

class RateLaw(ABC):
    RATE_CONSTANT_KEYS = {"kf", "kb", "k0_f", "k0_b", "Ea_f", "Ea_b"}

    def __init__(self, stoichiometry):
        self.stoichiometry = stoichiometry
        self.reactant_indices = []
        self.reactant_orders = []
        self.product_indices = []
        self.product_orders = []
        self.orders_f = {}
        self.orders_b = {}
        self.reversible = False
        self.R = 8.314

    def _validate_rate_constants(self, params):
        keys = set(params.keys())
        unknown = keys - self.RATE_CONSTANT_KEYS
        if unknown:
            raise ValueError(f"Unknown rate constant parameters: {unknown}")

        # Forward
        has_kf = params.get("kf") is not None
        has_arrhenius_f = params.get("k0_f") is not None and params.get("Ea_f") is not None
        if has_kf and has_arrhenius_f:
            raise ValueError("Provide either 'kf' or ('k0_f' and 'Ea_f'), not both")
        if not (has_kf or has_arrhenius_f):
            raise ValueError("Forward rate constant must be provided")

        # Backward
        has_kb = params.get("kb") is not None
        has_arrhenius_b = params.get("k0_b") is not None and params.get("Ea_b") is not None
        if any(k in params for k in ["kb", "k0_b", "Ea_b"]):
            if not (has_kb or has_arrhenius_b):
                raise ValueError("Backward rate constants must be provided correctly")
            if has_kb and has_arrhenius_b:
                raise ValueError("Provide either 'kb' or ('k0_b', 'Ea_b'), not both")
    
    def _compute_forward_k(self, T=None):
        if self.kf is not None:
            return self.kf
        elif self.k0_f is not None and self.Ea_f is not None:
            if T is None:
                raise ValueError("Temperature required for Arrhenius rate constant")
            return self.k0_f * np.exp(-self.Ea_f / (self.R * T))
        else:
            raise ValueError("No forward rate constant provided")

    def _compute_backward_k(self, T=None):
        if not self.reversible:
            return 0.0
        if self.kb is not None:
            return self.kb
        elif self.k0_b is not None and self.Ea_b is not None:
            if T is None:
                raise ValueError("Temperature required for backward Arrhenius rate constant")
            return self.k0_b * np.exp(-self.Ea_b / (self.R * T))
        else:
            raise ValueError("No backward rate constant provided for reversible reaction")

    def rate(self, C, T=None):
        rf = self._compute_forward_k(T)
        for idx, order in zip(self.reactant_indices, self.reactant_orders):
            rf *= C[idx] ** order
        if not self.reversible:
            return rf
        rb = self._compute_backward_k(T)
        for idx, order in zip(self.product_indices, self.product_orders):
            rb *= C[idx] ** order
        return rf - rb

@register_rate_law("mass_action")
class MassActionRateLaw(RateLaw):
    def __init__(self, stoichiometry, **params):
        super().__init__(stoichiometry)
        self._validate_rate_constants(params)
        self.reversible = any(k in params for k in ("kb", "k0_b", "Ea_b"))
        self.kf = params.get("kf")
        self.kb = params.get("kb")
        self.k0_f = params.get("k0_f")
        self.Ea_f = params.get("Ea_f")
        self.k0_b = params.get("k0_b")
        self.Ea_b = params.get("Ea_b")

        self.orders_f = {sp: -nu for sp, nu in stoichiometry.items() if nu < 0}
        self.orders_b = {sp: nu for sp, nu in stoichiometry.items() if nu > 0}

@register_rate_law("power")
class PowerLawRateLaw(RateLaw):
    def __init__(self, stoichiometry=None, **params):
        super().__init__(stoichiometry)

        allowed_keys = {"expression"} | self.RATE_CONSTANT_KEYS
        unknown = set(params.keys()) - allowed_keys

        if unknown:
            raise ValueError(f"Unknown parameters for PowerLaw: {unknown}")

        # Extract orders
        expression = params.get("expression")
        if not expression:
            raise ValueError("Power law requires 'expression'")
        self.orders_f, self.orders_b = self._parse_expression(expression)

        # Validate rate constants
        rate_constant_params = {k: v for k, v in params.items() if k != "expression"}
        self._validate_rate_constants(rate_constant_params)
        self.reversible = bool(self.orders_b)

        # Store rate constants
        self.kf = params.get("kf")
        self.kb = params.get("kb")
        self.k0_f = params.get("k0_f")
        self.Ea_f = params.get("Ea_f")
        self.k0_b = params.get("k0_b")
        self.Ea_b = params.get("Ea_b")

    def _parse_expression(self, expr: str):
        expr = expr.replace(" ", "").replace("**", "^")
        terms = expr.split("-")
        if len(terms) > 2:
            raise ValueError("Only forward and optional backward term allowed")
        return self._parse_term(terms[0]), self._parse_term(terms[1]) if len(terms) == 2 else {}

    def _parse_term(self, term: str):
        orders = {}
        for factor in term.split("*"):
            if factor.lower() in {"k", "kf", "kb"}:
                continue
            if "^" in factor:
                sp, order = factor.split("^")
                orders[sp.upper()] = float(order)
            else:
                orders[factor.upper()] = 1.0
        return orders

# @register_rate_law("michaelis_menten") #### To be implemented ####
# class MichaelisMentenRateLaw(RateLaw):
#     def __init__(self, stoichiometry, **params):
#         self.stoichiometry = stoichiometry
#         self.Vmax = params["Vmax"]
#         self.Km = params["Km"]

#     def rate(self, C, T=None):
#         idx = list(self.stoichiometry.keys())[0]  # assume single substrate
#         return self.Vmax * C[idx] / (self.Km + C[idx])