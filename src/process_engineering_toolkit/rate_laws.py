from abc import ABC, abstractmethod

import numpy as np

####################### RATE LAWS #######################

class RateLaw(ABC):

    def __init__(self, stoichiometry, species_list):
        self.stoichiometry = stoichiometry
        self.species_index = {s: i for i, s in enumerate(species_list)}
        self.species_list = species_list

        self.reactant_indices = []
        self.reactant_orders = []
        self.product_indices = []
        self.product_orders = []

        self.R = 8.314 # J/mol/K
        self.RATE_CONSTANT_KEYS = {"kf", "kb", "k0_f", "k0_b", "Ea_f", "Ea_b"}

    def _validate_rate_constant_params(self, params):

        keys = set(params.keys())
        allowed = self.RATE_CONSTANT_KEYS

        unknown = keys - allowed 
        if unknown:
            raise ValueError(f"Unknown parameters: {unknown}")

        # --- Forward ---
        has_kf = "kf" in params
        has_k0_f = "k0_f" in params
        has_Ea_f = "Ea_f" in params

        has_arrhenius_f = has_k0_f and has_Ea_f

        if has_kf and has_arrhenius_f:
            raise ValueError("Provide either 'kf' or ('k0_f' and 'Ea_f'), not both")

        if has_k0_f != has_Ea_f:
            raise ValueError("Forward Arrhenius requires both 'k0_f' and 'Ea_f'")

        if not (has_kf or has_arrhenius_f):
            raise ValueError("Must provide 'kf' or ('k0_f' and 'Ea_f')")

        # --- Backward ---
        has_kb = "kb" in params
        has_k0_b = "k0_b" in params
        has_Ea_b = "Ea_b" in params

        has_arrhenius_b = has_k0_b and has_Ea_b

        if has_kb and has_arrhenius_b:
            raise ValueError("Provide either 'kb' or ('k0_b' and 'Ea_b'), not both")

        if has_k0_b != has_Ea_b:
            raise ValueError("Backward Arrhenius requires both 'k0_b' and 'Ea_b'")

    @abstractmethod
    def _build(self):
        pass

    @abstractmethod
    def rate(self, concentrations, T: float) -> float:
        """
        Returns reaction rate
        """
        pass

class MassActionRateLaw(RateLaw):

    def __init__(self, stoichiometry, species_list, **params):
        super().__init__(stoichiometry, species_list)

        self._validate_rate_constant_params(params)
        self.reversible = any(k in params for k in ("kb", "k0_b", "Ea_b"))

        # Store raw parameters
        self.kf = params.get("kf")
        self.kb = params.get("kb")
        self.k0_f = params.get("k0_f")
        self.Ea_f = params.get("Ea_f")
        self.k0_b = params.get("k0_b")
        self.Ea_b = params.get("Ea_b")

        # Build
        self._build()
    
    def _build(self):
        for species, nu in self.stoichiometry.items():
            idx = self.species_index[species]
            if nu < 0:  # reactant
                self.reactant_indices.append(idx)
                self.reactant_orders.append(-nu)
            elif nu > 0:  # product
                self.product_indices.append(idx)
                self.product_orders.append(nu)

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

class PowerLawRateLaw(RateLaw):

    def __init__(self, stoichiometry, species_list, **params): 
        super().__init__(stoichiometry, species_list)

        # Get orders
        expression = params.get("expression")
        if not expression:
            raise ValueError("Power law requires 'expression'")
        self.orders_f, self.orders_b = self._parse_power_expression(expression)

        rate_constant_params = {k: v for k, v in params.items() if k in self.RATE_CONSTANT_KEYS}
        if not self.orders_b and any(k in rate_constant_params for k in ("kb", "k0_b", "Ea_b")):
            raise ValueError("Backward rate constants provided but expression has no backward term")
        
        self._validate_rate_constant_params(rate_constant_params)
        self.reversible = bool(self.orders_b)

        # Store raw parameters
        self.kf = params.get("kf")
        self.kb = params.get("kb")
        self.k0_f = params.get("k0_f")
        self.Ea_f = params.get("Ea_f")
        self.k0_b = params.get("k0_b")
        self.Ea_b = params.get("Ea_b")

        self._build()

    def _parse_power_expression(self, expression: str) -> tuple[dict, dict]:
        """
        Parses a power-law rate expression into forward and backward orders.

        Returns:
            orders_f: dict[str, float]
            orders_b: dict[str, float]
        """

        expr = expression.lower().replace(" ", "")
        expr = expr.replace("**", "^")

        terms = expr.split("-")
        if not terms[0]:
            raise ValueError("Forward term cannot be empty")

        if len(terms) > 2:
            raise ValueError("Only one forward term and one backward term allowed in power-law expression." \
            "\nMore than two terms total were found.")
        
        orders_f = self._parse_term(terms[0])
        orders_b = self._parse_term(terms[1]) if len(terms) == 2 else {}

        self._validate_species(orders_f)
        self._validate_species(orders_b)

        return orders_f, orders_b

    def _parse_term(self, term):
        orders = {}

        factors = term.split("*")
        for factor in factors:
            if factor in ["k", "kf", "k_f", "kb", "k_b"]:
                continue
            if factor.replace(".", "").isdigit():
                raise ValueError("Numeric pre-factors not allowed in power-law expression.")
            if "^" in factor:
                base, exponent = factor.split("^")
                # Clean exponent
                try:
                    exponent = exponent.strip("()")
                    order = float(exponent)
                except ValueError:
                    raise ValueError(f"Invalid exponent in factor {factor}")
            else:
                base = factor
                order = 1.0
            
            if base.startswith("c_"):
                base = base[2:]
            species = base.upper()

            if not species.isalpha():
                raise ValueError(f"Invalid species name: {species}")
            
            orders[species] = orders.get(species, 0.0) + order # 
        
        return orders
            
    def _validate_species(self, orders): 
        for species in orders:
            if species not in self.species_list:
                raise ValueError(f"Unknown species in rate expression {species}")
    
    def _build(self):

        # Forward term
        for species, order in self.orders_f.items():
            idx = self.species_index[species]
            self.reactant_indices.append(idx)
            self.reactant_orders.append(order)

        # Backward term
        for species, order in self.orders_b.items():
            idx = self.species_index[species]
            self.product_indices.append(idx)
            self.product_orders.append(order)

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