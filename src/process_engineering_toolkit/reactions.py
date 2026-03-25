from collections.abc import Iterable

import numpy as np

from .rate_laws import (
    RATE_LAWS,
    RateLaw,
)


class Reaction:
    def __init__(
        self,
        name: str,
        # reaction_string: str,
        stoichiometry: dict[str, float],
        rate_law: RateLaw,
    ):
        self.name = name
        # self.reaction_string = reaction_string
        self.stoichiometry = stoichiometry
        self.rate_law = rate_law
    
    def _compile(self, species_index):
        """
        Convert species names to indices once.
        Also validates user-provided species for power-law rate laws.
        """
        # Forward term (reactants)
        if hasattr(self.rate_law, "orders_f"):
            for sp in self.rate_law.orders_f:
                if sp not in species_index:
                    raise ValueError(f"Unknown species in rate expression: {sp}")

        # Backward term (products)
        if hasattr(self.rate_law, "orders_b"):
            for sp in self.rate_law.orders_b:
                if sp not in species_index:
                    raise ValueError(f"Unknown species in rate expression: {sp}")

        self.reactant_indices = [species_index[sp] for sp in self.rate_law.orders_f]
        self.reactant_orders = list(self.rate_law.orders_f.values())

        self.product_indices = [species_index[sp] for sp in self.rate_law.orders_b]
        self.product_orders = list(self.rate_law.orders_b.values())

        # Also store indices on rate law for fast evaluation
        self.rate_law.reactant_indices = self.reactant_indices
        self.rate_law.reactant_orders = self.reactant_orders
        self.rate_law.product_indices = self.product_indices
        self.rate_law.product_orders = self.product_orders
    
    def rate(self, C, T):
        return self.rate_law.rate(C, T)

class Reactions:
    """This class currently assumes elementary reactions"""
    def __init__(self, reactions):
        self.reaction_specs = self._normalize_input(reactions)
        self.reactions = self._parse_reactions(self.reaction_specs)
        self.species = self._build_species()
        self.SM = self._build_stoichiometric_matrix()
        self.reactant_orders = self._build_reactant_orders()
        self.product_orders = self._build_product_orders()
    
    def _build_species(self):
        species_set = set()

        for rxn in self.reactions:
            for sp in rxn.stoichiometry.keys():
                species_set.add(sp)
        
        return sorted(species_set)

    def _build_stoichiometric_matrix(self):
        n_species = len(self.species)
        n_reactions = len(self.reactions)

        SM = np.zeros((n_species, n_reactions))

        for j, rxn in enumerate(self.reactions):
            stoich = rxn.stoichiometry

            for i, sp in enumerate(self.species):
                SM[i, j] = stoich.get(sp, 0.0)

        return SM

    def _build_reactant_orders(self): #!# 
        return -np.minimum(self.SM, 0)
    
    def _build_product_orders(self): #!# 
        return np.maximum(self.SM, 0)
    
    def _normalize_input(self, reactions):
        """
        Accept:
          - single reaction spec
          - list/tuple of reaction specs

        Reaction spec is list or tuple of format: ### UPDATE THIS TO DICT AS IN NEW FORMAT
        ("reaction string", forward_params (k0_f and E_f in a list), optional backward_params (k0_b and E_b in a list))
        With k0 being the pre-exponential factor, E the activation energy and _f and _b indicating forward and backward reactions.
        
        Always return: list of reaction specs
        """

        # Case 1: single reaction spec
        if self._is_reaction_spec(reactions):
            return [reactions]

        # Case 2: iterable of reaction specs
        if isinstance(reactions, Iterable):
            if all(self._is_reaction_spec(r) for r in reactions):
                return list(reactions)

        raise TypeError(
            "Reactions must be a reaction spec or a list/tuple of reaction specs"
        )

    def _is_reaction_spec(self, obj):
        """
        Validate that obj matches the reaction spec format:

        {
            "equation": str,
            "rate_law": str,
            "parameters": dict,
            optional "name": str
        }
        """
        if not isinstance(obj, dict):
            return False

        # Required keys
        required_keys = {"equation", "rate_law", "parameters"}
        if not required_keys.issubset(obj.keys()):
            return False

        # Type checks
        if not isinstance(obj["equation"], str):
            return False

        if not isinstance(obj["rate_law"], str):
            return False

        if not isinstance(obj["parameters"], dict):
            return False

        # Optional name
        if "name" in obj and not isinstance(obj["name"], str):
            return False

        return True
    
    def _parse_reactions(self, reactions):
        parsed_reactions = []

        for i, reaction_spec in enumerate(reactions):
            # Extract reaction string and stoichiometry
            reaction_string = reaction_spec["equation"]
            stoichiometry = self._extract_stoichiometry(reaction_string)

            # Build rate-law object
            rate_law = self._build_rate_law(
                stoichiometry=stoichiometry,
                rate_law_type=reaction_spec["rate_law"],
                params=reaction_spec["parameters"],
            )

            # Determine reaction name
            name = reaction_spec.get("name", f"Reaction {i+1}")

            # Create Reaction object
            reaction = Reaction(
                name=name,
                stoichiometry=stoichiometry,
                rate_law=rate_law,
            )

            parsed_reactions.append(reaction)

        return parsed_reactions
    
    def _extract_stoichiometry(self, reaction):
        """Parses a reaction string into species and stoichiometric coefficients."""
        # Split reaction into reactants and products
        if "<->" in reaction:
            reactants, products = reaction.split("<->")
        elif "->" in reaction:
            reactants, products = reaction.split("->")
        else:
            raise ValueError("Reaction must contain '->' or '<->'")

        # Split on plus signs and remove all whitespace
        reactants = [x.strip() for x in reactants.split("+")]
        products  = [x.strip() for x in products.split("+")]
        reactants = [x.replace(" ", "") for x in reactants]
        products  = [x.replace(" ", "") for x in products]

        stoichiometry = {}
        for i, term in enumerate(reactants):
            coeff, species = self._split_coeff_species(term)
            stoichiometry[species] = -coeff

        for i, term in enumerate(products):
            coeff, species = self._split_coeff_species(term)
            stoichiometry[species] = coeff

        return stoichiometry
    
    def _split_coeff_species(self, term):
        """Splits a term into coefficient and species."""
        for i, ch in enumerate(term):
            if ch.isalpha():
                coeff = int(term[:i]) if term[:i] else 1
                species = term[i:]
                return coeff, species
    
    def _build_rate_law(self, stoichiometry, rate_law_type, params):
        rate_law_cls = RATE_LAWS.get(rate_law_type.lower())
        if rate_law_cls is None:
            raise ValueError(f"Unknown rate law: {rate_law_type}")
        return rate_law_cls(stoichiometry=stoichiometry, **params)

    def list_reactions(self):
        for reaction in self.reactions:
            print(reaction.name, reaction.stoichiometry)