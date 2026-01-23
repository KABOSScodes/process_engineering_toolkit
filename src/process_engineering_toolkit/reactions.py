from fractions import Fraction
from collections.abc import Iterable
from abc import ABC, abstractmethod
import numpy as np

class Reactions:
    """This class currently assumes elementary reactions"""
    def __init__(self, reactions):
        self.reaction_specs = self._normalize_input(reactions)
        self.reactions = self._parse_reactions(self.reaction_specs)
    
    def _normalize_input(self, reactions):
        """
        Accept:
          - single reaction spec
          - list/tuple of reaction specs

        Reaction spec is list or tuple of format: 
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
        Heuristic check:
        ("reaction string", forward_params (k0_f and E_f in a list), optional backward_params (Same format as forward_params))
        """
        if not isinstance(obj, (list, tuple)):
            return False

        if len(obj) not in (2, 3):
            return False

        if not isinstance(obj[0], str):
            return False

        if not isinstance(obj[1], (list, tuple)):
            return False

        if len(obj) == 3 and not isinstance(obj[2], (list, tuple)):
            return False

        return True
    
    def _parse_reactions(self, reactions):
        parsed_reactions = []

        for i, reaction_spec in enumerate(reactions):
            reaction_string = reaction_spec[0]
            params = reaction_spec[1:]

            # Parse stoichiometry
            stoichiometry = self._extract_stoichiometry(reaction_string)

            # Build rate-law object
            rate_law = self._build_rate_law(
                stoichiometry=stoichiometry,
                params=params,
                reversible="<->" in reaction_string,
            )

            # Create Reaction object
            reaction = Reaction(
                name=f"Reaction {i+1}",
                reaction_string=reaction_string,
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
            stoichiometry[species] = -Fraction(coeff)

        for i, term in enumerate(products):
            coeff, species = self._split_coeff_species(term)
            stoichiometry[species] = Fraction(coeff)

        return stoichiometry
    
    def _split_coeff_species(self, term):
        """Splits a term into coefficient and species."""
        for i, ch in enumerate(term):
            if ch.isalpha():
                coeff = term[:i] or "1"
                species = term[i:]
                return coeff, species
    
    def _build_rate_law(self, stoichiometry, params, reversible):
        """To be implemented"""

    def list_reactions(self):
        for reaction in self.reactions:
            print(list(reaction.values())[0])


class RateLaw(ABC):
    """Abstract base class for a reaction rate law."""

    @abstractmethod
    def evaluate(self, concentrations: dict[str, float], T: float = None) -> float:
        """
        Evaluate the reaction rate.
        :param concentrations: dict of species concentrations {species: C}
        :param T: Temperature in K (if needed)
        :return: reaction rate 
        """
        pass

class ElementaryRateLaw(RateLaw):
    def __init__(
        self,
        stoichiometry: dict[str, float],
        k0_f: float,
        Ea_f: float,
        k0_b: float | None = None,
        Ea_b: float | None = None,
    ):
        
        self.stoichiometry = stoichiometry
        self.k0_f = k0_f
        self.Ea_f = Ea_f
        self.k0_b = k0_b
        self.Ea_b = Ea_b

    def k_forward(self, T: float) -> float:
        R = 8.314
        return self.k0_f * np.exp(-self.Ea_f / (R * T))

    def k_backward(self, T: float) -> float:
        if self.k0_b is None:
            return 0.0
        R = 8.314
        return self.k0_b * np.exp(-self.Ea_b / (R * T))

    def evaluate(self, concentrations: dict[str, float], T: float) -> float:
        """
        Elementary mass-action rate:
        r = k_f * Π C_i^{-ν_i}  (reactants)
          - k_b * Π C_j^{ν_j}   (products)
        """
        # Forward rate
        r_f = self.k_forward(T)
        for species, coeff in self.stoichiometry.items():
            if coeff < 0:
                if species not in concentrations:
                    raise KeyError(f"Missing concentration for species '{species}'")
                r_f *= concentrations[species] ** (-coeff)

        # Backward rate (if reversible)
        r_b = 0.0
        if self.k0_b is not None:
            r_b = self.k_backward(T)
            for species, coeff in self.stoichiometry.items():
                if coeff > 0:  # products
                    if species not in concentrations:
                        raise KeyError(f"Missing concentration for species '{species}'")
                    r_b *= concentrations[species] ** coeff

        return r_f - r_b
    

class Reaction:
    def __init__(
        self,
        name: str,
        reaction_string: str,
        stoichiometry: dict[str, float],
        rate_law: RateLaw,
    ):
        self.name = name
        self.reaction_string = reaction_string
        self.stoichiometry = stoichiometry
        self.rate_law = rate_law

    def rate(self, concentrations, T):
        return self.rate_law.evaluate(concentrations, T)
