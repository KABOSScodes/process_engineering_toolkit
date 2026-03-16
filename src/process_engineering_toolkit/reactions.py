from collections.abc import Iterable

import numpy as np

from .rate_laws import (
    MassActionRateLaw,
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
        """
        reactant_indices = []
        reactant_orders = []

        product_indices = []
        product_orders = []

        for sp, coeff in self.stoichiometry.items():

            idx = species_index[sp]

            if coeff < 0:
                reactant_indices.append(idx)
                reactant_orders.append(abs(coeff))

            elif coeff > 0:
                product_indices.append(idx)
                product_orders.append(coeff)

        # Store on reaction
        self.reactant_indices = reactant_indices
        self.reactant_orders = reactant_orders
        self.product_indices = product_indices
        self.product_orders = product_orders

        # ALSO store on rate law
        self.rate_law.reactant_indices = reactant_indices
        self.rate_law.reactant_orders = reactant_orders
        self.rate_law.product_indices = product_indices
        self.rate_law.product_orders = product_orders
    
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
    
    def _build_rate_law(self, stoichiometry: dict[str, float], rate_law_type: str, params: dict):
        
        # Verify rate law type
        allowed_rate_laws = {"elementary", "power", "michaelis-menten"}
        if rate_law_type.lower() not in allowed_rate_laws:
            raise ValueError(f"Unknown rate law type: {rate_law_type}. Allowed types are: {allowed_rate_laws}")

        if rate_law_type.lower() == "elementary":
            # Verify parameters
            allowed_keys = {"kf", "kb", "k0_f", "k0_b", "Ea_f", "Ea_b"}
            unknown_keys = set(params.keys()) - allowed_keys
            if unknown_keys:
                raise ValueError(f"Unknown parameters for MassActionRateLaw: {unknown_keys}")
            
            # Verify forward rate specification
            has_direct_kf = "kf" in params
            has_arrhenius_f = "k0_f" in params and "Ea_f" in params
            if not (has_direct_kf or has_arrhenius_f):
                raise ValueError(
                    "Forward rate must be specified either as "
                    "'kf' or as both 'k0_f' and 'Ea_f'."
                    )

            # Verify backward rate specification
            has_kb = "kb" in params
            has_arrhenius_b = "k0_b" in params and "Ea_b" in params
            has_any_backward = any(key in params for key in ("kb", "k0_b", "Ea_b"))

            if has_any_backward:

                # Must be either kb OR Arrhenius pair
                if not (has_kb or has_arrhenius_b):
                    raise ValueError(
                        "Backward rate must be specified either as "
                        "'kb' or as both 'k0_b' and 'Ea_b'."
                        )

                # Cannot specify both forms
                if has_kb and has_arrhenius_b:
                    raise ValueError(
                        "Backward rate cannot be specified both as "
                        "'kb' and as ('k0_b', 'Ea_b')."
                        )
            
            # Construct object #### This does not yet fit into properly into framework - I think #### 
            return MassActionRateLaw(
                stoichiometry=stoichiometry,
                kf=params.get("kf"),
                kb=params.get("kb"),
                k0_f=params.get("k0_f"),
                Ea_f=params.get("Ea_f"),
                k0_b=params.get("k0_b"),
                Ea_b=params.get("Ea_b"),
                )
        
        # elif rate_law_type.lower() == "power": # To be implmented
        #     # User-defined orders
        #     return PowerLawRateLaw(
        #         orders=params.get("orders", {}),
        #         k=params.get("kf", 0)
        #     )
        
        # elif rate_law_type.lower() == "michaelis-menten": # To be implemented
        #     return MichaelisMentenRateLaw(
        #         Vmax=params.get("Vmax", 0),
        #         Km=params.get("Km", 1)
        #     )
        
        else:
            raise ValueError(f"Unknown rate law type: {rate_law_type}")

    def list_reactions(self):
        for reaction in self.reactions:
            print(reaction.name, reaction.stoichiometry)