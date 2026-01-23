from fractions import Fraction
from collections.abc import Iterable

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
        for i, reaction in enumerate(reactions):
            # Hvad hvis man vil give reaktionsparametre ved andet end pre-exponential factor og aktiveringsenergi?
            # Temp er konstant og man ved hvad k_f og k_b er direkte.
            # Måske man gerne vil give funktion for reaktionshastighed, -r_A, direkte?
            # Hvad hvis man vil reaktion ikke er elementær og man gerne vil give det (samme som ovenover)?
            # Måske det vil være en ide i stedet at bruge funktion for reaktionshastighed direkte som input?
            # Dette kan man lave en anden beregner til. For skal vi måske bare fortsætte.

            parsed_reaction = {
                f"Reaction{i+1}": reaction[0],
                "stoichiometry": self._extract_stoichiometry(reaction[0]),
                "forward_rate_parameters": reaction[1],
                "backward_rate_parameters": reaction[2] if len(reaction) > 2 else None,
                "rate_law": "mass_action",
                }

            parsed_reactions.append(parsed_reaction)
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

    def list_reactions(self):
        for reaction in self.reactions:
            print(list(reaction.values())[0])
