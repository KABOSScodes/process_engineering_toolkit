from fractions import Fraction

class Reactions:
    def __init__(self, reactions):
        self.reactions = self._parse_reactions(reactions)
    
    def _parse_reactions(self, reactions):
        parsed_reactions = []
        for i, reaction in enumerate(reactions):
            stoichiometry = self._extract_stoichiometry(reaction[0])
            params = reaction[1:]

            parsed_reaction = {
                f"Reaction {i+1}": reaction[0],
                "stoichiometry": stoichiometry,
                "forward_rate_constant": params[0],
                "backward_rate_constant": params[1] if len(params) > 1 else None,
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
