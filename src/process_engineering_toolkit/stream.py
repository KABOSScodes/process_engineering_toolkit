class Stream:
    def __init__(self, name, total_flow, composition):
        """
        name: Str - stream name
        total_flow: float - total flow rate in mol/min
        composition: dict of species and mole fractions in stream - must sum to 1
        """
        self.name = name
        self.total_flow = total_flow
        if any(x < 0 for x in composition.values()):
            raise ValueError("Mole fractions must be non-negative")
        if abs(sum(composition.values()) - 1.0) > 1e-6:
            raise ValueError("Composition must sum to 1")
        self.composition = composition

    @property
    def flows(self):
        return {species: self.total_flow * molar_fraction for species, molar_fraction in self.composition.items()}