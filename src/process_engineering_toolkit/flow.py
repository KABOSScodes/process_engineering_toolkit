# Stream = one physical flow with composition
class Stream:
    def __init__(self, name, total_flow, composition):
        """
        total_flow: mol/s
        composition: dict of mole fractions, must sum to 1
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
        return {sp: self.total_flow * x for sp, x in self.composition.items()}


# FlowBlock = collection of streams
class Flow:
    def __init__(self, streams=None):
        self.streams = streams or []

    def add_stream(self, stream):
        self.streams.append(stream)

    def total_flows(self):
        total = {}
        for stream in self.streams:
            for sp, f in stream.flows.items():
                total[sp] = total.get(sp, 0.0) + f
        return total
