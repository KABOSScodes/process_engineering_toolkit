# units.py
from pint import UnitRegistry
import numpy as np

ureg = UnitRegistry()
Q_ = ureg.Quantity

def standardize_units(value, expected_unit): # Pint can parse strings, so allow expected_unit to be be str
    """
    Convert input to SI units.

    - If value is a Quantity → convert to expected_unit
    - If value is a float → assume already in expected_unit
    """

    if isinstance(value, Q_):
        return value.to(expected_unit).magnitude
    else:
        return value  # assume already in SI

def _attach_units(profile, units_map): # Should remove _
    def recurse(value, unit):
        if value is None:
            return None
        elif isinstance(value, dict):
            return {k: recurse(v, unit) for k, v in value.items()}
        else:
            return Q_(value, unit)

    return {
        key: recurse(value, units_map[key])
        for key, value in profile.items()
    }