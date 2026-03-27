# units.py
from pint import UnitRegistry

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