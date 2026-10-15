"""Canonical rational text helpers without Python's decimal digit guard."""
from __future__ import annotations

from fractions import Fraction
import re


_CANONICAL_INTEGER = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")
_DECIMAL_LEAF_DIGITS = 256


def integer_from_decimal_text(token: str) -> int:
    if not isinstance(token, str) or not _CANONICAL_INTEGER.fullmatch(token):
        raise ValueError("integer text is not canonical decimal")
    negative = token.startswith("-")
    digits = token[1:] if negative else token
    powers = {}

    def convert(text):
        if len(text) <= _DECIMAL_LEAF_DIGITS:
            return int(text)
        split = len(text)//2
        width = len(text)-split
        if width not in powers:
            powers[width] = 10**width
        return convert(text[:split])*powers[width] + convert(text[split:])

    value = convert(digits)
    return -value if negative else value


def _integer_to_decimal_text(value: int) -> str:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("expected an integer")
    if value == 0:
        return "0"
    negative = value < 0
    value = -value if negative else value
    powers = {}

    def convert(number):
        # Leaves stay below even Python's smallest enabled digit guard (640).
        if number.bit_length() <= 512:
            return str(number)
        # This is only a balanced splitting hint, not a numeric approximation.
        # Exact divmod and zero padding reconstruct every digit unchanged.
        width = ((number.bit_length()*30103)//100000 + 1)//2
        if width not in powers:
            powers[width] = 10**width
        high, low = divmod(number, powers[width])
        return convert(high) + convert(low).zfill(width)

    text = convert(value)
    return "-" + text if negative else text


def fraction_from_rational_text(token: str, *, require_canonical: bool = False) -> Fraction:
    if not isinstance(token, str) or token.count("/") != 1:
        raise ValueError("rational text must contain one slash")
    numerator_text, denominator_text = token.split("/")
    numerator = integer_from_decimal_text(numerator_text)
    denominator = integer_from_decimal_text(denominator_text)
    if denominator <= 0:
        raise ValueError("rational denominator must be positive")
    parsed = Fraction(numerator, denominator)
    # Lexical validation above plus unchanged reduced numerator/denominator
    # is exactly the canonical decimal contract. Do not serialize arbitrarily
    # large integers merely to compare their text with the input again.
    if require_canonical and (numerator_text == '-0' or
                              parsed.numerator != numerator or
                              parsed.denominator != denominator):
        raise ValueError('authority value must be a canonical rational')
    return parsed


def fraction_to_rational_text(value: Fraction) -> str:
    if not isinstance(value, Fraction):
        raise TypeError("expected Fraction")
    return (_integer_to_decimal_text(value.numerator) + "/" +
            _integer_to_decimal_text(value.denominator))
