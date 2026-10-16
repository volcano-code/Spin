"""Regression for native-owned NumPy index arrays. No native-call mocks."""
import numpy as np
import pytest
from make_it_stand import spin

FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))

@pytest.mark.parametrize('dtype', [np.int8, np.int16, np.int32, np.int64,
                                  np.uint8, np.uint16, np.uint32, np.uint64])
def test_accepts_lossless_numpy_integer_faces(dtype):
    got = spin._faces(np.array(FACES, dtype=dtype), 'faces')
    assert got == FACES
    assert all(type(x) is int for row in got for x in row)

def test_preserves_builtin_integer_input():
    assert spin._faces(FACES, 'faces') == FACES

@pytest.mark.parametrize('bad', [True, np.bool_(False), 1.0, np.float64(1), -1,
                                 2**31, np.uint64(2**63)])
def test_rejects_noninteger_or_nonrepresentable_indices(bad):
    rows = [list(row) for row in FACES]
    rows[0][0] = bad
    with pytest.raises(ValueError):
        spin._faces(rows, 'faces')

def test_keeps_triangle_count_guard():
    with pytest.raises(ValueError, match='four faces'):
        spin._faces(FACES[:3], 'faces')
