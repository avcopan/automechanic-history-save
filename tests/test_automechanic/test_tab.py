""" test the automechanic.tab module
"""
import tempfile
import numpy
from automechanic import tab

A_KEY = 'a'
B_KEYS = ('b1', 'b2', 'b3')
C_KEY = 'c'
VALS = numpy.arange(5*5).reshape((5, 5))
VALS[:, 4] = numpy.sum(VALS[:, 1:4], axis=1)
KEYS = (A_KEY,) + tuple(B_KEYS) + (C_KEY,)
TYPS = (float, int, int, int, float)
TBL = tab.from_records(VALS, KEYS, typs=TYPS)


def test__fancy_iter():
    """ test tab.fancy_iter
    """
    tbl = TBL[list((A_KEY,) + B_KEYS)]

    tbl = tbl[tbl[A_KEY] > 5]

    fancy_keys = (A_KEY, B_KEYS)
    fancy_vals = tuple(tab.fancy_iter(tbl, fancy_keys))
    assert fancy_vals == ((10.0, (11, 12, 13)), (15.0, (16, 17, 18)),
                          (20.0, (21, 22, 23)))

    fancy_vals = tuple(tab.fancy_enum(tbl, fancy_keys))
    assert fancy_vals == ((2, (10.0, (11, 12, 13))), (3, (15.0, (16, 17, 18))),
                          (4, (20.0, (21, 22, 23))))


def test__from_records():
    """ test tab.from_records
    """
    tbl = TBL[list((A_KEY,) + B_KEYS)]

    fancy_keys = (A_KEY, B_KEYS)
    fancy_vals = tuple(tab.fancy_iter(tbl, fancy_keys))
    group_typs = (float, int)
    tbl_ = tab.from_records(fancy_vals, fancy_keys, group_typs)

    assert tab.equal(tbl, tbl_)


def test__read_csv():
    """ test tab.read_csv
    """
    tbl = TBL[list((A_KEY,) + B_KEYS)]
    tbl2 = tab.from_starmap(tbl[tbl[A_KEY].isin((5, 15, 20))],
                            sum, [B_KEYS], ['c'], [int])

    _, tbl2_fle = tempfile.mkstemp()
    tab.write_csv(tbl2, tbl2_fle)
    tbl2_ = tab.read_csv(tbl2_fle)
    assert tab.equal(tbl2, tbl2_)


def test__update():
    """ test tab.update
    """
    tbl = TBL[list((A_KEY,) + B_KEYS)]
    tbl2 = tab.from_starmap(tbl[tbl[A_KEY].isin((5, 15, 20))],
                            sum, [B_KEYS], ['c'], [int])

    tbl3 = tab.update(tbl, tbl2)

    tbl3_ = tab.from_records(vals=tab.vals_(TBL),
                             keys=tab.keys_(TBL),
                             typs=tab.typs_(TBL),
                             idxs=tab.idxs_(TBL))
    tbl3_.loc[~tbl3_[A_KEY].isin((5, 15, 20)), C_KEY] = tab.NAN
    assert (numpy.array_equal(tab.keys_(tbl3), tab.keys_(tbl3_)) and
            numpy.array_equal(tab.typs_(tbl3), tab.typs_(tbl3_)) and
            numpy.array_equal(tab.idxs_(tbl3), tab.idxs_(tbl3_)) and
            numpy.array_equal(tab.vals_(tbl3.fillna(value=-99)),
                              tab.vals_(tbl3_.fillna(value=-99))))


def test__set_typs():
    """ test tab.change_type
    """
    tbl_ref = tab.from_records(
        VALS, KEYS, typs=(int, float, float, float, float))
    tbl = tab.set_typs(TBL, (A_KEY, B_KEYS), (int, float))
    assert tab.equal(tbl, tbl_ref)


def test__has_keys():
    """ test tab.has_keys
    """
    assert tab.has_keys(TBL, (A_KEY, B_KEYS)) is True
    assert tab.has_keys(TBL, (A_KEY, B_KEYS, 'd')) is False


if __name__ == '__main__':
    test__fancy_iter()
    test__from_records()
    test__read_csv()
    test__update()
    test__set_typs()
    test__has_keys()
