""" Test the automech CLI
"""
import os
import sys
import tempfile
import subprocess

AUTOMECH_CMD = 'automech'
PATH = os.path.dirname(os.path.realpath(__file__))
HEPTANE_PATH = os.path.join(PATH, '../examples/heptane')


def test__help():
    """ test `automech -h`
    """
    subprocess.check_call([AUTOMECH_CMD, '-h'])


def test__chemkin__help():
    """ test `automech chemkin -h`
    """
    subprocess.check_call([AUTOMECH_CMD, 'chemkin', '-h'])


def test__chemkin__to_csv():
    """ test `automech chemkin to_csv`
    """
    subprocess.check_call([AUTOMECH_CMD, 'chemkin', 'to_csv', '-h'])

    calling_dir = os.getcwd()
    routine_dir = tempfile.mkdtemp()

    mech_txt = os.path.join(HEPTANE_PATH, 'mechanism.txt')
    ther_txt = os.path.join(HEPTANE_PATH, 'thermo_data.txt')

    sys.stdout.write("{:s}\n".format(routine_dir))

    os.chdir(routine_dir)
    subprocess.check_call([AUTOMECH_CMD, 'chemkin', 'to_csv',
                           mech_txt, ther_txt, '-p'])
    os.chdir(calling_dir)


def test__species__help():
    """ test `automech species -h`
    """
    subprocess.check_call([AUTOMECH_CMD, 'species', '-h'])


def test__species__expand_stereo():
    """ test `automech species expand_stereo`
    """
    subprocess.check_call([AUTOMECH_CMD, 'species', 'expand_stereo', '-h'])

    calling_dir = os.getcwd()
    routine_dir = tempfile.mkdtemp()

    geom_type_key = 'smi_'
    spc_csv = os.path.join(HEPTANE_PATH, 'smiles.csv')

    sys.stdout.write("{:s}\n".format(routine_dir))

    os.chdir(routine_dir)
    subprocess.check_call([AUTOMECH_CMD, 'species', 'expand_stereo',
                           geom_type_key, spc_csv, '-p'])
    os.chdir(calling_dir)


if __name__ == '__main__':
    # test__help()
    # test__chemkin__help()
    # test__chemkin__to_csv()
    # test__species__help()
    test__species__expand_stereo()
