""" Test scripts
"""
import os
import subprocess
import tempfile
import shutil

PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
BUTANE_PATH = os.path.join(PATH, '../examples/butane')


def test__automech__butane_from_rmg():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    rmg_mech_json = os.path.join(BUTANE_PATH, 'mechanism.json')
    subprocess.check_call(['automech', 'init',
                           '-j', rmg_mech_json,
                           '-P', tmp_path,
                           '-S', 'species.csv',
                           '-R', 'reactions.csv',
                           '-G', 'geoms'])


def test__automech__butane_from_chemkin():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    chemkin_mech_txt = os.path.join(BUTANE_PATH, 'mechanism.txt')
    spc_csv = os.path.join(BUTANE_PATH, 'species.csv')
    subprocess.check_call(['automech', 'init',
                           '-m', chemkin_mech_txt,
                           '-s', spc_csv,
                           '-P', tmp_path,
                           '-S', 'species.csv',
                           '-R', 'reactions.csv',
                           '-G', 'geoms'])


def test__automech():
    """ test the automech_init script
    """
    mech_file_path = os.path.join(DATA_PATH, 'mechanism.txt')
    spc_csv_path = os.path.join(DATA_PATH, 'species.csv')
    tmp_path = tempfile.mkdtemp()
    spc_csv_out_path = os.path.join(tmp_path, 'species.csv')
    geom_dir_path = os.path.join(tmp_path, 'geom')
    rxn_csv_out_path = os.path.join(tmp_path, 'reactions.csv')
    subprocess.check_call(['automech_init',
                           '-m', mech_file_path,
                           '-s', spc_csv_path,
                           '-S', spc_csv_out_path,
                           '-G', geom_dir_path,
                           '-x', 'OHV', 'CHV',
                           '-R', rxn_csv_out_path])

    """ test the automech_abstr_init script
    """
    rxn_csv_path = os.path.join(DATA_PATH, 'reactions-sample.csv')
    spc_csv_path = spc_csv_out_path
    abstr_csv_out_path = os.path.join(tmp_path, 'abstractions.csv')
    abstr_dir_path = os.path.join(tmp_path, 'abstr')
    subprocess.check_call(['automech_abstr_init',
                           '-r', rxn_csv_path,
                           '-s', spc_csv_path,
                           '-R', rxn_csv_out_path,
                           '-A', abstr_csv_out_path,
                           '-D', abstr_dir_path])

    """ test the automech_abstr_run script
    """
    abstr_csv_path = abstr_csv_out_path
    tmp_file_path = os.path.join(DATA_PATH, 'template.txt')
    subprocess.check_call(['automech_abstr_run',
                           '-a', abstr_csv_path,
                           '-t', tmp_file_path,
                           '-n', 'b444', 'b445',
                           '-x', '0-3', '5', '7-9', '17',
                           'cmd', 'ls', '-la'])
    shutil.rmtree(tmp_path)


if __name__ == '__main__':
    test__automech__butane_from_rmg()
    test__automech__butane_from_chemkin()
