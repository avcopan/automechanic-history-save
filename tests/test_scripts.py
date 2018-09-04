""" Test scripts
"""
import os
import subprocess
import tempfile

PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(PATH, '../examples/templates')
NATGAS_PATH = os.path.join(PATH, '../examples/natgas')
SYNGAS_PATH = os.path.join(PATH, '../examples/syngas')


def test__automech__natgas_from_rmg():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    rmg_mech_json = os.path.join(NATGAS_PATH, 'mechanism.json')
    subprocess.check_call(['automech', 'init',
                           '-j', rmg_mech_json,
                           '-P', tmp_path])


def test__automech__syngas_from_chemkin():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    subprocess.check_call(['automech', 'init',
                           '-m', os.path.join(SYNGAS_PATH, 'mechanism.txt'),
                           '-s', os.path.join(SYNGAS_PATH, 'species.csv'),
                           '-P', tmp_path])
    subprocess.check_call(['automech', 'abstractions', 'init',
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions')])
    subprocess.check_call(['automech', 'additions', 'init',
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'additions')])
    subprocess.check_call(['automech', 'migrations', 'init',
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'migrations')])
    subprocess.check_call(['automech', 'abstractions', 'run',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'abstractions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-y', 'nodes:d',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'additions', 'run',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'addition.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'additions',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'additions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-y', 'nodes:d',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'migrations', 'run',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'migration.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'migrations',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'migrations',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-y', 'nodes:d',
                           'cmd', 'ls'])


if __name__ == '__main__':
    test__automech__syngas_from_chemkin()
