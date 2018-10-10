""" Test scripts
"""
import os
import sys
import subprocess
import tempfile

PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(PATH, '../examples/templates')
SYNGAS_PATH = os.path.join(PATH, '../examples/syngas')


def test__automech__syngas_from_rmg():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    sys.stdout.write("{:s}\n".format(tmp_path))

    mech_json = os.path.join(SYNGAS_PATH, 'mechanism.json')
    spc_json = os.path.join(SYNGAS_PATH, 'species.json')
    subprocess.check_call(['automech', 'init_from_rmg',
                           mech_json,
                           spc_json,
                           '-P', tmp_path,
                           '-p'])
    subprocess.check_call(['automech', 'abstractions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'additions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p'])
    subprocess.check_call(['automech', 'migrations', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'additions'),
                           os.path.join(tmp_path, 'additions',
                                        'reactions.csv')])
    subprocess.check_call(['automech', 'additions', 'run',
                           os.path.join(tmp_path, 'additions',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'addition.txt'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'migrations'),
                           os.path.join(tmp_path, 'migrations',
                                        'reactions.csv')])
    subprocess.check_call(['automech', 'migrations', 'run',
                           '-P', os.path.join(tmp_path, 'migrations'),
                           os.path.join(tmp_path, 'migrations',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'migration.txt'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p',
                           'cmd', 'ls'])

    # test divided runs
    subprocess.check_call(['automech', 'divide',
                           'rad-rad', 'rad-rad', 'rad-mol',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'divide',
                           'high-spin', 'high-spin', 'low-spin',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad'),
                           '-p'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol'),
                           os.path.join(tmp_path, 'abstractions', 'rad-mol',
                                        'reactions.csv')])
    subprocess.check_call(['automech', 'abstractions', 'run',
                           os.path.join(tmp_path, 'abstractions', 'rad-mol',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin'),
                           os.path.join(tmp_path, 'abstractions', 'rad-rad',
                                        'high-spin', 'reactions.csv')])
    subprocess.check_call(['automech', 'abstractions', 'run',
                           os.path.join(tmp_path, 'abstractions', 'rad-rad',
                                        'high-spin', 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin'),
                           '-p',
                           'cmd', 'ls'])


def test__automech__syngas_from_chemkin():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    sys.stdout.write("{:s}\n".format(tmp_path))

    subprocess.check_call(['automech', 'init',
                           os.path.join(SYNGAS_PATH, 'mechanism.txt'),
                           os.path.join(SYNGAS_PATH, 'species.csv'),
                           '-P', tmp_path])
    subprocess.check_call(['automech', 'abstractions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'additions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p'])
    subprocess.check_call(['automech', 'migrations', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'additions'),
                           os.path.join(tmp_path, 'additions',
                                        'reactions.csv')])
    subprocess.check_call(['automech', 'additions', 'run',
                           os.path.join(tmp_path, 'additions',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'addition.txt'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'migrations'),
                           os.path.join(tmp_path, 'migrations',
                                        'reactions.csv')])
    subprocess.check_call(['automech', 'migrations', 'run',
                           os.path.join(tmp_path, 'migrations',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'migration.txt'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p',
                           'cmd', 'ls'])

    # test divided runs
    subprocess.check_call(['automech', 'divide',
                           'rad-rad', 'rad-rad', 'rad-mol',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'divide',
                           'high-spin', 'high-spin', 'low-spin',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad'),
                           '-p'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol'),
                           os.path.join(tmp_path, 'abstractions', 'rad-mol',
                                        'reactions.csv')])
    subprocess.check_call(['automech', 'abstractions', 'run',
                           os.path.join(tmp_path, 'abstractions', 'rad-mol',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'csv', 'reindex',
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin'),
                           os.path.join(tmp_path, 'abstractions', 'rad-rad',
                                        'high-spin', 'reactions.csv')])
    subprocess.check_call(['automech', 'abstractions', 'run',
                           os.path.join(tmp_path, 'abstractions', 'rad-rad',
                                        'high-spin', 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin'),
                           '-p',
                           'cmd', 'ls'])


if __name__ == '__main__':
    test__automech__syngas_from_rmg()
    test__automech__syngas_from_chemkin()
