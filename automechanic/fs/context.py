""" filesystem context management
"""
import os


def enter(root_pth):
    """ enter the filesystem, creating its root directory if needed
    """
    if not os.path.isdir(root_pth):
        os.mkdir(root_pth)
    return _EnterFilesystem(root_pth=root_pth)


class _EnterFilesystem():

    def __init__(self, root_pth):
        self.root_pth = root_pth
        self.work_pth = os.getcwd()

    def __enter__(self):
        os.chdir(self.root_pth)

    def __exit__(self, _exc_type, _exc_value, _traceback):
        os.chdir(self.work_pth)
