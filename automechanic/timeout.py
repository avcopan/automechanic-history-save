""" defines a timeout decorator
"""
from functools import wraps
import errno
import os
import signal


class TimeoutError(Exception):
    """ TimeoutError exception class
    """
    pass


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    """ decorator to time out after a given interval
    """
    def _decorator(func):

        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def _wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(_wrapper)

    return _decorator
