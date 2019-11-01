import sys
from io import StringIO

class CaptureOutErr(object):
    """Context manager to capture the content of stdout and stderr.

    Example:
        >>> with CaptureOutErr() as cm:
        >>> ...run_code()
        >>> print(cm)
    """

    def __enter__(self):
        self.stdout = []
        self.stderr = []
        self._out = StringIO()
        self._err = StringIO()
        sys.stdout = self._out
        sys.stderr = self._err
        return self

    def __exit__(self, *args):
        self.stdout.extend(self._out.getvalue().splitlines())
        self.stderr.extend(self._err.getvalue().splitlines())
        # self.stdout.extend(self._out.getvalue())
        # self.stderr.extend(self._err.getvalue())
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__