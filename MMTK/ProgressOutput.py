# This module implements a TrajectoryAction that outputs the current progress.
#
# Written by Dmitri Iouchtchenko
#

"""
Progress bars for simulations

If the `progressbar <https://pypi.python.org/pypi/progressbar>`_ package
(version 2.3 or greater) is installed, :class:`ProgressOutput` outputs a
progress bar during a simulation::

    37% |###############                           | 0:00:17 / ETA:  0:00:29

Otherwise it only prints out the percentage of steps done so far.
"""

__docformat__ = 'restructuredtext'

from distutils.version import LooseVersion
from warnings import warn

pb_avail = True

try:
    import progressbar as pb

    if LooseVersion(pb.__version__) < LooseVersion('2.3'):
        pb_avail = False
        warn('Package "progressbar" must be at least version 2.3 for fancy progress bars.')
except ImportError:
    # Resort to the plain output.
    pb_avail = False


from sys import stdout

from MMTK.Trajectory import LogOutput


class _ProgressIO(object):
    """
    A sufficiently file-like object that can be passed to
    :class:`MMTK.Trajectory.LogOutput` to get the current step number of the
    trajectory.
    """

    def __init__(self, steps, skip):
        self._steps = steps
        self._skip = skip

        self._started = False
        self._finished = False

    def write(self, line, *args, **kwargs):
        """
        Pretend that we're a file object that can be written to.
        """

        if self._finished:
            return

        # Here we rely on the fact that PyTrajectory_Output always prints the
        # step number in this particular format.
        if line.startswith('Step '):
            step = int(line[5:])

            if not self._started:
                self._start()
                self._started = True

            if self._steps - step <= self._skip:
                # This is the last time we get to update, so fast forward to
                # the end.
                self._update(self._steps)
                self._finish()
                self._finished = True
            else:
                self._update(step)

    def _start(self):
        pass

    def _update(self, step):
        pass

    def _finish(self):
        pass


class _ProgressBarIO(_ProgressIO):
    """
    Use progressbar to display a fancy progress bar.
    """

    def __init__(self, *args, **kwargs):
        _ProgressIO.__init__(self, *args, **kwargs)

        widgets = [pb.Percentage(), ' ', pb.Bar(), ' ', pb.Timer('%s'), ' / ', pb.ETA()]
        self._pbar = pb.ProgressBar(widgets=widgets, maxval=self._steps)

    def _start(self):
        self._pbar.start()

    def _update(self, step):
        self._pbar.update(step)

    def _finish(self):
        self._pbar.finish()


class _PercentIO(_ProgressIO):
    """
    Display the current progress as a percentage.
    """

    def __init__(self, *args, **kwargs):
        _ProgressIO.__init__(self, *args, **kwargs)

    def _update(self, step):
        print '\r{0}%'.format(100 * step // self._steps),
        stdout.flush()

    def _finish(self):
        print


class ProgressOutput(LogOutput):
    """
    A :class:`~MMTK.Trajectory.TrajectoryAction` that prints out the trajectory
    progress.
    """

    def __init__(self, steps, skip=None):
        """
        :param steps: number of steps in the trajectory
        :type name: int
        :param skip: number of steps to skip between updates
        :type name: int
        """

        if skip is None:
            # By default request about 10 updates per percent.
            skip = steps // 1000

        if skip < 1:
            skip = 1

        if pb_avail:
            f = _ProgressBarIO(steps, skip)
        else:
            f = _PercentIO(steps, skip)

        LogOutput.__init__(self, f, data=(), skip=skip)
