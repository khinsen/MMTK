# Thread manager for long-running threads (MD etc.)
#
# Written by Konrad Hinsen
# last revision: 2010-7-22
#

_undocumented = 1

try:
    import threading
    if not hasattr(threading, 'Thread'):
        threading = None
except ImportError:
    threading = None

_threads = []

def registerThread(thread):
    _threads.append(thread)
    _cleanup()

def activeThreads():
    _cleanup()
    return _threads

def waitForThreads():
    while _threads:
        _threads[0].join()
        _cleanup()

def _cleanup():
    i = 0
    while i < len(_threads):
        if _threads[i].isAlive():
            i = i + 1
        else:
            del _threads[i]

if threading:

    import MMTK_state_accessor

    class TrajectoryGeneratorThread(threading.Thread):

        def __init__(self, universe, target, args, state_accessor):
            threading.Thread.__init__(self, group = None,
                                      name = 'MMTK thread',
                                      target = target,
                                      args = args)
            self.universe = universe
            self.function = (target, args)
            self.state_accessor = state_accessor
            self.start()
            registerThread(self)

        def run(self):
            target, args = self.function
            self.universe.acquireConfigurationChangeLock()
            target(*args)
            self.universe.releaseConfigurationChangeLock()

        def copyState(self):
            if self.state_accessor is None:
                return None
            else:
                return self.state_accessor.copyState()

else:

    # a fake thread class that raises an exception
    class TrajectoryGeneratorThread(object):

        def __init__(self, universe, target, args, kwargs):
            raise OSError("background processing not available")
