from enum import Enum, auto


class JobStatus(Enum):
    """Enum for job status"""

    INITIALIZED = auto()
    PREPARED = auto()
    RUNNING = auto()
    FINISHED = auto()
    FAILED = auto()
    CANCELLED = auto()
    UNKNOWN = auto()

    def __str__(self):
        return self.name.lower()


class CalcType(Enum):
    RESTRICTED = 1
    UNRESTRICTED_SPINPOL = 2
    UNRESTRICTED_NO_SPINPOL = 3


class CalcProgram(Enum):
    AMS = 1
    ORCA = 2
    GAUSSIAN = 3
