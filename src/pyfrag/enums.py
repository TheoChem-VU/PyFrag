from enum import Enum


class JobStatus(Enum):
    """Enum for job status"""

    INITIALIZED = 0
    RUNNING = 1
    FINISHED = 2
    FAILED = 3
    CANCELLED = 4
    UNKNOWN = 5

    def __str__(self):
        return self.name.lower()
