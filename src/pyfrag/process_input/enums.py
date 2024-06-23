from enum import Enum


class CalcType(Enum):
    RESTRICTED = 1
    UNRESTRICTED_SPINPOL = 2
    UNRESTRICTED_NO_SPINPOL = 3


class CalcProgram(Enum):
    AMS = 1
    ORCA = 2
    GAUSSIAN = 3
