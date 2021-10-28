from ...core.errors import PlamsError

__all__ = ["AMSPipeError", "AMSPipeDecodeError", "AMSPipeLogicError", "AMSPipeRuntimeError", "AMSPipeUnknownVersionError", "AMSPipeUnknownMethodError", "AMSPipeUnknownArgumentError", "AMSPipeInvalidArgumentError"]

class AMSPipeError(PlamsError):
    """Base class for exceptions mapped to errors defined by the pipe protocol.

    Not to be used directly.
    """
    def __init__(self, message, argument=None):
        self.status = _class2code[self.__class__]
        self.message = message
        self.method = None
        self.argument = argument

    def __str__(self):
        return ": ".join(s for s in (self.method, self.argument, self.message) if s is not None)

    @staticmethod
    def from_message(msg):
        exc = _code2class[msg["status"]](msg.get("message"), msg.get("argument"))
        exc.method = msg.get("method")

        return exc

    def to_message(self):
        msg = {"status": self.status}
        if self.message is not None:
            msg["message"] = self.message
        if self.method is not None:
            msg["method"] = self.method
        if self.argument is not None:
            msg["argument"] = self.argument

        return msg



class AMSPipeDecodeError(AMSPipeError):
    pass

class AMSPipeLogicError(AMSPipeError):
    pass

class AMSPipeRuntimeError(AMSPipeError):
    pass

class AMSPipeUnknownVersionError(AMSPipeError):
    pass

class AMSPipeUnknownMethodError(AMSPipeError):
    pass

class AMSPipeUnknownArgumentError(AMSPipeError):
    pass

class AMSPipeInvalidArgumentError(AMSPipeError):
    pass


_code2class = {
    1: AMSPipeDecodeError,
    2: AMSPipeLogicError,
    3: AMSPipeRuntimeError,
    4: AMSPipeUnknownVersionError,
    5: AMSPipeUnknownMethodError,
    6: AMSPipeUnknownArgumentError,
    7: AMSPipeInvalidArgumentError
}
_class2code = {c: i for i, c in _code2class.items()}
