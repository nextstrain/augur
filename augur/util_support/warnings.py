import warnings

# Set up warnings & exceptions
_deprecationWarningsEmitted = False

def deprecated(message):
    warnings.warn(message, DeprecationWarning, stacklevel=2)
    global _deprecationWarningsEmitted
    _deprecationWarningsEmitted=True

def warn(message):
    warnings.warn(message, UserWarning, stacklevel=2)

def configure_warnings():
    # This function must be run before any calls to `warn()` and `deprecated()`, however
    # don't execute the code here since that would apply this configuration to _all_
    # augur commands due to the way all commands are pulled in by the augur runner (augur/__init__.py)
    #
    # Intended usage: call this function at the top of the entry function for the augur (sub)command, e.g.
    #       def run(args):
    #           configure_warnings()
    #           ...
    def customformatwarning(message, category, filename, lineno, line=None):
        if category.__name__ == "UserWarning":
            return "WARNING: {}\n\n".format(message)
        if category.__name__ == "DeprecationWarning":
            return "DEPRECATED: {}\n\n".format(message)
        return "{}\n".format(message)

    warnings.formatwarning = customformatwarning
    warnings.simplefilter("default") # show DeprecationWarnings by default

def deprecationWarningsEmitted():
    return _deprecationWarningsEmitted