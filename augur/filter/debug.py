from . import constants


def print_debug(message):
    if constants.RUNTIME_DEBUG:
        print(f"DEBUG: {message}")
