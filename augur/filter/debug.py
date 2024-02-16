from functools import wraps
import time
from typing import Callable

from . import constants


def print_debug(message):
    if constants.RUNTIME_DEBUG:
        print(f"DEBUG: {message}")


def add_debugging(func: Callable):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if constants.RUNTIME_DEBUG:
            start_time = time.perf_counter()
            print_debug(f'Starting {func.__name__}.')
            result = func(*args, **kwargs)
            end_time = time.perf_counter()
            total_time = end_time - start_time
            print_debug(f'Function {func.__name__} finished in {total_time:.4f} seconds.')
            return result
        else:
            return func(*args, **kwargs)
    return wrapper
