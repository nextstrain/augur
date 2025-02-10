import random

def get_random_unicode(length):
    """Adapted from <https://stackoverflow.com/a/21666621>"""

    try:
        get_char = unichr
    except NameError:
        get_char = chr

    # code point ranges to be sampled
    include_ranges = [
        # ASCII non-control characters excluding single quote (0x27) and backslash (0x5c)
        (0x20, 0x26), (0x28, 0x5b), (0x5d, 0x7e)
    ]

    alphabet = [
        get_char(code_point) for current_range in include_ranges
            for code_point in range(current_range[0], current_range[1] + 1)
    ]
    return ''.join(random.choice(alphabet) for _ in range(length))


if __name__ == "__main__":
    for i in range(10):
        print('>' + get_random_unicode(5))
        print("ATGC")

