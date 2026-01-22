"""
This file, augur/io/json.py, started as a copy of lib/id3c/json.py and the
functions shorten_left(), contextualize_char(), and mark_char() from
lib/id3c/utils.py in the https://github.com/seattleflu/id3c repo as of commit
911e7d7bfccc4d050e63e6f73d7a7e59fa1a80e8, licensed under the MIT License.

Subsequent modifications (as recorded in Augur's own version control history)
are licensed under the same terms as the rest of Augur, the GNU AGPL 3.0.

The LICENSE file included in ID3C's repo is copied below verbatim::

    MIT License

    Copyright (c) 2018 Brotman Baty Institute

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
"""
import json
import numpy as np
import os
import pandas as pd
from collections import OrderedDict
from datetime import date, datetime, time, timedelta
from io import RawIOBase
from isodate import duration_isoformat
from typing import Iterable
from uuid import UUID
from augur.io.file import open_file


def as_json(value):
    """
    Converts *value* to a JSON string using our custom :class:`JsonEncoder`.

    The custom encoder supports serialization of :class:`~datetime.date` objects:

    >>> as_json(date(year=2024, month=7, day=17))
    '"2024-07-17"'

    :class:`~datetime.datetime` objects:

    >>> as_json(datetime(year=2024, month=7, day=17, hour=11, minute=38))
    '"2024-07-17T11:38:00"'

    :class:`~datetime.time` objects:

    >>> as_json(time(hour=11, minute=38))
    '"11:38:00"'

    :class:`~datetime.timedelta` objects:

    >>> as_json(timedelta(days=42))
    '"P42D"'

    and :class:`~uuid.UUID` objects:

    >>> as_json(UUID(int=147952133113722764103424939352979237618))
    '"6f4e8b5a-8500-4928-b7ae-dc098a256af2"'
    """
    return json.dumps(value, allow_nan = False, cls = JsonEncoder)


def load_json(value):
    """
    Converts *value* from a JSON string with better error messages.
    Raises an :exc:`augur.io.json.JSONDecodeError` which provides improved error
    messaging, compared to :exc:`json.JSONDecodeError`, when stringified.
    """
    try:
        return json.loads(value)
    except json.JSONDecodeError as e:
        raise JSONDecodeError(e) from e


MINIFY_THRESHOLD_MB = 5


def write_json(data, file, minify=None, minify_threshold_mb=MINIFY_THRESHOLD_MB, indent=2):
    """
    Write ``data`` as JSON to the given ``file``, creating parent directories
    if necessary.

    Parameters
    ----------
    data : dict
        data to write out to JSON
    file
        file path or handle to write to
    minify : bool or None, optional
        Control output minification. ``True`` forces minified output, ``False``
        forces non-minified output, ``None`` (default) uses auto-detection based
        on *minify_threshold_mb*.
        A truthy value in the environment variable :envvar:`AUGUR_MINIFY_JSON`
        also forces minified output.
    minify_threshold_mb : int or float, optional
        Threshold in megabytes above which output is automatically minified.
        Only applies when *minify* is None.
    indent : int or None, optional
        JSON indentation level when not minifying.

    Raises
    ------
    OSError
    """
    if isinstance(file, (str, os.PathLike)):
        #in case parent folder does not exist yet
        parent_directory = os.path.dirname(file)
        if parent_directory and not os.path.exists(parent_directory):
            try:
                os.makedirs(parent_directory)
            except OSError: #Guard against race condition
                if not os.path.isdir(parent_directory):
                    raise

    # Should output be minified?
    # Order of precedence:
    # 1. 'minify' parameter
    # 2. 'AUGUR_MINIFY_JSON' environment variable
    # 3. Automatically determine based on size of data
    if minify is True:
        effective_indent = None
    elif minify is False:
        effective_indent = indent
    elif os.environ.get("AUGUR_MINIFY_JSON"):
        effective_indent = None
    elif json_size(data) > minify_threshold_mb * 10**6:
        effective_indent = None
    else:
        effective_indent = indent

    with open_file(file, 'w', encoding='utf-8') as handle:
        sort_keys = False if isinstance(data, OrderedDict) else True
        json.dump(data, handle, indent=effective_indent, sort_keys=sort_keys, cls=AugurJSONEncoder)


class AugurJSONEncoder(json.JSONEncoder):
    """
    A custom JSONEncoder subclass to serialize data types used for various data
    stored in dictionary format.
    """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pd.Series):
            return obj.tolist()
        return super().default(obj)


class BytesWrittenCounterIO(RawIOBase):
    """Binary stream to count the number of bytes sent via write()."""
    def __init__(self):
        self.written = 0
        """Number of bytes written."""

    def write(self, b):
        n = len(b)
        self.written += n
        return n


def json_size(data):
    """Return size in bytes of a Python object in JSON string form."""
    with BytesWrittenCounterIO() as counter:
        write_json(data, counter, minify=False)
    return counter.written


def dump_ndjson(iterable: Iterable) -> None:
    """
    :func:`print` *iterable* as a set of newline-delimited JSON records.
    """
    for item in iterable:
        print(as_json(item))


def load_ndjson(file: Iterable[str], ignore_empty_lines = True) -> Iterable:
    """
    Load newline-delimited JSON records from *file*. Ignore empty lines
    in the file by default.
    """
    for line in file:
        # Skip empty lines when requested
        if ignore_empty_lines and not line.strip():
            continue

        yield load_json(line)


class JsonEncoder(json.JSONEncoder):
    """
    Encodes Python values into JSON for non-standard objects.
    """

    def __init__(self, *args, **kwargs):
        """
        Disallows the floating point values NaN, Infinity, and -Infinity.
        Python's :class:`json` allows them by default because they work with
        JSON-as-JavaScript, but they don't work with spec-compliant JSON
        parsers.
        """
        kwargs["allow_nan"] = False
        super().__init__(*args, **kwargs)

    def default(self, value):
        """
        Returns *value* as JSON or raises a TypeError.
        Serializes:
        * :class:`~datetime.date` using :meth:`~datetime.date.isoformat()`
        * :class:`~datetime.datetime` using :meth:`~datetime.datetime.isoformat()`
        * :class:`~datetime.time` using :meth:`~datetime.time.isoformat()`
        * :class:`~datetime.timedelta` using ``isodate.duration_isoformat()``
        * :class:`~uuid.UUID` using ``str()``
        """
        if isinstance(value, (date, datetime, time)):
            return value.isoformat()

        elif isinstance(value, timedelta):
            return duration_isoformat(value)

        elif isinstance(value, UUID):
            return str(value)

        else:
            # Let the base class raise the TypeError
            return super().default(value)


class JSONDecodeError(json.JSONDecodeError):
    """
    Subclass of :class:`json.JSONDecodeError` which contextualizes the
    stringified error message by including a snippet of the JSON source input.
    Typically you won't need to ever reference this class directly.  It will be
    raised by :func:`load_json` and be caught by except blocks which catch the
    standard :class:`json.JSONDecodeError`.

    Examples
    --------
    >>> load_json('{foo: "bar"}')
    Traceback (most recent call last):
        ...
    augur.io.json.JSONDecodeError: Expecting property name enclosed in double quotes: line 1 column 2 (char 1): '{▸▸▸f◂◂◂oo: "bar"}'
    >>> load_json('not json')
    Traceback (most recent call last):
        ...
    augur.io.json.JSONDecodeError: Expecting value: line 1 column 1 (char 0): 'not json'
    >>> load_json("[0, 1, 2, 3, 4, 5")
    Traceback (most recent call last):
        ...
    augur.io.json.JSONDecodeError: Expecting ',' delimiter: line 1 column 18 (char 17): unexpected end of document: '…, 3, 4, 5'
    >>> load_json("[\\n")
    Traceback (most recent call last):
        ...
    augur.io.json.JSONDecodeError: Expecting value: line 2 column 1 (char 2): unexpected end of document: '[\\n'
    >>> load_json("\\n")
    Traceback (most recent call last):
        ...
    augur.io.json.JSONDecodeError: Expecting value: line 2 column 1 (char 1): unexpected end of document: '\\n'
    >>> load_json('')
    Traceback (most recent call last):
        ...
    augur.io.json.JSONDecodeError: Expecting value: line 1 column 1 (char 0): (empty source document)
    """
    CONTEXT_LENGTH = 10

    def __init__(self, exc: json.JSONDecodeError):
        super().__init__(exc.msg, exc.doc, exc.pos)

    def __str__(self):
        error = super().__str__()

        if self.doc:
            if self.pos == 0 and self.msg == "Expecting value":
                # Most likely not a JSON document at all, so show the whole thing.
                context = repr(self.doc)
            elif self.pos > 0 and self.pos == len(self.doc):
                context = "unexpected end of document: " + repr(shorten_left(self.doc, self.CONTEXT_LENGTH, "…"))
            else:
                context = repr(contextualize_char(self.doc, self.pos, self.CONTEXT_LENGTH))
        else:
            context = "(empty source document)"

        return f"{error}: {context}"


def shorten_as_json(value, length, placeholder) -> str:
    """
    Converts *value* to JSON with :py:func:`as_json` and then truncates it to a
    maximum *length* (if necessary), indicating truncation with the given
    *placeholder*.

    >>> shorten_as_json({'hello': 'world', 'x': 42}, 100, '…')
    '{"hello": "world", "x": 42}'
    >>> shorten_as_json({'hello': 'world', 'x': 42}, 21, '…')
    '{"hello": "world", …}'

    For readability, the outermost JSON value delimiter at the right hand side,
    i.e. ``}``, ``]``, or ``"``, is preserved, such that *placeholder* will put
    placed "inside" *value*'s JSON representation.

    >>> shorten_as_json([1,2,3,'four',5,6], 20, '…')
    '[1, 2, 3, "four", …]'
    >>> shorten_as_json([1,2,3,'four',5,6], 15, '…')
    '[1, 2, 3, "fo…]'

    The maximum *length* must be at least two characters longer than the length
    of the *placeholder*.

    >>> shorten_as_json({'foo': 'bar'}, 4, '...')
    Traceback (most recent call last):
        ...
    ValueError: maximum length (4) must be two greater than length of placeholder (3), i.e. at least 5
    >>> shorten_as_json({'foo': 'bar'}, 5, '...')
    '{...}'
    """
    min_length = len(placeholder) + 2
    if length < min_length:
        raise ValueError(f"maximum length ({length}) must be two greater than length of placeholder ({len(placeholder)}), i.e. at least {min_length}")

    json_value = as_json(value)

    if len(json_value) > length:
        return json_value[0:length - len(placeholder) - 1] + placeholder + json_value[-1:]
    else:
        return json_value


def shorten_left(text, length, placeholder):
    """
    Truncate the left end of *text* to a maximum *length* (if necessary),
    indicating truncation with the given *placeholder*.

    The maximum *length* must be longer than the length of the *placeholder*.

    Behaviour is slightly different than :py:func:`textwrap.shorten` which is
    intended for shortening sentences and works at the word, not character,
    level.

    Examples
    --------
    >>> shorten_left("foobar", 6, "...")
    'foobar'
    >>> shorten_left("foobarbaz", 6, "...")
    '...baz'
    >>> shorten_left("foobar", 3, "...")
    Traceback (most recent call last):
        ...
    ValueError: maximum length (3) must be greater than length of placeholder (3)
    """
    if length <= len(placeholder):
        raise ValueError(f"maximum length ({length}) must be greater than length of placeholder ({len(placeholder)})")

    if len(text) > length:
        return placeholder + text[-(length - len(placeholder)):]
    else:
        return text


def contextualize_char(text, idx, context = 10):
    """
    Marks the *idx* char in *text* and snips out a surrounding amount of
    *context*.

    Avoids making a copy of *text* before snipping, in case *text* is very
    large.

    Examples
    --------
    >>> contextualize_char('hello world', 0, context = 4)
    '▸▸▸h◂◂◂ello…'
    >>> contextualize_char('hello world', 5, context = 3)
    '…llo▸▸▸ ◂◂◂wor…'
    >>> contextualize_char('hello world', 5, context = 100)
    'hello▸▸▸ ◂◂◂world'
    >>> contextualize_char('hello world', 10)
    'hello worl▸▸▸d◂◂◂'
    >>> contextualize_char('hello world', 2, context = 0)
    '…▸▸▸l◂◂◂…'

    >>> contextualize_char('hello world', 11)
    Traceback (most recent call last):
        ...
    IndexError: string index out of range
    """
    if context < 0:
        raise ValueError("context must be positive")

    start = max(0, idx - context)
    end   = min(len(text), idx + context + 1)
    idx   = min(idx, context)

    start_placeholder = "…" if start > 0         else ""
    end_placeholder   = "…" if end   < len(text) else ""

    return start_placeholder + mark_char(text[start:end], idx) + end_placeholder


def mark_char(text, idx):
    """
    Prominently marks the *idx* char in *text*.

    Examples
    --------
    >>> mark_char('hello world', 0)
    '▸▸▸h◂◂◂ello world'
    >>> mark_char('hello world', 2)
    'he▸▸▸l◂◂◂lo world'
    >>> mark_char('hello world', 10)
    'hello worl▸▸▸d◂◂◂'

    >>> mark_char('hello world', 11)
    Traceback (most recent call last):
        ...
    IndexError: string index out of range

    >>> mark_char('', 0)
    Traceback (most recent call last):
        ...
    IndexError: string index out of range
    """
    return text[0:idx] + '▸▸▸' + text[idx] + '◂◂◂' + text[idx+1:]
