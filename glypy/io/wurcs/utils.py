import string
from glypy.io.file_utils import ParserError


def base52(x):
    code = []
    if x == 0:
        return string.ascii_letters[0]
    while x > 0:
        i = x % 52
        code.append(i)
        x //= 52
    code = code[::-1]
    n = len(code)
    if n == 1:
        return ''.join([string.ascii_letters[c] for j, c in enumerate(code)])
    else:
        return ''.join([string.ascii_letters[c - 1 if j != n - 1 else c] for j, c in enumerate(code)])


class WURCSError(ParserError):
    pass


class WURCSFeatureNotSupported(WURCSError):
    pass
