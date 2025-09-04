class PhysError(Exception):
    """Generates an error if a physical constant is missing.
    """
    def __init__(self, name, msg=None):
        if msg is None:
            msg = f'Missing physical quantity: {name}'
        super(PhysError, self).__init__(msg)


class GaussianVersionError(Exception):
    """Generates an error for an unsupported Gaussian version.
    """
    def __init__(self, name, msg=None):
        if msg is None:
            msg = f'Unsupported version of Gaussian: {name}'
        super(GaussianVersionError, self).__init__(msg)


class GaussianFileError(Exception):
    """Generates an error for an unrecognized Gaussian file.
    """
    def __init__(self, fname, msg=None):
        if msg is None:
            msg = f'Unrecognized Gaussian file: {fname}'
        super(GaussianFileError, self).__init__(msg)


class KeywordError(Exception):
    """Generates an error if keyword not found.
    """
    def __init__(self, name, msg=None):
        if msg is None:
            msg = f'Unrecognized keyword: {name}'
        super(KeywordError, self).__init__(msg)
