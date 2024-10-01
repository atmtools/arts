from setuptools import Distribution, setup


class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name"""

    def has_ext_modules(foo):
        return True


setup(distclass=BinaryDistribution)
