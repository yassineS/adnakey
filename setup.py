from distutils.core import setup
from setuptools import find_packages

README = open('README.rst').read()

setup(name='ReichKey',
      version='0.1',
      description = "Next-generation Sequencing Calling Pipeline, tailored to Cteam at Reich Lab.",
      license='Non-commercial',
      long_description=README,
      packages=find_packages(),
      scripts=['bin/cteamkey'],
      install_requires=['ipdb']
)
