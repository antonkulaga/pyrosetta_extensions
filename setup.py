from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.2'
DESCRIPTION = 'Pyroseta extensions'
LONG_DESCRIPTION = 'A package with helper functions for pyrosetta'

# Setting up
setup(
    name="pyrosetta_extensions",
    version=VERSION,
    author="antonkulaga (Anton Kualga)",
    author_email="<antonkulaga@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=['pycomfort'],
    keywords=['python', 'utils', 'files', "pyrosetta_extensions"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
