import pathlib
import os
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# specify requirements of your package here
REQUIREMENTS = ['biopython', 'numpy', 'pandas']

setup(name='stacksPairwise',
      version='0.0.0',
      description='Calculate pairwise divergence (pairwise pi) from Stacks `samples.fa` output fle',
      long_description=README,
      long_description_content_type="text/markdown",
      url='https://github.com/gibsonmatt/stacks-pairwise',
      author='Matt Gibson',
      author_email='matthewjsgibson@gmail.com',
      license='MIT',
      packages=['stacksPairwise'],
      install_requires=REQUIREMENTS,
      entry_points={
        "console_scripts": [
            "stacksPairwise=stacksPairwise.__main__:main"
        ]
    },
      keywords='genetics genotyping sequencing Stacks'
)
