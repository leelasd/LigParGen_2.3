from setuptools import setup, find_packages

setup(
  name = 'LigParGen',
  packages = ['LigParGen'], # this must be the same as the name above
  version = '3.0',
  description = 'Python script to provide BOSS generated OPLS-AA/CM1A(-LBCC) parameters for organic molecules and ligands.',
  author = 'Leela S. Dodda, Matthew C. Robinson',
  author_email = 'leela.dodda@yale.edu,matthew.robinson@yale.edu',
  license='MIT',
  url='https://bitbucket.org/leelasd/ligpargen_2017_sep18',
  keywords = ['computational chemistry', 'force fields', 'molecular dynamics'],
  classifiers = [
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering',
        ],
  python_requires='>=3.11,<3.13',
  install_requires=[
        'numpy==2.4.6',
        'pandas==2.3.3',
        'rdkit==2025.9.6',
        'networkx==3.5',
        'openbabel==3.2.1',
        ],
  entry_points={
        'console_scripts': [
            'LigParGen=LigParGen.Converter:main',
        ],
    },

)
