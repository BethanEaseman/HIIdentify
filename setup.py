'''Setup for the module'''

__author__ = 'Bethan Easeman'
__version__ = '0.0.2'

#import sys
from os import path
from setuptools import setup, find_packages

def install():
    '''The installer'''

    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as file_id:
        long_description = file_id.read()
    short_description = ('Identifies HII regions from H alpha emission line maps.')

    setup(name='HIIdentify',
          version=__version__,
          license="GPLv3",
          description=short_description,
          long_description=long_description,
          long_description_content_type='text/markdown',
          author=__author__,
          author_email='be329@bath.ac.uk',
          packages=find_packages(),
          # keywords=[],
          # zip_safe=False,
          url='https://github.com/BethanEaseman/HIIdentify',
          # project_urls={},
          # classifiers=[],
          install_requires=['numpy', 'astropy', 'matplotlib'],
          python_requires=">=3.6",
          # entry_points={}
    )

if __name__ == "__main__":
    install()
