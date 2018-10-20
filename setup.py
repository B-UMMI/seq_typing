import pathotyping

from setuptools import setup

VERSION = pathotyping.__version__

with open('README.md') as fh:
    README = fh.read()

setup(
    name='seq_typing',
    version='{}'.format(VERSION),
    packages=['seqtyping',
              'seqtyping.modules'],
    package_dir={'seqtyping': 'seqtyping'},
    package_data={'pathotyping': ['../.git/*', '../.git/*/*', '../.git/*/*/*',
                                  'reference_sequences/*/*']},
    include_package_data=True,
    data_files=[('', ['LICENSE'])],
    install_requires=[
        'ReMatCh', 'biopython'
    ],
    description='Determines which reference sequence is more likely to be present in a given sample',
    long_description=README,
    long_description_content_type='text/markdown',
    keywords=['reference mapping', 'sequence Blast search', 'typing'],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Environment :: Console',
        'Operating System :: Unix',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    url='https://github.com/B-UMMI/seq_typing',
    author='Miguel P. Machado',
    author_email='mpmachado@medicina.ulisboa.pt',
    license='GPL3',
    # To use entry_points with .py the first folder cannot have the same name of the script
    entry_points={
        'console_scripts': [
            'seq_typing.py = seqtyping.seq_typing:main',
            'ecoli_stx_subtyping.py = seqtyping.ecoli_stx_subtyping:main',
            'get_stx_db.py = seqtyping.modules.get_stx_db:main'
        ]
    },
    python_requires='>=3.4'
)
