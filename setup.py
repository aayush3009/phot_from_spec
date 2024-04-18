from setuptools import setup, find_packages

setup(
    name='phot_from_spec',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'seaborn',
        'astropy',
        'glob',
        'pandas',
        'multiprocessing',
        'os',
    ],
    author='Aayush Saxena',
    author_email='aayush.saxena@physics.ox.ac.uk',
    description='A Python package for obtaining photometry in HST and JWST filters from NIRSpec PRISM spectra',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/aayush3009/phot_from_spec',
    license='MIT',
)