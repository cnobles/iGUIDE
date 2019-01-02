from distutils.core import setup
from setuptools import find_packages

setup(
    name="iguide",
    use_scm_version = True,
    setup_requires=['setuptools_scm'],
    #packages=["iguidelib"],
    packages=find_packages(),
    include_package_data=True,
    package_data={"iguidelib": ["iguidelib/data/*.yml"]},
    entry_points={'console_scripts': [
        'iguide = iguidelib.scripts.command:main',
    ]},
    classifiers=[
        'Programming Language :: Python :: 3.4'
    ]
)
