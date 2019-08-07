import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()

setuptools.setup(
    name='chemdescriptor',
    version='0.0.5',
    author="DRP Project",
    author_email="darkreactionproject@haverford.edu",
    description="A standalone module to help generate molecular descriptors from various cheminformatics software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/darkreactions/chemdescriptor",
    packages=['chemdescriptor'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['pandas'],
    include_package_data=True,
    entry_points={
        'console_scripts': ['chemdescriptor-cx=chemdescriptor.command_line:main'],
    }
)
