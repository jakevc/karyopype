from pathlib import Path
from setuptools import setup
from setuptools import find_packages

__version__ = "0.1.0"

HERE = Path(__file__).parent

README = (HERE / "README.md").read_text()

install_requires = [
    "matplotlib", "pandas", "numpy",
    ]

setup(
    name="karyopype",
    packages=find_packages(exclude=("tests",)),
    package_dir={__name__: __name__},
    package_data={
        'karyopype': ['data/chromsizes/*.chrom.sizes']
    },
    include_dirs=["."],
    # include_package_data=True,
    version=__version__,
    description="Chromosomal visualization in Python.",
    long_description=README,
    long_description_content_type='text/markdown',
    author="Jake VanCampen",
    author_email="jake.vancampen7@gmail.com",
    url="http://github.com/jakevc/karyopype",
    keywords=["Bioinformatics"],
    license="MIT",
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta", "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        'License :: OSI Approved :: MIT License',
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
)