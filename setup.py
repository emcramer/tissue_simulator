from setuptools import setup, find_packages

setup(
    name="tissue_simulator",
    version="0.1.2",
    description="3D simulated biological tissue section generator",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        "PyQt5>=5.15.0",
        "scipy>=1.7.0",
        "networkx>=2.6.0",
        "pandas>=1.3.0",
        "statsmodels>=0.13.0",
    ],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    entry_points={
        "console_scripts": [
            "tissue-simulator=tissue_simulator.gui:main",
        ],
    },
)
