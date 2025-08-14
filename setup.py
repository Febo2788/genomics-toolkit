from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="genomics-toolkit",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive bioinformatics pipeline for DNA/RNA sequence analysis and variant calling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/genomics-toolkit",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "genomics-toolkit=genomics_toolkit.cli:main",
            "genomics-toolkit-gui=genomics_toolkit.gui:main",
        ],
    },
    include_package_data=True,
    package_data={
        "genomics_toolkit": ["data/*", "templates/*"],
    },
)