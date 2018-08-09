import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cathpy",
    version="0.0.1",
    author="Ian Sillitoe",
    author_email="i.sillitoe@ucl.ac.uk",
    description="Protein Bioinformatics related helpers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sillitoe/cathpy",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
