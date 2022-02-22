import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="desclamp",
    version="0.0.1",
    author="Remy Joseph",
    author_email="remy.joseph11@gmail.com",
    description="DESC Strong Lensing postage stamp extraction and lensed source injection pipeline.",
    long_description="This section of DESC SL pipeline provides tools to  extract deep coadd images of strong gravitational lens candidates as well as light curves from DC2 and from the Rubin data. The preselection of lens candidates is done at the catalog levels. Postage stamps and light curves can then be used to conduct detailed strong gravitational lens searches, that usually rely on Machine Learning techniques. This pipeline section also provides tools to insert synthetic lensed sources in Rubin (and DC2) images. Synthetic injection will be used to train and characterise lens finders.",
    long_description_content_type="text/markdown",
    url="https://github.com/herjy/DESC-Lamp/",
    project_urls={
        "Bug Tracker": "https://github.com/herjy/DESC-Lamp/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: BSD 3-Clause License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "./"},
    packages=setuptools.find_packages(where="./"),
    python_requires=">=3.6",
)
