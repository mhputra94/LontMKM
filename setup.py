from setuptools import setup, find_packages

setup(
    name="LontMKM",
    version="0.1.0",
    description="Light-Oriented Nano-catalysis Toolkit for MicroKinetic Modeling",
    author="Your Name",
    author_email="you@example.com",
    url="https://github.com/<username>/LontMKM",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "scipy",
        "pyyaml",
        "matplotlib"
    ],
    entry_points={
        "console_scripts": [
            "lontmkm=LontMKM.cli:main",
        ],
    },
    python_requires=">=3.8",
)

