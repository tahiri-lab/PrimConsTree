from setuptools import find_packages, setup

setup(
    name="primconstree",
    version="0.1.0",
    description="PrimConsTree consensus tree algorithm",
    author="Elio Torquet",
    packages=find_packages(),
    install_requires=["ete3==3.1.3", "networkx==3.3", "matplotlib==3.9.0"],
)
