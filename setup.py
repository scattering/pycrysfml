from setuptools import setup, find_packages

packages = find_packages()
print("packages", packages)

setup(
    name='pycrysfml', 
    version='0.0.1', 
    packages=packages, 
    install_requires=[],
)
