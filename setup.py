from setuptools import setup

setup(
    name="cfanalysis",
    version="0.0.1",
    description="Clusters functional analysis",
    url="https://github.com/patryk-jarnot/gbsc_clusters_functional_analysis",
    author="Joanna Ziemska-Legiecka, Aleksandra Gruca, Patryk Jarnot",
    author_email="patryk.jarnot@gmail.com",
    license="MIT",
    packages=["cfanalysis", "cfanalysis.src"],
    install_requires=[
        "statsmodels",
        # "sqlite3",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
    ],
)
