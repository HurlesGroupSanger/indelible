from setuptools import setup

setup(
    name="Indelible",
    version="0.9",
    install_requires=[
        "pysam"
    ],
    author = "Alejandro Sifrim",
    author_email = "as33@sanger.ac.uk",
    description = "A structural variant caller for genomic data",
    license = "GPL",
    keywords = "SV genomics variant caller",
    url = "https://github.com/asifrim/indelible"   # project home page, if any
)