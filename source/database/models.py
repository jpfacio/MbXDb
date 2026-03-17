from sqlalchemy.orm import declarative_base
from sqlalchemy import Column, Integer, String

Base = declarative_base()

class Genes(Base):
    __tablename__ = "genes"
    Gene_Tag = Column(String, primary_key=True)
    Contig = Column(String)
    Start = Column(Integer)
    Stop = Column(Integer)
    Product = Column(String)
    Databases = Column(String)
    Product_Sequence = Column(String)
    Bin = Column(String)

class Bins(Base):
    __tablename__ = "bins"
    Bin = Column(String, primary_key=True)
    Sample = Column(String)
    Project = Column(String)
    Date = Column(String)
    Study_ID = Column(String)

class Studies(Base):
    __tablename__ = "studies"
    Study_ID = Column(String, primary_key=True)
    Title = Column(String)
    Keywords = Column(String)
    Authors = Column(String)
    Journal = Column(String)
    Year = Column(String)
    Affiliations = Column(String)