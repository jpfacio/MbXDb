from sqlalchemy.orm import declarative_base
from sqlalchemy import Column, Integer, String

Base = declarative_base()

class Genes(Base):
    __tablename__ = 'genes'
    
    Gene_Tag = Column(String, primary_key=True)
    Contig = Column(String)
    Start = Column(Integer)
    Stop = Column(Integer)
    Product = Column(String)
    Databases = Column(String)
    Product_Sequence = Column(String)
    Bin = Column(String)