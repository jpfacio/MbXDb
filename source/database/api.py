from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from models import Base, Genes

engine = create_engine('sqlite:///../../mxd.db')
SessionLocal = sessionmaker(bind=engine)

app = FastAPI()

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
        
@app.get("/genes/{gene_tag}")
def get_gene(gene_tag: str, db: Session = Depends(get_db)):
    gene = db.query(Genes).filter(Genes.Gene_Tag == gene_tag).first()
    
    if not gene:
        raise HTTPException(status_code=404, detail='Gene not found')
    
    return {
        "Gene_Tag": gene.Gene_Tag,
        "Contig": gene.Contig,
        "Start": gene.Start,
        "Stop": gene.Stop,
        "Product": gene.Product,
        "Databases": gene.Databases,
        "Bin": gene.Bin
    }
    
@app.get("/genes/")
def list_genes(limit: int = 10, skip: int = 0, db: Session = Depends(get_db)):
    genes = db.query(Genes).offset(skip).limit(limit).all()
    return genes

@app.get("/genes/product2gene/{keyword}")
def product2gene(keyword: str, limit: int = 50, skip: int = 0, db: Session = Depends(get_db)):
    query = db.query(Genes).filter(
        Genes.Product.ilike(f"%{keyword}%")
    )
    
    total = query.count()
    genes = query.offset(skip).limit(limit).all()
    
    return {
        'keyword': keyword,
        'total': total,
        'limit': limit,
        'skip': skip,
        'data': genes
    }