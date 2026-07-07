from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker
from models import Base, Genes
import csv

engine = create_engine('sqlite:///mxd.db')
Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)
session = Session()

print("Populating Genes entity")
batch_size = 1000
batch = []

with open('Data/Entities/genes.csv', 'r', encoding='utf-8') as file:
    csv_reader = csv.DictReader(file)
    
    for i, row in enumerate(csv_reader):
        gene = Genes(Gene_Tag=row['Gene_Tag'],
                     Contig=row['Contig'],
                     Start=row['Start'],
                     Stop=row['Stop'],
                     Product=row['Product'],
                     Databases=row['Databases'],
                     Product_Sequence=row['Product_Sequence'],
                     Bin=row['Bin']
        )
        batch.append(gene)
        
        if len(batch) >= batch_size:
            session.bulk_save_objects(batch)
            session.commit()
            batch = []
            
            
    if batch:
        session.bulk_save_objects(batch)
        session.commit()
        

session.close()

print('Done!')


