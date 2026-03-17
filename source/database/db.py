import os
import sys
import csv
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Set up paths so this works from mxd.sh
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, "../../"))
if root_dir not in sys.path:
    sys.path.append(root_dir)
    
from source.database.models import Base, Genes, Bins, Studies

# Setup Database
db_path = os.path.join(root_dir, "mxd.db")
engine = create_engine(f"sqlite:///{db_path}")
Session = sessionmaker(bind=engine)

def populate():
    Base.metadata.create_all(engine)
    session = Session()
    
    # Path to your CSVs
    entities_path = os.path.join(root_dir, "Data", "Entities")
    
    files_to_models = {
        "studies.csv": Studies,
        "bins.csv": Bins,
        "genes.csv": Genes
    }

    for file_name, model in files_to_models.items():
        full_path = os.path.join(entities_path, file_name)
        if os.path.exists(full_path):
            with open(full_path, 'r') as f:
                reader = csv.DictReader(f)
                # This takes the CSV row and maps it to your class columns
                session.bulk_insert_mappings(model, list(reader))
            print(f"Done loading {file_name}")
    
    session.commit()
    session.close()

if __name__ == "__main__":
    populate()