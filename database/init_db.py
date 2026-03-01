import sqlite3
from pathlib import Path

DB_PATH = Path(__file__).parent / "mxd.db"
SCHEMA_PATH = Path(__file__).parent / "schema.sql"

def create_database():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    with open(SCHEMA_PATH, "r") as f:
        cursor.executescript(f.read())
        
    conn.commit()
    conn.close()
    
    print("Database successfully created")
    
if __name__ == "__main__":
    create_database()