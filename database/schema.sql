CREATE TABLE studies (
    study_id TEXT PRIMARY KEY,
    title TEXT,
    keywords TEXT,
    authors TEXT,
    journal TEXT,
    year TEXT,
    affiliations TEXT
);

CREATE TABLE bins (
    bin TEXT PRIMARY KEY, 
    sample TEXT,
    project TEXT,
    date TEXT,
    study_id TEXT,
    FOREIGN KEY (study_id) REFERENCES studies(study_id)
);

CREATE TABLE genes (
    gene_tag TEXT PRIMARY KEY,
    contig TEXT, 
    start INTEGER,
    stop INTEGER,
    product TEXT,
    databases TEXT, 
    product_sequence TEXT,
    bin TEXT,
    FOREIGN KEY (bin) REFERENCES bins(bin)
);