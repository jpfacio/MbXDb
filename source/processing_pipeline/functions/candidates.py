import pandas as pd
import subprocess
import os

def against_pah(query, base_dir):
    
    dbs = []
    
    for i in os.listdir(base_dir):
        if i.endswith('.faa'):
            dbs.append(i[:-4])
            
    for i in dbs:
        
        os.system(f'diamond makedb --in support_files/pah_db_v2.4/{i}.faa --db support_files/pah_db_v2.4/{i}')
        
        genes_load = pd.read_csv('Data/Entities/genes.csv')
        
        sequences = genes_load[['Gene_Tag', 'Product_Sequence']]
        
        with open('Data/Raw/sequences.faa', 'w') as f:
            for idx, row in sequences.iterrows():
                f.write(f'>{row["Gene_Tag"]}\n')
                f.write(f'{row["Product_Sequence"]}\n')
        
        
        subprocess.run(['mkdir', 'tmp/pah_log'])
        
        for i in dbs:
            
            diamond_result = subprocess.run([
                'diamond', 'blastp', '--query', 'candidates.faa', '--db', f'support_files/pah_db_v2.4/{i}.dmnd', 
                '--out', f'tmp/pah_log/results_{i}.txt', '--outfmt', '6'
            ], capture_output=True, text=True)
            
            with open(f'tmp/pah_log/stderr_{i}.log', 'w') as f:
                f.write(diamond_result.stderr)
                
            with open(f"tmp/pah_log/stdout_{i}.log", 'w') as f:
                f.write(diamond_result.stdout)
                
            print(f'{i} Done!')
        
                
        
        
    

