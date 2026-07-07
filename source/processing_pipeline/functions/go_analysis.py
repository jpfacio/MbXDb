from goatools.obo_parser import GODag
import pandas as pd

go_basic = "support_files/go-basic.obo"
godag = GODag(go_basic)

def analyze_go_terms(genes):
    genes = genes.dropna(subset=['go']).copy()
    
    def to_list(x):
        if isinstance(x, str):
            return eval(x)
        return x if isinstance(x, list) else []
    
    genes['go'] = genes['go'].apply(to_list)
    genes = genes[genes['go'].apply(lambda x: len(x) > 0)]
    
    rows = []
    
    for idx, row in genes.iterrows():
        for go_term in row['go']:
            if go_term in godag:
                term = godag[go_term]
                rows.append({
                    'Gene_Tag': row['Gene_Tag'],
                    'go_term': go_term,
                    'name': term.name,
                    'namespace': term.namespace,
                    'level': term.level,
                    'depth': term.depth
                })
    
    return pd.DataFrame(rows)

def get_hierarchy(go_df):
    
    gene_grouped = go_df.groupby('Gene_Tag').agg(list).reset_index()
    
    rows_to_keep = []
    
    for idx, row in gene_grouped.iterrows():
        namespaces = row['namespace']
        levels = row['level']
        names = row['name']
        go = row['go_term']
        
        best_index = {}
        
        for i, ns in enumerate(namespaces):
            if ns not in best_index:
                best_index[ns] = i
            else:
                current_best_level = levels[best_index[ns]]
                current_level = levels[i]
                
                if current_level > current_best_level:
                    best_index[ns] = i
                elif current_level == current_best_level:
                    pass
                
        selected_index = list(best_index.values())
        
        for idx_sel in selected_index:
            rows_to_keep.append({
                'Gene_Tag': row['Gene_Tag'],
                'go_term': go[idx_sel],
                'name': names[idx_sel],
                'namespace': namespaces[idx_sel],
                'level': levels[idx_sel],
                'source': 'original'
            })
    
    go_info_clean = pd.DataFrame(rows_to_keep)
    
    all_rows = go_info_clean.to_dict('records')
    
    for idx, row in go_info_clean.iterrows():
        go_id = row['go_term']
        gene_tag = row['Gene_Tag']
    
        
        if go_id not in godag:
            continue
            
        term = godag[go_id]
        
        for ancestor_id in term.get_all_parents():
                ancestor = godag[ancestor_id]
                if ancestor.level > 0:
                    all_rows.append({
                        'Gene_Tag': gene_tag,
                        'go_term': ancestor.id,
                        'name': ancestor.name,
                        'namespace': ancestor.namespace,
                        'level': ancestor.level,
                        'source': 'derivative'
                    })
    
    go_info_full =  pd.DataFrame(all_rows)
    go_info_full = go_info_full.sort_values(['Gene_Tag', 'level']).reset_index(drop=True)
    
    return go_info_full
    
    
                    
            
            
            
        
        
    
                    
                
        
        
    
    
    



















                    
            