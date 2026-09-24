import subprocess
from pathlib import Path
import json
import functions as f

#Defining directories and files

tmp = Path("tmp")
out_jsonl = tmp / "gulf.jsonl"

# PRJNA870083

biopr = ["PRJNA870083"]

gcf_info = f.acq.get_gcf_info(biopr[0])
        
# PRJNA340003

srrs = ["SRR4342129", "SRR4342130", "SRR4342133", "SRR4342134", "SRR4342135", "SRR4342136"]

srr_info = f.acq.get_srr_info(srrs[0])










