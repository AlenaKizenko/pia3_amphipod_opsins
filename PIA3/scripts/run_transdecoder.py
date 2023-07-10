import subprocess
import sys
import os


input_file = sys.argv[1]
wild_query = sys.argv[2]


os.chdir(wild_query)

longorfs = subprocess.run([f'TransDecoder.LongOrfs -t {input_file}'], shell = True)

predict = subprocess.run([f'TransDecoder.Predict -t {input_file} --single_best_only'], shell = True)

if predict.returncode != 0:
    predict_no_refine = subprocess.run(f'TransDecoder.Predict -t {input_file} --single_best_only --no_refine_starts', shell = True)

sed = subprocess.run([f"sed 's/[|,~,/,\,]/_/g' {wild_query}.fasta.transdecoder.pep > {wild_query}.fasta.transdecoder_underscores.pep"], shell = True)