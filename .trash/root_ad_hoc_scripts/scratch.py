import pandas as pd
import warnings
import gc

def generate_data(path, chunk_size):
    reader = pd.read_csv(path, chunksize=chunk_size)
    try:
        for f in reader:
            yield f
    finally:
        reader.close()

# Create dummy
with open('dummy.csv', 'w') as f:
    f.write('A\n')
    for i in range(100):
        f.write(f'{i}\n')

import sys
warnings.simplefilter('always', ResourceWarning)

gen = generate_data('dummy.csv', 10)
for chunk in gen:
    break
    
print("generator stopped")
gen.close() # explicitly close or left for GC
print("generator closed")

