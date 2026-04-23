import concurrent.futures
import time

def worker(i):
    def gen():
        for j in range(3):
            time.sleep(i * 0.1)
            yield f"worker {i} chunk {j}"
    return gen()

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    futures = [executor.submit(worker, i) for i in [5, 1]]
    for f in concurrent.futures.as_completed(futures):
        for chunk in f.result():
            print(chunk)
