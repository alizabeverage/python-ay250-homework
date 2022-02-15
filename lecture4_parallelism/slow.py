from time import sleep, time

from concurrent.futures import ProcessPoolExecutor
e = ProcessPoolExecutor() 

def slowfunc(x):
    sleep(1)
    return(x+1)

if __name__ == "__main__":
    s = time()
    results = list(e.map(slowfunc, range(8)))
    print(f"Finished in {time() - s:0.3f} sec")
    print(results)
    
    e.shutdown()
