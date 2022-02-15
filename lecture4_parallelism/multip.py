
import multiprocessing
import time
import os
import sys
import logging
import random


root = logging.getLogger()
root.handlers = []
logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                    format='(%(threadName)-9s-pid %(process)d) %(message)s',)


def worker(num):
    """thread worker function"""
    
    sleep_time = random.randint(1, 8)
    logging.debug(f'worker: {num} sleeping for {sleep_time} s, pid: {os.getpid()}')
    time.sleep(sleep_time)
    logging.debug(f'done: worker: {num}')
    return

if __name__ == "__main__":
    procs = []
    for i in range(2):
        p = multiprocessing.Process(target=worker, args=(i,))
        procs.append(p)
        p.start()
