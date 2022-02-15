from time import sleep

def slowfunc(x, y, delay=1):
    sleep(delay)
    return(x + y)
