
import requests
from bs4 import BeautifulSoup
from concurrent.futures import ProcessPoolExecutor

url = "https://en.wikipedia.org/wiki/Special:Random"

def run():
    """worker function"""
    
    resp = requests.get(url).text
    title = BeautifulSoup(resp, 'html.parser').title.string.split('- Wikipedia')[0]
    return len(title)


e = ProcessPoolExecutor()
lens = list(e.map(run,range(10)))
print(lens)
