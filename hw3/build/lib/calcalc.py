import numexpr as ne
import urllib.request

appID = 'JWKAQ4-U2LKHL8YEE'


def eval_wolfram(s, return_float=True):
    s = s.replace(' ','+')
    f = urllib.request.urlopen(f'http://api.wolframalpha.com/v1/result?appid={appID}&i={s}%3F&units=metric')
    result = str(f.peek(),"utf-8")
    result = result.replace('about ','') # remove the `about`

    # convert result to a float
    if return_float: 
        result = result.replace(' times 10 to the ','e') # put in scientific notation
        # replace common words with numerical value
        sub_words = {' thousand':'e3', ' million':'e6', ' billion':'e8', ' trillion':'e9'}
        for key in sub_words.keys():
            result = result.replace(key, sub_words[key])
        return(float(result.split(' ')[0]))

    # or return the string
    else: return(result)

def calculate(s, run_python=False, run_wolfram=False, return_float=True):
    '''
    use numexpr to evaluate string expression
    return result in desired format
    '''
    if run_python:
        try: 
            return(ne.evaluate(s).item())
        except Exception as e:
            print("Oops!", e.__class__, "occurred.")
            return("Perhaps try running with Wolfram (-w flag)")
    
    elif run_wolfram:
        return eval_wolfram(s, return_float)
    
    # if for some reason both python and wolfram flags are set, try evaluating
    # in python. If test fails, evaluate in wolfram
    else: 
        try: 
            return(ne.evaluate(s).item())
        except:
            return eval_wolfram(s, return_float)

    
if __name__ == '__main__':
    import argparse
    
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Evaluate a string.')
    
    # define input argument
    parser.add_argument('string',
                       help='String to be evaulated')
    
    # -s will try to evaluate using python
    parser.add_argument('-s', action='store_true', default=False,
                        dest='run_python',
                        help='Run in python?')
    
    # -w will try to evaluate using wolfram
    parser.add_argument('-w', action='store_true', default=False,
                        dest='run_wolfram',
                        help='Send to wolfram?')
    # add versioning
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    results = parser.parse_args()
    print(calculate(results.string, run_python=results.run_python, run_wolfram=results.run_wolfram))
    
