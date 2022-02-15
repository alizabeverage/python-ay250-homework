def main():
    import argparse
    from calcalc.CalCalc import calculate
    
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
