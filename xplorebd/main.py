"""
Main module for Explore_Biological_Databases.
"""

import sys
from .cli import main as cli_main

def main():
    return cli_main()

if __name__ == "__main__":
    sys.exit(main())
