"""
pytest configuration file
"""

import sys
import os

# Add FIDASIM lib/python to Python path for all tests
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../lib/python'))
