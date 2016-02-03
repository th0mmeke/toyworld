"""
Created on 9/02/2013

@author: thom
"""

import os

DataDir = os.environ.get('TOYWORLDDATA') or os.getcwd()
TestDir = os.path.join(DataDir, 'test')
EnergyTolerance = 2E-2
